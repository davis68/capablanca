/** statistics.cpp
 *  29 Oct 2009--21 Oct 2010
 *  Neal Davis
 *  
 */

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <vector>

#include "definitions.h"
#include "error.h"
#include "particle.h"
#include "process.h"
#include "statistics.h"

using namespace std;

//  External global variables from global.cpp
extern bool verbose;

extern uint numStates,
            dissolnStates,
            maxNN;

extern double   T,
                phi;

extern uint outputInterval,
            tmax,
            ranseed;

extern int  rank,
            size;

extern char *progCL;

extern ParticleMap      pmap;

extern list<Particle*>  surfaceA,
                        surfaceB;

extern uint *oldParticleCount,
            *newParticleCount,
            *myParticleCount;
extern Vector   rates;
extern uint     initialTotalParticles;

//  Global variables
static uint totalAtoms,
            dissolnAtoms;

/*  bool fileExists()
 *  
 *  Returns whether or not a file exists.
 *  From http://www.techbytes.ca/techbyte103.html
 */
bool fileExists(char* fileName)
{   //  Attempt to get the file attributes.
    struct stat stFileInfo; 
    int intStat = stat(fileName, &stFileInfo); 
    
    if (!intStat)   return true;    //  The file exists.
    else            return false; } //  The file does not exist, or we do not have permission to access it.

/** outputSurface()
 *  
 *  Output a particle data to a file.
 */
void outputSurface(list<Particle*> particles, uint t)
{   fstream outDataFile;
    char    outDataFileName[36];
    sprintf(outDataFileName, "surface-t=%d-P=%d.xyz", t, rank);
    outDataFile.open(outDataFileName, fstream::out);
    if (!outDataFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", outDataFileName);
        error(err); }
    
    try
    {   //  Output particle positions to outDataFile.
        outDataFile << particles.size() << endl;
        for (list<Particle*>::iterator iter = particles.begin(); iter != particles.end(); iter++)
        {   outDataFile << (*iter)->state << "\t" << (*iter)->x << "\t"
                        << (*iter)->y     << "\t" << (*iter)->z << endl; }
        
        outDataFile.close(); }
    
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to write output file %s.", outDataFileName);
        error(err); } }

/** outputToFile()
 *  
 *  Output all data to a file by time step and process.
 */
void outputToFile(uint t)
{   /// Output particle positions and states to file.
    fstream outDataFile;
    char    outDataFileName[36];
    sprintf(outDataFileName, "t=%d-P=%d.xyz_temp", t, rank);
    outDataFile.open(outDataFileName, fstream::out);
    if (!outDataFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", outDataFileName);
        error(err); }
    
    try
    {   //  Output non-dissolved particle positions to files which will later be collected.
        if (!rank) outDataFile << totalAtoms << endl;
        for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
        {   if (iter->second.state >= dissolnStates) continue;
            outDataFile << iter->second.state << "\t" << iter->second.x << "\t"
                        << iter->second.y     << "\t" << iter->second.z << endl; }
        
        outDataFile.close(); }
    
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to write output file %s.", outDataFileName);
        error(err); }
    
    /// Output system statistics (state counts, etc.) to file.
    static fstream  statFile;
    static char     statFileName[36];
    if (!rank)
    {   sprintf(statFileName, "N=%d.csv", initialTotalParticles);
        
        //  Check to make sure the file doesn't already exist, as we'll append to it.
        if (t == 0 && fileExists(statFileName))
        {   try
            {   remove(statFileName); }
            
            catch(...)
            {   char err[64];
                sprintf(err, "⚠ Unable to delete statistics file %s from previous program run.", statFileName);
                warning(err); } }
        
        statFile.open(statFileName, fstream::out | fstream::app);
        if (!statFile)
        {   char err[64];
            sprintf(err, "⚠ Unable to load file %s", statFileName);
            error(err); }
        
        try
        {   if (!t) statFile << "random seed," << ranseed << endl;
            
            statFile << t << ",rates,";
            for (Vector::iterator iter = rates.begin(); iter != rates.end(); iter++)
            {   statFile << *iter << ","; }
            statFile << "surfaces," << surfaceA.size() << "," << surfaceB.size() << endl;
            
            statFile.close(); }
        
        catch (...)
        {   char err[64];
            sprintf(err, "⚠ Unable to write statistics file %s.", statFileName);
            error(err); } } }

/** collateStatistics()
 *  
 *  Output the system configuration and empirical calculations to disk.
 */
void collateStatistics(uint t)
{   //  Clear the variables.
    for (uint i = 0; i < numStates; i++)
    {   if (!rank)  oldParticleCount[i] = newParticleCount[i];
        myParticleCount[i] = 0; }
    
    //  Count the total of particles in each state.
    for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
    {   myParticleCount[iter->second.state]++; }
    MPI_Reduce(myParticleCount, newParticleCount, numStates, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //  Get the total number of atoms in the system.
    totalAtoms   = 0;
    dissolnAtoms = 0;
    
    for (uint i = 0; i < dissolnStates; i++)
    {   totalAtoms += newParticleCount[i]; }
    for (uint i = 0; i < numStates; i++)
    {   dissolnAtoms += newParticleCount[i]; }
    
    //  Calculate rates if data exist.
    for (uint i = 0; i < numStates; i++)
    {   rates[i] = (double) newParticleCount[i] - (double) oldParticleCount[i]; }
    
    //  Write the system configuration and statistics to disk.
    outputToFile(t); }

/*  collectStatFiles()
 *  
 *  Load all files for each time step and combine them together into one file
 *  for each time step.
 */
void collectStatFiles()
{   try
    {   uint    totalParticles = 0;
        state_t state;
        coord_t x, y, z;
        fstream allFile, inFile;
        char    allFileName[36], inFileName[36];
        
        for (uint t = 0; t < tmax; t += outputInterval)
        {   for (uint i = 0; i < (uint) size; i++)
            {   //  Open each input file for the time step.
                sprintf(inFileName, "t=%d-P=%d.xyz_temp", t, i);
                inFile.open(inFileName, ifstream::in);
                if (!inFile)
                {   char err[64];
                    sprintf(err, "Unable to load file %s", inFileName);
                    error(err); }
                
                //  Open output file for the time step.
                sprintf(allFileName, "N=%d-t=%d.xyz", initialTotalParticles, t);
                allFile.open(allFileName, ofstream::out);
                if (!allFile)
                {   char err[64];
                    sprintf(err, "Unable to create file %s", allFileName);
                    error(err); }
                
                //  Copy total number of particles in this time step.
                if (!i)
                {   inFile >> totalParticles;
                    allFile << "# " << progCL << endl;
                    allFile << totalParticles << endl; }
                
                //  Read the contents into allFile.
                while (!inFile.eof())
                {   inFile >> state >> x >> y >> z;
                    allFile << state << "\t" << x << "\t" << y << "\t" << z << endl; }
                
                //  Clean up.
                inFile.close();
                allFile.close();
                remove(inFileName); } } }
    
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to collate position data to file.");
        error(err); } }

