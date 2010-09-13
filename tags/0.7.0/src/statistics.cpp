/** statistics.cpp
 *  29 Oct 2009--02 Jun 2010
 *  Neal Davis
 *  
 */

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numeric> //FIXME
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

extern ParticleMap      pmap;

extern list<Particle*>  surfaceA;
extern list<Particle*>  surfaceB;

extern uint *oldParticleCount;
extern uint *newParticleCount;
extern uint *myParticleCount;
extern Vector   rates;
extern uint     initialTotalParticles;

//  Global variables
static uint oldTotalAtoms, newTotalAtoms;

/*  bool fileExists()
 *  
 *  Returns whether or not a file exists.  From http://www.techbytes.ca/techbyte103.html
 */
bool fileExists(char* fileName)
{   //  Attempt to get the file attributes.
    struct stat stFileInfo; 
    int intStat = stat(fileName, &stFileInfo); 
    if(intStat == 0)
    //  The file exists.
    {   return true; }
    else
    //  The file does not exist, or we do not have permission to access it.
    {   return false; } }

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
void outputToFile(uint t, uint totalAtoms)
{   fstream outDataFile;
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
        uint numNN;
        for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
        {   if (iter->second.state >= dissolnStates) continue;
            numNN = accumulate(iter->second.countN.begin(), iter->second.countN.begin() + dissolnStates, 0);
            outDataFile << iter->second.state << "\t" << iter->second.x << "\t"
                        << iter->second.y     << "\t" << iter->second.z << endl; }
        outDataFile.close(); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to write output file %s.", outDataFileName);
        error(err); }
    
    static fstream  statFile;
    static char     statFileName[36];
    if (!rank)
    {   sprintf(statFileName, "N=%d.dat", initialTotalParticles);
        
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
        {   if (t == 0)
            {   statFile << "random seed," << ranseed << endl; }
            statFile << t << ",";
            for (Vector::iterator iter = rates.begin(); iter != rates.end(); iter++)
            {   statFile << *iter << ","; }
            statFile << surfaceA.size() << "," << surfaceB.size() << endl;
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
{   //  Get number of each particles in each state across this and all processes.
    for (uint i = 0; i < numStates; i++)
    {   if (rank == 0) oldParticleCount[i] = newParticleCount[i];
        myParticleCount[i] = 0; }
    
    //  Count the total of particles in each state.
    for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
    {   myParticleCount[iter->second.state]++; }
    MPI_Reduce(myParticleCount, newParticleCount, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //  Get the total number of (non-dissolved) atoms in the system.
    oldTotalAtoms = newTotalAtoms;
    newTotalAtoms = 0;
    for (uint i = 0; i < dissolnStates; i++)
    {   newTotalAtoms += newParticleCount[i]; }
    
    //  Calculate rates if data exist.
    for (uint i = 0; i < numStates; i++)
    {   rates[i] = (double) newParticleCount[i] - (double) oldParticleCount[i]; }
    //FIXME:why rates seem to be wrong or zero always
    
    //  Write the system configuration and statistics to disk.
    outputToFile(t, newTotalAtoms); }

/*  collectStatFiles()
 *  
 *  Load all files for each time step and put them together into one.
 */
void collectStatFiles()
{   fstream allFile;
    char    allFileName[36];
    sprintf(allFileName, "N=%d.xyz", initialTotalParticles);
    allFile.open(allFileName, ofstream::out);
    if (!allFile)
    {   char err[64];
        sprintf(err, "Unable to create file %s", allFileName);
        error(err); }
    
    try
    {   uint    totalParticles = 0;
        state_t state;
        coord_t x, y, z;
        fstream inFile;
        char    inFileName[36];
        for (uint t = 0; t < tmax; t += outputInterval)
        {   for (uint i = 0; i < (uint) size; i++)
            {   //  Open each file, read its contents into allfile, and close it.
                sprintf(inFileName, "t=%d-P=%d.xyz_temp", t, i);
                inFile.open(inFileName, ifstream::in);
                if (!inFile)
                {   char err[64];
                    sprintf(err, "Unable to load file %s", inFileName);
                    error(err); }
                
                //  Copy total number of particles in this time step.
                if (!i)
                {   inFile >> totalParticles;
                    allFile << totalParticles << endl; }
                
                while (!inFile.eof())
                {   inFile >> state >> x >> y >> z;
                    allFile << state << "\t" << x << "\t" << y << "\t" << z << endl; }
                
                inFile.close();
                remove(inFileName); } }
        
        allFile.close(); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to collate position data to file %s.", allFileName);
        error(err); } }
