//  TODO:  Output total configuration, run time, seeds, etc.

#include <fstream>
#include <iostream> //DEBUG
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "mpi.h"
#include "particle.h"
#include "process.h"
#include "process_helpers.h"
#include "statistics.h"

using namespace std;

extern bool     verbose;

extern uint numStates,
            dissolnStates,
            maxNN;

extern double   T,
                phi;

extern uint outputInterval,
            tmax;

extern int  rank,
            size;

extern ParticleMap          pmap;

extern vector<uint>     oldParticleCount;
extern vector<uint>     myParticleCount;
extern vector<double>   rates;
extern uint             initialTotalParticles;

/*  outputToFile(int)
 *  
 *  Output all data to a file by time step and process.
 */
void outputToFile(int t)
{   //  Get the total number of (non-dissolved) atoms in the system.
    uint    myAtoms = 0,
            totalAtoms = 0;
    for (uint i = 0; i < numStates; i++)
    {   myAtoms += myParticleCount[i]; }
    MPI_Reduce(&myAtoms, &totalAtoms, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
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
        for (pmapiter iter = pmap.begin(); iter != pmap.end(); iter++)
        {   if (iter->second.state < dissolnStates)
            {   outDataFile << iter->second.state << "\t" << iter->second.x << "\t"
                            << iter->second.y     << "\t" << iter->second.z << endl; } }
        outDataFile.close(); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to write output file %s.", outDataFileName);
        error(err); }
    
    if (!rank)
    {   fstream statFile;
        char    statFileName[36];
        sprintf(statFileName, "t=%d-P=%d.dat_temp", t, rank);
        statFile.open(statFileName, fstream::out);
        if (!statFile)
        {   char err[64];
            sprintf(err, "⚠ Unable to load file %s", statFileName);
            error(err); }
        
        try
        {   for (Vector::iterator iter = rates.begin(); iter != rates.end(); iter++)
            {   statFile << *iter << ","; }
                statFile.close(); }
        catch (...)
        {   char err[64];
            sprintf(err, "⚠ Unable to write statistics file %s.", statFileName);
            error(err); } } }

/*  collateStatistics(int)
 *  
 *  Output the system configuration and empirical calculations to disk.
 */
void collateStatistics(int t)
{   //  Calculate rates if data exist.
    for (uint i = 0; i < numStates; i++)
    {   rates[i] = (double) myParticleCount[i] - (double) oldParticleCount[i]; }
    
    //  Get totals of particles stored locally.
    myParticleCount.assign(numStates, 0);
    for(pmapiter i = pmap.begin(); i != pmap.end(); i++)
    {   myParticleCount[i->second.state]++; }
    
    //  Output to a file the system configuration.
    outputToFile(t);
    
    //  Store myParticleCount in oldCount.
    oldParticleCount = myParticleCount; }

/*  collectStateFiles()
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
        sprintf(err, "Unable to load file %s", allFileName);
        error(err); }
    
    try
    {   uint    totalParticles = 0;
        state_t state;
        coord_t x, y, z;
        fstream inFile;
        char    inFileName[36];
        for (uint t = 0; t < tmax; t += outputInterval)
        {   for (uint i = 0; i < size; i++)
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

