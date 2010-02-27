//  TODO:  Output total configuration, run time, seeds, etc.

#include <fstream>
#include <iostream> //DEBUG
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "mpi.h"
#include "particle.h"
#include "process_helpers.h"
#include "statistics.h"

using namespace std;

extern uint numStates;
extern uint outputInterval,
            tmax;

extern int  rank,
            size;

extern ParticleList particles;
extern ParticleMap  pmap[NUM_OF_THREADS];

//  Output & statistics variables.***
extern uint   *oldCount;           //  Total particles of last step.
extern uint   *myParticleCount;    //  Total particles this step.
extern double *rates;              //  Rate at which each reaction is proceeding (empirical).
extern uint    initialTotalParticles;

//  Output all data to a file by time step and processor.
void outputToFile(int iter, double coverage)
{   //  Get the total number of (non-dissolved) atoms in the system.
    uint    myAtoms = 0,
            totalAtoms = 0;
    for (uint i = 0; i < numStates; i++)
    {   myAtoms += myParticleCount[i]; }
    MPI_Reduce(&myAtoms, &totalAtoms, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //  Output non-dissolved particle positions to files which will later be collected.
    fstream datafile;
    char    datafilename[36];
    sprintf(datafilename, "t=%d-P=%d.xyz_temp", iter, rank);
    datafile.open(datafilename, fstream::out);
    if (!datafile)
    {   cerr << "Unable to create output file " << datafilename << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
    if (!rank) datafile << totalAtoms << endl;
    for (uint i = 0; i < particles.size(); i++)
    {   if (particles[i].state != 6)
        {   datafile << particles[i].state << "\t" << particles[i].x << "\t"
                     << particles[i].y << "\t" << particles[i].z << endl; } }
    datafile.close();
    
    if (!rank)
    {   fstream statfile;
        char    statfilename[36];
        sprintf(statfilename, "t=%d-P=%d.txt", iter, rank);
        statfile.open(statfilename, fstream::out);
        if (!statfile)
        {   cerr << "Unable to create output file " << statfilename << "\n";
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
        
        for (int i = 0; i < numStates; i++)
        {   statfile << rates[i] << ","; }
        statfile << "\n" << coverage << endl;
        statfile.close(); } }

void collateStatistics(int iter)
{   //  Calculate rates if data exist.
    for (uint i = 0; i < numStates; i++)
    {   rates[i] = (double) myParticleCount[i] - (double) oldCount[i];
        #ifdef STAT_DEBUG
        cerr << "rates[" << i << "]=" << rates[i] << ";";
        #endif
    }
    
    //  Calculate the surface coverage of Fe(III).
    //  Sum up all surface particles, and then divide the count of Fe(III) by that.
    uint     mySurfaceCount = 0;
    for(uint i = 0; i < NUM_OF_THREADS; i++)
     mySurfaceCount += surface[i].size();
    uint     totalSurfaceCount;
    MPI_Reduce(&mySurfaceCount, &totalSurfaceCount, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    
    for(uint i = 0; i < numStates; i++) myParticleCount[i] = 0;
    
    #ifdef STAT_DEBUG
    cerr << pmap[0].size() << "::" << endl;
    #endif
    
    for(uint t = 0; t < NUM_OF_THREADS; t++)
    {   for(pmapiter i = pmap[t].begin(); i != pmap[t].end(); i++)
        {   myParticleCount[i->second.state]++;
            #ifdef STAT_DEBUG
            cerr << i->second.state << " | ";
            #endif
        }
    }
    
    #ifdef STAT_DEBUG
    cerr << rank << ":  " << myParticleCount[0] << "; " << myParticleCount[1] << "; " << myParticleCount[2] << "; " << myParticleCount[3] << "; " << myParticleCount[4] << "; " << myParticleCount[5] << "; " << myParticleCount[6] << endl;
    #endif
    
    uint totalParticleCount[numStates];
    MPI_Reduce(myParticleCount, totalParticleCount, numStates, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    uint     myFeIIICount = myParticleCount[3];
    uint     totalFeIIICount;
    MPI_Reduce(&myFeIIICount, &totalFeIIICount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    double  coverage;
    coverage = (double) totalFeIIICount / (double) totalSurfaceCount;
    #ifdef STAT_DEBUG
    cerr << "rank=" << rank << "; " << totalFeIIICount << ", " << totalSurfaceCount;
    cerr << "; coverage=" << coverage << endl;
    #endif

    //  Output to a file the system configuration.
    outputToFile(iter, coverage);

    //  Store myParticleCount in oldCount.
    for (uint i = 0; i < numStates; i++)
    {   oldCount[i] = myParticleCount[i]; } }

//  Load all files for each time step and put them together into one.
void collectStatFiles()
{   fstream allFile;
    char    allFileName[36];
    sprintf(allFileName, "N=%d.xyz", initialTotalParticles);
    allFile.open(allFileName, ofstream::out);
    if (!allFile)
    {   char err[64];
        sprintf(err, "Unable to load file %s", allFileName);
        error(err); }
    
    uint    totalParticles;
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
            remove(inFileName);
        }
    }
    allFile.close();
}
