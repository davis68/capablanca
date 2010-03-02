/*  CAPABLANCA Transition Machine
 *  Oct 09--
 *  0.1 (ppproj2009) Initial codebase; calculates hard-coded Fe(0)-->Fe(VI)
 *      transitions; output broken.
 *  0.2 (incomplete) Added CLI parameters.
 *  0.5 (incomplete) Rewrite of basic functions for expanded rule interactions.
 *  TODO:  blocking into conf file?
 *  TODO:  currently the assertion that all particles are loaded fails sometimes;
 *         the nn calc needs to be validated
 */

#include <cstdio> //***
#include <fstream> //***
#include <iostream>
#include "mpi.h"
#include <cstdlib>
#include "definitions.h"
#include "error.h"
#include "input.h"
#include "neighbors.h"
#include "particle.h"
#include "process.h"
#include "readxyz.h"

using namespace std;

extern int  rank,
            size;

extern ParticleList    particles;  //  One vector of particles per node.
extern uint numStates,      //  Init in input.cpp loadRules()
            dissolnStates,  //  Init in input.cpp loadRules()
            maxNN;          //  Init in input.cpp loadRules()

void dumpParticles()
{   fstream statfile;
    char    statfilename[36];
    sprintf(statfilename, "%d.txt", rank);
    statfile.open(statfilename, fstream::out);
    if (!statfile)
    {   cerr << "⚠ Unable to create output file " << statfilename << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
    
    statfile << particles.size() << endl;
    for (int i = 0; i < particles.size(); i++)
    {   statfile << particles[i].id << "\t"
                 << particles[i].x  << "\t"
                 << particles[i].y  << "\t"
                 << particles[i].z  << endl; }
    statfile.close(); }

int main(int argc, char** argv)
{   //  Perform MPI-required initialization.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    //  Load arguments and initialize parameters.
    parseInput(argc, argv);
    loadConfig();
    loadRules();
    
    //  Read in data for each processor.
    double readTime = -MPI_Wtime();
    readXYZ();
    MPI_Barrier(MPI_COMM_WORLD);
    readTime += MPI_Wtime();
    if (!rank)  cout << "◔ Particles read in " << readTime << " s.\n";
    
    //  Calculate nearest neighbors of each particle.
    double nnTime = -MPI_Wtime();
    calculateNeighbors();
    MPI_Barrier(MPI_COMM_WORLD);
    nnTime += MPI_Wtime();
    if (!rank)  cout << "◑ Nearest neighbors calculated in " << nnTime << " s.\n";
    
    //  Perform system calculation.
    double procTime = -MPI_Wtime();
    process();
    MPI_Barrier(MPI_COMM_WORLD);
    procTime += MPI_Wtime();
    if (!rank)  cout << "◕ Processing complete in " << procTime << " s.\n";
    
    //  Finalize and exit.
    MPI_Finalize();
    return EXIT_SUCCESS; }

