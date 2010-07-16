/** main.cpp
 *  29 Oct 2009--07 Jun 2010
 *  Neal Davis
 *  
 *  capablanca Transition Machine
 *  Neal Davis with contributions by Minas Charalambides, Bill Tuohy, Simon
 *      Jenkins, and Jinquan Zhuang
 *  
 *  0.1 (ppproj2009) Initial codebase; calculates hard-coded Fe(0)-->Fe(VI)
 *      transitions; output broken.
 *  0.2 Added CLI parameters.
 *  0.5 Rewrite of basic functions for expanded rule interactions.
 *  0.6 First completely working version.
 *  0.7 Dissolved (deposition) surface and rules added.
 */

#include <iostream>
#include <mpi.h>
#include <cstdlib>

#include "definitions.h"
#include "error.h"
#include "input.h"
#include "neighbors.h"
#include "particle.h"
#include "process.h"
#include "readxyz.h"
#include "statistics.h"

using namespace std;

//  External global variables from global.cpp
extern int  rank,
            size;

extern Reaction    *rxnA,
                   *rxnB;

extern char *progName,
            *progVers,
            *confFileName,
            *ruleFileName,
            *dataFileName;

int main(int argc, char** argv)
{   //  Perform MPI-required initialization.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    //  Load arguments and initialize parameters.
    double progTime = -MPI_Wtime();
    parseInput(argc, argv);
    loadConfig();
    loadRules();
    
    //  Read in data for each processor.
    double readTime = -MPI_Wtime();
    readXYZ();
    MPI_Barrier(MPI_COMM_WORLD);
    readTime += MPI_Wtime();
    if (!rank)  cout << "▦ Particles read in " << readTime << " s.\n";
    
    //  Calculate nearest neighbors of each particle.
    double nnTime = -MPI_Wtime();
    calculateNeighbors();
    MPI_Barrier(MPI_COMM_WORLD);
    nnTime += MPI_Wtime();
    if (!rank)  cout << "☍ Nearest neighbors calculated in " << nnTime << " s.\n";
    
    //  Perform system calculation.
    double procTime = -MPI_Wtime();
    process();
    MPI_Barrier(MPI_COMM_WORLD);
    procTime += MPI_Wtime();
    if (!rank)  cout << "✻ Processing complete in " << procTime << " s.\n";
    
    //  Collate data.
    if (!rank)
    {   double collTime = -MPI_Wtime();
        FIXME:collectStatFiles();
        collTime += MPI_Wtime();
        cout << "✇ Data output collated and finalized in " << collTime << " s.\n"; }
    
    //  Finalize, free allocated memory and return successfully.
    progTime += MPI_Wtime();
    if (!rank)  cout << "★ Program complete in " << progTime << " s.\n";
    
    delete[] rxnA;
    delete[] rxnB;
    
    cout.flush();
    MPI_Finalize();
    return EXIT_SUCCESS; }

