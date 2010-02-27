#include <iostream>
#include <cstdlib>
#include "error.h"
#include "mpi.h"

using namespace std;

void error(const char* msg)
{   cerr << msg << endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }

void warning(const char* msg)
{   cerr << msg << endl;
    cerr.flush(); }
