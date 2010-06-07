/** error.cpp
 *  27 Feb 2010--02 Jun 2010
 *  Neal Davis
 */

#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include "error.h"

using namespace std;

void error(const char* msg)
{   cerr << msg << endl; cerr.flush();
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }

void warning(const char* msg)
{   cerr << msg << endl; cerr.flush();
    cerr.flush(); }
