/*
    The chief problem right now is that only when a particle dissolves is its state
    abroad updated.  We need a plan to update the information every time a particle
    changes at all instead.  States is not currently transferred.

    process_helpers.h uses MPI buffered send mode, and so is incompatible with other
    functions which may also use that mode.
*/
#ifndef _PROCESS_HELPERS_H
#define _PROCESS_HELPERS_H

#include <mpi.h>
#include <vector>
#include <map>
#include "particle.h"
#include "process.h"
#include <string.h>
#include <cstdlib>
#include <iostream> // debug
#include <cassert>
#include "definitions.h"
#include "error.h"

using namespace std;

extern bool verbose;

extern uint dissolnStates;

extern int  rank,
            size;

extern ParticleMap          pmap;
extern vector<ParticlePtr>  surface;

static id_t     interNodeChangesBuffer[MAX_ARRAY_SIZE]; //*** best approach?
static uint     interNodeChangesBufferSize;
static States   interNodeStatesBuffer[MAX_ARRAY_SIZE];
static uint     interNodeStatesBufferSize;
static id_t     interNodeSendBuffer[MAX_ARRAY_SIZE];

static vector<ParticlePtr>  internalChangesBuffer;
static vector<States>       internalStatesBuffer;

static coord_t      maxCoord;
static coord_t      globalMaxCoord;

static uint        *numToRecv;
static id_t         receiveBuffer[MAX_ARRAY_SIZE];
static uint         receiveCount;
static MPI_Status   status;

/*  findSurface()
 *  
 *  Determine the particles whose z-coordinate is less than SURFACE_CUTOFF from
 *  the global maximum z-coordinate.
 */
inline void findSurface()
{   //  Determine local (and thence global) maximum z-coordinates.
    maxCoord = -1e8;
    for (pmapiter iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if (iter->second.z > maxCoord)   maxCoord = iter->second.z; }
    MPI_Allreduce(&maxCoord, &globalMaxCoord, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Add particles within SURFACE_CUTOFF of the top to the surface.
    for (pmapiter iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if (globalMaxCoord - iter->second.z <= SURFACE_CUTOFF) surface.push_back(&(iter->second)); } }

/*  updateOnBoundary(ParticlePtr ptr)
 *  
 *  Add the particle indicated by ptr to the interNodeChangesBuffer.
 */
inline void updateOnBoundary(ParticlePtr ptr)
{   uint initialNCount = accumulate(ptr->initialN.begin(),
                                    ptr->initialN.begin() + dissolnStates, 0);
    
    memcpy(&(interNodeChangesBuffer[interNodeChangesBufferSize]),
           ptr->getNList(), initialNCount * sizeof(id_t));
    interNodeChangesBufferSize += initialNCount;
    /***memcpy(&(interNodeStatesBuffer[interNodeStatesBufferSize]),
           ptr->getNList(), initialNCount * sizeof(id_t));
    interNodeStatesBufferSize += initialNCount;*/ }

/*  updateOffBoundary(ParticlePtr ptr)
 *
 *  Add the particle indicated by ptr to the internalChangesBuffer.
 */
inline void updateOffBoundary(ParticlePtr ptr)
{   internalChangesBuffer.push_back(ptr);
    States st;
    st.oldState = ptr->state;
    st.newState = ptr->state; //*** based on rule transition?
    internalStatesBuffer.push_back(st); }

/*  recv(const int src)
 *
 *  Receive all changes from src.
 */
inline void recv(const int src)
{   if (src == MPI_PROC_NULL) return;
    
    MPI_Recv(receiveBuffer + receiveCount, numToRecv[src], MPI_UNSIGNED,
             src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    //  Retrieve the number of particles sent from the message tag.
    receiveCount += status.MPI_TAG; }

/*  send(const int dest)
 *  
 *  Send all changes to dest.  (Note that capablanca opts to send all data to
 *  all neighboring processes in order to not have to carry around extra data in
 *  the Particle class.)***
 */
inline void send(const int dest)
{   if (dest == MPI_PROC_NULL) return;
    
    MPI_Ssend(interNodeChangesBuffer, numToRecv[rank], MPI_UNSIGNED,
              dest, numToRecv[rank], MPI_COMM_WORLD); }

/*  exchangeInterNodeChanges()
 *  
 *  Initiate the exchange of inter-node data by even processors first, then odd
 *  ones.
 */
inline void exchangeInterNodeChanges()
{   int rankLeft, rankRight;
    
    //  Let each process know how many ids to expect to receive.
    numToRecv   = new uint[size];
    MPI_Allgather(&interNodeChangesBufferSize, 1, MPI_UNSIGNED,
                  numToRecv,                   1, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    //  Determine neighboring processors.
    if (rank == size - 1)   rankRight = MPI_PROC_NULL;
    else                    rankRight = rank + 1;
    if (rank == 0)          rankLeft  = MPI_PROC_NULL;
    else                    rankLeft  = rank - 1;
    
    //  Send and receive particle ids.
    if (rank % 2)
    {   recv(rankLeft);     send(rankRight);
        recv(rankRight);    send(rankLeft); }
    else
    {   send(rankRight);    recv(rankLeft);
        send(rankLeft);     recv(rankRight); } }

/*  resetVariables()
 *
 *  Reset the values of variables unique to each iteration.
 */
inline void resetVariables()
{   interNodeChangesBufferSize  = 0;
    interNodeStatesBufferSize   = 0;
    receiveCount = 0;
    internalChangesBuffer.clear(); }

/*  hasDissolved(ParticlePtr ptr)
 *  
 *  Determine if the particle indicated by ptr has dissolved.
 */
inline bool hasDissolved(ParticlePtr ptr)
{   return (ptr->state >= dissolnStates); }

/*  placeOnSurface(Particle& p)
 *  
 *  Add a particle to the surface and affect the number of nearest neighbors
 *  accordingly. *** shouldn't do the latter here w/o more info (who dissolved?)
 */
inline void placeOnSurface(Particle& p)
{   if (!hasDissolved(&p) && !p.onSurface)
    {   surface.push_back(&p);
        p.onSurface = true; } }

/*  updateLocal()
 *  
 *  Update the state of the neighbors of each local dissolved particle.
 */
inline void updateLocal()
{   uint count = 0, initialN;
    
    for (uint i = 0; i < internalChangesBuffer.size(); i++)
    {   initialN = accumulate(internalChangesBuffer[i]->initialN.begin(),
                              internalChangesBuffer[i]->initialN.begin() + dissolnStates, 0);
        for (uint n = 0; n < initialN; n++)
        {   pmapiter iter = pmap.find(internalChangesBuffer[i]->getNList()[n]);
            if (iter != pmap.end())
            {   placeOnSurface(iter->second);
                iter->second.countN[internalStatesBuffer[i].oldState]--;
                iter->second.countN[internalStatesBuffer[i].newState]++;
                count++; } } } }

/*  updateExternal()
 *
 *  Update the state of the neighbors of each remote dissolved particle.
 */
inline void updateExternal()
{   uint count = 0; //***
    if (verbose && count > 0) cout << rank << ": received " << receiveCount
                                   << " external particles.\n";
    
    for (uint i = 0; i < receiveCount; i++)
    {   pmapiter iter = pmap.find(receiveBuffer[i]);
        if (iter != pmap.end())
        {   placeOnSurface(iter->second);
            iter->second.countN[interNodeStatesBuffer[i].oldState]--;
            iter->second.countN[interNodeStatesBuffer[i].newState]++;
            count++; } } }

#endif  /* _PROCESS_HELPERS_H */

