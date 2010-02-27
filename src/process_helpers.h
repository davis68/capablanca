/*
The chief problem right now is that only when a particle dissolves is its state
abroad updated.  We need a plan to update the information every time a particle
changes at all instead.
*/
#ifndef _PROCESS_HELPERS_H
#define _PROCESS_HELPERS_H

#include <omp.h>
#include <mpi.h>
#include <vector>
#include <map>
#include "particle.h"
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

struct States
{   state_t oldState, newState; };

typedef ParticleMap::iterator pmapiter;

static vector<ParticlePtr>  surface;

static id_t     interNodeChangesBuffer[MAX_ARRAY_SIZE]; //*** best approach?
static uint     interNodeChangesBufferSize;
static States   interNodeStatesBuffer[MAX_ARRAY_SIZE];
static uint     interNodeStatesBufferSize;
static id_t     interNodeSendBuffer[MAX_ARRAY_SIZE];

static vector<ParticlePtr>  internalChangesBuffer;
static vector<States>       internalStatesBuffer;

static id_t receiveBuffer[MAX_ARRAY_SIZE];
static uint receiveCount;

static coord_t  maxCoord;
static coord_t  globalMaxCoord;
static uint     NUM_OF_NEIGHBORS;

static uint    *recvCountList,
               *sendCountList;

/*  hasDissolved(ParticlePtr ptr)
 *  
 *  Determine if the particle indicated by ptr has dissolved.
 */
inline bool hasDissolved(ParticlePtr ptr)
{   return (ptr->state >= dissolnStates); }

/*  findSurface()
 *  
 *  Determine the particles whose z-coordinate is less than SURFACE_CUTOFF from
 *  the global maximum z-coordinate.
 */
inline void findSurface()
{   //  Determine local (and thence global) maximum z-coordinates.
    maxCoord = -1e8;
    for(pmapiter iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if(iter->second.z > maxCoord)   maxCoord = iter->second.z; }
    MPI_Allreduce(&maxCoord, &globalMaxCoord, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Add particles within SURFACE_CUTOFF of the top to the surface.
    for(pmapiter iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if(globalMaxCoord - iter->second.z <= SURFACE_CUTOFF) surface.push_back(&(iter->second)); }
}

/*  updateOnBoundary(ParticlePtr ptr)
 *  
 *  Add the particle indicated by ptr to the interNodeChangesBuffer.
 */
inline void updateOnBoundary(ParticlePtr ptr)
{   uint initialN = accumulate(ptr->initialN.begin(), ptr->initialN.begin() + dissolnStates, 0);
    memcpy(&(interNodeChangesBuffer[interNodeChangesBufferSize]),
           ptr->getNList(), initialN * sizeof(id_t));
    interNodeChangesBufferSize += initialN;
    memcpy(&(interNodeStatesBuffer[interNodeStatesBufferSize]),
           ptr->getNList(), initialN * sizeof(id_t));
    interNodeStatesBufferSize += initialN; }

/*  updateOffBoundary(ParticlePtr ptr)
 *
 *  Add the particle indicated by ptr to the internalChangesBuffer.
 */
inline void updateOffBoundary(ParticlePtr ptr)
{   internalChangesBuffer.push_back(ptr);
    States st;
    st.oldState = ptr->state;   st.newState = ptr->state; //*** based on rule transition?
    internalStatesBuffer.push_back(st); }

/*  recv()
 *
 *  Receive all changes from other processors.
 */
inline void recv()
{   static MPI_Status status;
    MPI_Recv(&(receiveBuffer[receiveCount]), MAX_ARRAY_SIZE, MPI_UNSIGNED,
             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    receiveCount += status.MPI_TAG; //  The tag is the number of particles sent.
}

/*  send(const int dest)
 *  
 *  Send all changes for dest to dest.
 */
inline void send(const int dest)
{   static id_t sendBuffer[MAX_ARRAY_SIZE];
    uint sendBufferSize = 0;
    memcpy(&(sendBuffer[sendBufferSize]), &(interNodeChangesBuffer[0]),
           interNodeChangesBufferSize * sizeof(id_t));
    sendBufferSize += interNodeChangesBufferSize;
    MPI_Send(sendBuffer, sendBufferSize, MPI_UNSIGNED,
             dest, sendBufferSize, MPI_COMM_WORLD);
}

/*  exchangeInterNodeChanges()
 *  
 *  Initiate the exchange of inter-node data by even processors first, then odd
 *  ones. *** make non-blocking (issend?)
 */
inline void exchangeInterNodeChanges()
{   
    /*MPI_Alltoall(sendCountList, size, MPI_INT, recvCountList, size, MPI_INT,
                 MPI_COMM_WORLD);*/
    
    if(rank % 2)
    {   recv();
        send(rank - 1);
        if(rank < size - 1)
        {   recv();
            send(rank + 1); } }
    else
    {   if(rank < size - 1)
        {   send(rank + 1);
            recv(); }
        if(rank > 0)
        {   send(rank - 1);
            recv(); } } }

/*  placeOnSurface(Particle& p)
 *  
 *  Add a particle to the surface and affect the number of nearest neighbors
 *  accordingly. *** shouldn't do the latter here w/o more info (who dissolved?)
 */
/***inline void placeOnSurface(Particle& p)
{   if (!hasDissolved(&p))
    {   p.countN--;
        if (!p.onSurface) { surface.push_back(&p); p.onSurface = true; } } }*/
inline void placeOnSurface(Particle& p)
{   if (!hasDissolved(&p) && !p.onSurface)
    {   surface.push_back(&p);
        p.onSurface = true; } }

/*  updateLocal()
 *  
 *  Update the state of the neighbors of each local dissolved particle.
 */
inline void updateLocal()
{   uint    count = 0,
            initialN;
    
    for (uint i = 0; i < internalChangesBuffer.size(); i++)
    {   initialN = accumulate(internalChangesBuffer[i]->initialN.begin(),
                              internalChangesBuffer[i]->initialN.begin() + dissolnStates, 0);
        for (uint n = 0; n < initialN; n++)
        {   pmapiter iter = pmap.find(internalChangesBuffer[i]->getNList()[n]);
            if (iter != pmap.end())
            {   placeOnSurface(iter->second);
                iter->second.countN[internalStatesBuffer[i].oldState]--;
                iter->second.countN[internalStatesBuffer[i].newState]++;
                count++; } } }  //***
    #ifdef PROC_DEBUG
    cerr<<"rank: "<<rank<<" updated "<<count<<" thread internal particles"<<endl;
    #endif
}

/*  updateExternal()
 *
 *  Update the state of the neighbors of each remote dissolved particle.
 */
inline void updateExternal()
{   uint count=0; //***
    if (verbose && count > 0) cout << rank << ": received " << receiveCount
                                   << " external particles." << endl;
    
    for (uint i = 0; i < receiveCount; i++)
    {   pmapiter iter = pmap.find(receiveBuffer[i]);
        if (iter != pmap.end())
        {   placeOnSurface(iter->second);
            iter->second.countN[interNodeStatesBuffer[i].oldState]--;
            iter->second.countN[interNodeStatesBuffer[i].newState]++;
            count++; } }  //***
    #ifdef PROC_DEBUG
    cerr<<"rank: "<<rank<<" updated "<<count<<" external particles"<<endl;
    #endif
}

#endif  /* _PROCESS_HELPERS_H */

