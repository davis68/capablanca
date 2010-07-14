/** process.cpp
 *  29 Oct 2009--14 Jun 2010
 *  Minas Charalambides and Neal Davis
 *  
 */

#include <iostream>
#include <list>
#include <map>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "definitions.h"
#include "error.h"
#include "particle.h"
#include "process.h"
#include "statistics.h" //  FIXME debug

using namespace std;

//  External global variables from global.cpp
extern bool verbose;

extern uint numStates,
            dissolnStates,
            maxNN;

extern double       T,
                    phi;
extern Reaction    *rxn;

extern uint outputInterval,
            tmax,
            ranseed;
extern bool deposition;

extern int  rank,
            size;

extern coord_t myMaxZ;

extern ParticleMap              pmap;
extern vector<vector<double> >  probA;
extern vector<vector<double> >  probB;

extern list<Particle*> surfaceA;
extern list<Particle*> surfaceB;

extern uint    *oldParticleCount;
extern uint    *newParticleCount;
extern uint    *myParticleCount;
extern Vector   rates;

//  Global variables
static id_t     sendChangesBuffer[MAX_ARRAY_SIZE];
static state_t  sendStatesBuffer[2 * MAX_ARRAY_SIZE];
static uint     sendBufferSize;

static vector<Particle*>    internalChangesBuffer;
static vector<state_t>      internalStatesBuffer;
static vector<Particle*>    boundaryChangesBuffer;
static vector<state_t>      boundaryStatesBuffer;

static id_t     recvChangesBuffer[MAX_ARRAY_SIZE];
static state_t  recvStatesBuffer[2 * MAX_ARRAY_SIZE];
static uint     recvChangesCount;

/** process()
 *  
 *  Cycle through the iterations of the system, updating each particle both
 *  locally and remotely.
 */
void process()
{   calcProbs();
    findSurface();
    
    ranseed = 43; //FIXME:(unsigned int) time(0);
    oldParticleCount= new uint[numStates];
    newParticleCount= new uint[numStates];
    myParticleCount = new uint[numStates];
    for (uint i = 0; i < numStates; i++)
    {   oldParticleCount[i] = 0;
        newParticleCount[i] = 0;
        myParticleCount[i]  = 0; }
    rates.resize(numStates, 0);
    resetVariables();
    
    for (uint t = 0; t <= tmax; t++)
    {   //  Periodically output calculation status and system statistics.
        if (!(t % outputInterval)) collateStatistics(t);
        if (!(t % outputInterval)) outputSurface(surfaceA, t); //FIXME
        
        //  Apply rules to system particle-by-particle.
        for (list<Particle*>::iterator iter = surfaceA.begin(); iter != surfaceA.end(); iter++)
        {   transitionParticle(*iter); }
        if (deposition)
        {   for (list<Particle*>::iterator iter = surfaceB.begin(); iter != surfaceB.end(); iter++)
            {   transitionParticle(*iter); } }
        
        //  Pack the surface down and indicate changes to all processes.
        setParticleStates();
        updateLocalNeighbors();
        exchangeInterNodeChanges();
        updateExternalNeighbors();
        packSurface();
        resetVariables();
        
        if (verbose && !(t % outputInterval) && !rank)
        {   cout << "  " << 100 * (double) t / (double) tmax << "% complete.\n";
            cout.flush(); } }
    
    //  Output the data from the last time step.
    collateStatistics(tmax);
    
    delete[] oldParticleCount;
    delete[] newParticleCount;
    delete[] myParticleCount; }

/** calcProbs()
 *  
 *  As there are only maxNN discrete states permissible for any given
 *  temperature, pre-calculate the probabilities for each arrangement.
 *  See for details:  Alekseev Yu V, Alekseev G Yu, Bityurin V A. (2002).
 *  Atomic--topological and statistical background for a correct dissolution
 *  theory of crystal substances. Protection of Metals.  38(6):517-529.
 *  
 *  The probabilities of transition are stored in the globally-linked variables
 *  probA and probB, with the first index iterating over the number of states
 *  and the second index iterating over the number of neighbors (up to maxNN).
 */
inline void calcProbs()
{   /// Dissolution:
    probA.resize(numStates, vector<double>(maxNN, 0.0));
    //  Calculate values for each state and nn condition.
    for (uint i = 0; i < numStates; i++)
    {   for (uint j = 0; j < maxNN; j++)
        {   probA[i][j] = rxn[i].prefactor
                        * exp(-rxn[i].alpha * j * rxn[i].E_s      / (k_B * T))
                        * exp(-rxn[i].alpha     * rxn[i].E_r      / (k_B * T))
                        * exp(-rxn[i].alpha * rxn[i].z * ec * phi / (k_B * T)); } }
    
    /// Deposition:
    probB.resize(numStates, vector<double>(maxNN, 0.0));
    //  Calculate values for each state and nn condition.
    for (uint i = 0; i < numStates; i++)
    {   for (uint j = 0; j < maxNN; j++)
        {   probB[i][j] = rxn[i].prefactor
                        * exp(rxn[i].beta * j * rxn[i].E_s      / (k_B * T))
                        * exp(rxn[i].beta     * rxn[i].E_r      / (k_B * T))
                        * exp(rxn[i].beta * rxn[i].z * ec * phi / (k_B * T)); } } }

/** findSurface()
 *  
 *  Determine the particles contiguous with those whose z-coordinate is less
 *  than SURFACE_CUTOFF from the global maximum z-coordinate and add them to the
 *  appropriate surface.
 *  
 *  Vacancy defects within SURFACE_CUTOFF of the top will be incorrectly added
 *  as part of the surface.
 */
void findSurface()
{   //  Determine global maximum z-coordinate.
    coord_t  maxZ;
    MPI_Allreduce(&myMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Add particles within SURFACE_CUTOFF of the global maximum z-coordinate
    //  to the appropriate surface.  (At this point, the surfaces are empty, so
    //  we don't need to worry about removing them from any surfaces first.)
    /*FIXME:static uint numNN;
    for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if (maxZ - iter->second.z > SURFACE_CUTOFF) continue;
        
        if (!hasDissolved(&(iter->second))) placeOnSurfaceA(iter->second);
        else                                placeOnSurfaceB(iter->second); }*/
    ParticleMap::iterator mapIter;
    //ParticleMap::iterator iter1 = pmap.begin();
    mapIter = pmap.begin();
    /*for (iter1 = pmap.begin(); iter1 != pmap.end(); iter1++)
    {   if (iter1->second.) break; }*/ //TODO
    placeOnSurfaceA(mapIter->second);
    
    //  Now expand the surface from these starting points.
    for (list<Particle*>::iterator iter = surfaceA.begin(); iter != surfaceA.end(); iter++)
    {   //  Attempt to add the neighboring particles to the surface.
        for (int i = 0; i < (*iter)->neighbors.size(); i++)
        {   mapIter = pmap.find((*iter)->neighbors[i]);
            if (mapIter == pmap.end()) continue;
            
            expandSurfaces(mapIter->second); } }
    for (list<Particle*>::iterator iter = surfaceB.begin(); iter != surfaceB.end(); iter++)
    {   //  Attempt to add the neighboring particles to the surface.
        for (int i = 0; i < (*iter)->neighbors.size(); i++)
        {   mapIter = pmap.find((*iter)->neighbors[i]);
            if (mapIter == pmap.end()) continue;
            
            expandSurfaces(mapIter->second); } } }

/** resetVariables()
 *
 *  Reset the values of iteration-dependent variables.
 */
void resetVariables()
{   sendBufferSize = 0;
    internalChangesBuffer.clear();
    internalStatesBuffer.clear();
    boundaryChangesBuffer.clear();
    boundaryStatesBuffer.clear(); }

/** transitionParticle(Particle* ptr)
 *  
 *  Calculate the probability that the particle indicated by ptr has its state
 *  altered by the currently applicable rule and update the state of that
 *  particle.  Indicate that the data of its neighbors also need to be updated.
 */
void transitionParticle(Particle* ptr)
{   state_t  oldState,
             newState;
    uint     numNN = accumulate(ptr->countN.begin(), ptr->countN.begin() + dissolnStates, 0);
    
    if (!hasDissolved(ptr))
    {   /// Dissolution:
        if (uniformRand() < probA[ptr->state][numNN])
        {   //  Set the particle states.
            oldState    = ptr->state;
            newState    = rxn[ptr->state].newState;
            ptr->state  = newState;
            
            //  Queue the particle so that the neighbors can be updated.
            if (ptr->onBoundary) updateBoundaryParticle(ptr, oldState, newState);
            else                 updateInternalParticle(ptr, oldState, newState); } }
    else
    {   /// Deposition:
        if (uniformRand() < probB[ptr->state][numNN])
        {   //  Set the particle states.
            oldState    = ptr->state;
            newState    = rxn[ptr->state].newState;
            ptr->state  = newState;
            
            //  Queue the particle so that the neighbors can be updated.
            if (ptr->onBoundary) updateBoundaryParticle(ptr, oldState, newState);
            else                 updateInternalParticle(ptr, oldState, newState); } } }

/** updateInternalParticle()
 *  
 *  Add data for the particle indicated by ptr to the internal buffers.
 */
inline void updateInternalParticle(Particle* ptr, state_t oldState, state_t newState)
{   //  Queue the particle for the state transition in setParticleStates().
    internalChangesBuffer.push_back(ptr);
    internalStatesBuffer.push_back(oldState);
    internalStatesBuffer.push_back(newState); }

/** updateBoundaryParticle()
 *  
 *  Add data for the particle indicated by ptr to the external buffers.  A
 *  particle which is on the boundary will have its neighbors' data updated
 *  when all external particles are updated, in updateExternal().
 */
void updateBoundaryParticle(Particle* ptr, state_t oldState, state_t newState)
{   //  Queue the particle for the state transition in setParticleStates().
    boundaryChangesBuffer.push_back(ptr);
    boundaryStatesBuffer.push_back(oldState);
    boundaryStatesBuffer.push_back(newState);
    
    //  Write the particle's neighbors, which must be updated, to sendChangesBuffer.
    for (uint i = 0; i < ptr->neighbors.size(); i++)
    {   sendChangesBuffer[sendBufferSize + i] = ptr->neighbors[i];
        
        //  For each neighbor of the particle indicated by ptr, we must indicate the
        //  previous and current state of the particle.
        sendStatesBuffer[2 * sendBufferSize + 2 * i    ] = oldState;
        sendStatesBuffer[2 * sendBufferSize + 2 * i + 1] = newState; }
    
    sendBufferSize += ptr->neighbors.size(); }

/** setParticleStates()
 *  
 *  Update old particle states to new ones.
 */
void setParticleStates()
{   //  Update the state of local particles.
    ParticleMap::iterator mapIter;
    for (int i = 0; i < internalChangesBuffer.size(); i++)
    {   mapIter = pmap.find(internalChangesBuffer[i]->id);
        mapIter->second.state = internalStatesBuffer[2 * i + 1]; }
    
    //  Update the state of boundary particles.
    for (int i = 0; i < boundaryChangesBuffer.size(); i++)
    {   mapIter = pmap.find(boundaryChangesBuffer[i]->id);
        mapIter->second.state = boundaryStatesBuffer[2 * i + 1]; } }

/** exchangeInterNodeChanges()
 *  
 *  Exchange data in order to properly update boundary particles.
 */
void exchangeInterNodeChanges()
{   //  Let each process know how many particles to expect to receive data for.
    uint *numToRecv;
    numToRecv   = new uint[size];
    MPI_Allgather(&sendBufferSize, 1, MPI_UNSIGNED,
                  numToRecv,       1, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    //  Create sendCount and recvCount, sendDispls and recvDispls buffers.
    //  These are for the MPI_Alltoallv() calls below, in which the count and
    //  displacements of the data are required.  We reallocate each time instead
    //  of once statically to prevent memory leaks.
    int *sendCount,  *recvCount,
        *sendDispls, *recvDispls,
         sendCt = 0,  recvCt = 0;
    sendCount  = new int[size]; recvCount  = new int[size];
    sendDispls = new int[size]; recvDispls = new int[size];
    
    //  Calculate the count and displacements for the buffers based on the
    //  number of expected particle data.  Thus xxxxCt acts as a running total
    //  to correctly indicate the number already expected in xxxxCount/xxxxDispls.
    for (int i = 0; i < size; i++)
    {   sendCount[i]  = 0;
        sendDispls[i] = 0;
        //  Only processes immediately adjacent to this one have particles bordering it.
        if (i >= rank - 1 && i <= rank + 1)
        {   sendCount[i] = sendBufferSize;
            sendCt += sendBufferSize; } }
    for (int i = 0; i < size; i++)
    {   recvCount[i]  = 0;
        recvDispls[i] = recvCt;
        //  Only processes immediately adjacent to this one have particles bordering it.
        if (i >= rank - 1 && i <= rank + 1)
        {   recvCount[i] = numToRecv[i];
            recvCt += numToRecv[i]; } }
    
    //  Update the global count of changes received (the proper buffer size).
    recvChangesCount = 0;
    for (int i = 0; i < size; i++)
    {   recvChangesCount += recvCount[i]; }
    
    //  Send and receive lists of particle neighbor ids.
    MPI_Alltoallv(sendChangesBuffer, sendCount, sendDispls, MPI_UNSIGNED,
                  recvChangesBuffer, recvCount, recvDispls, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
    
    //  Send and receive the corresponding particle state data (of which there
    //  are twice as many entries).
    for (int i = 0; i < size; i++)
    {   sendCount[i]  *= 2; recvCount[i]  *= 2;
        sendDispls[i] *= 2; recvDispls[i] *= 2; }
    MPI_Alltoallv(sendStatesBuffer, sendCount, sendDispls, MPI_UNSIGNED,
                  recvStatesBuffer, recvCount, recvDispls, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
    
    delete[] numToRecv;
    delete[] sendCount;  delete[] recvCount;
    delete[] sendDispls; delete[] recvDispls; }

/** updateLocalNeighbors()
 *  
 *  Update the state of the neighbors of each local dissolved particle by
 *  querying the list of neighbors and altering their nearest neighbor counts.
 *  
 *  The particle's state is changed in setParticleState(); ramifications of
 *  that transition are considered here.  Neighbors of boundary particles are
 *  updated in updateExternalNeighbors().
 */
void updateLocalNeighbors()
{   ParticleMap::iterator mapIter;
    
    //  Update the neighbors' neighbor counts for non-boundary particles.
    for (int i = 0; i < internalChangesBuffer.size(); i++)
    {   for (vector<uint>::iterator iter1 = internalChangesBuffer[i]->neighbors.begin(); iter1 != internalChangesBuffer[i]->neighbors.end(); iter1++)
        {   //  Find the current neighbor of interest.  (If this returns
            //  pmap.end() we will throw an exception, but this is a local update,
            //  so if it returns pmap.end() there is a problem with the buffer.)
            mapIter = pmap.find(*iter1);
            
            //  Update the neighbor count from the state buffer.  (We increment
            //  first to prevent underflow if the states are identical.)
            (mapIter->second.countN[internalStatesBuffer[2 * i + 1]])++;
            (mapIter->second.countN[internalStatesBuffer[2 * i    ]])--; } } }

/** updateExternalNeighbors()
 *  
 *  Update the neighbor count of the neighbors of each received dissolved
 *  particle by directly altering their nearest neighbor counts.  (This includes
 *  both remote and local boundary particles.)
 */
void updateExternalNeighbors()
{   ParticleMap::iterator mapIter;
    
    for (uint i = 0; i < recvChangesCount; i++)
    {   //  Find the current particle of interest.
        mapIter = pmap.find(recvChangesBuffer[i]);
        if (mapIter == pmap.end()) continue;    //  The remote particle may not be found on this process.
        
        //  Update the neighbor count from the state buffer.  (We increment
        //  first to prevent underflow if the states are identical.)
        (mapIter->second.countN[recvStatesBuffer[2 * i + 1]])++;
        (mapIter->second.countN[recvStatesBuffer[2 * i    ]])--; } }

/** packSurface()
 *  
 *  Update the surface, removing detached or covered particles and expanding
 *  into new contiguous areas.  Definitionally, if a particle's state has
 *  changed in this time step, it is part of the surface.  An STL list of
 *  particles is used as the surfaces because of this function, which would
 *  otherwise invalidate the pointers and iterators as the alterations are made.
 *  
 *  Throughout, I use detached to mean a particle with zero dissolved
 *  neighbors and covered to mean a particle with maxNN non-dissolved neighbors.
 */
void packSurface()
{   //  Local particles which could react are by definition part of the surface
    //  already, but we need to check the neighbors if the particle dissolved
    //  or deposited.
    ParticleMap::iterator mapIter;
    for (vector<Particle*>::iterator iter = internalChangesBuffer.begin(); iter != internalChangesBuffer.end(); iter++)
    {   for (vector<uint>::iterator iter1 = (*iter)->neighbors.begin(); iter1 != (*iter)->neighbors.end(); iter1++)
        {   mapIter = pmap.find(*iter1);
            if (mapIter == pmap.end()) continue;
            
            expandSurfaces(mapIter->second); } }
    
    //  Add the neighbors of each boundary particle to the appropriate surface.
    for (uint i = 0; i < recvChangesCount; i++)
    {   //  Find the current particle of interest.
        mapIter = pmap.find(recvChangesBuffer[i]);
        if (mapIter == pmap.end()) continue;
        
        expandSurfaces(mapIter->second); }
    
    //  Remove dissolved or covered particles from surfaceA.
    uint numNN;
    for (list<Particle*>::iterator iter = surfaceA.begin(); iter != surfaceA.end(); iter++)
    {   numNN = accumulate((*iter)->countN.begin(), (*iter)->countN.begin() + dissolnStates, 0);
        if ( hasDissolved(*iter) || numNN == maxNN)
        {   (*iter)->onSurface = notOnSurface;
            iter = surfaceA.erase(iter);
            iter--; } } //  I've tested this separately, and it doesn't cause trouble at list::begin().
    
    //  Remove deposited or detached particles from surfaceB.
    for (list<Particle*>::iterator iter = surfaceB.begin(); iter != surfaceB.end(); iter++)
    {   numNN = accumulate((*iter)->countN.begin(), (*iter)->countN.begin() + dissolnStates, 0);
        if (!hasDissolved(*iter) || numNN == 0)
        {   (*iter)->onSurface = notOnSurface;
            iter = surfaceB.erase(iter);
            iter--; } } }

/** expandSurfaces()
 *  
 *  Attempt to add the particle to each surface.
 */
inline void expandSurfaces(Particle &p)
{   expandSurfaceA(p);
    expandSurfaceB(p); }

/** expandSurfaceA()
 *  
 *  Add a particle to surfaceA and recursively affect the nearest neighbors
 *  accordingly.  The particle should have previously been removed from any
 *  other surface.
 *  
 *  surfaceA consists of particles which fulfill three criteria:
 *      1.  they are not in a dissolved state;
 *      2.  they have an incomplete complement of nearest neighbors; and
 *      3.  they are contiguous with the original surfaceA or a piece of it.
 */
void expandSurfaceA(Particle& p)
{   uint numNN;
    
    //  If a particle is ineligible for addition to surfaceA, whether by already
    //  belonging to surfaceA or not fulfilling the three criteria, then return
    //  without adding it to surfaceA.
    if (p.onSurface == onSurfaceA)  return;
    if ( hasDissolved(&p))          return;
    numNN = accumulate(p.countN.begin(), p.countN.begin() + dissolnStates, 0);
    if (numNN == maxNN)             return;
    
    placeOnSurfaceA(p);
    
    //  Attempt to add the neighboring particles to the surfaces.
    ParticleMap::iterator mapIter;
    for (int i = 0; i < p.neighbors.size(); i++)
    {   mapIter = pmap.find(p.neighbors[i]);
        if (mapIter == pmap.end()) continue;
        
        expandSurfaces(mapIter->second); } }

/** expandSurfaceB()
 *  
 *  Add a particle to surfaceB and recursively affect the nearest neighbors
 *  accordingly.  The particle should have previously been removed from any
 *  other surface.
 *  
 *  surfaceB consists of particles which fulfill three criteria:
 *      1.  they are in a dissolved state;
 *      2.  they have a nonzero complement of nearest neighbors; and
 *      3.  they have dissolved from the original surfaceA or a piece of it.
 */
void expandSurfaceB(Particle& p)
{   uint numNN;
    
    //  If a particle is ineligible for addition to surfaceA, whether by already
    //  belonging to surfaceB or not fulfilling the three criteria, then return
    //  without adding it to surfaceB.
    if (p.onSurface == onSurfaceB)  return;
    if (!hasDissolved(&p))          return;
    numNN = accumulate(p.countN.begin(), p.countN.begin() + dissolnStates, 0);
    if (numNN == 0)                 return;
    
    placeOnSurfaceB(p);
    
    //  Attempt to add the neighboring particles to the surfaces.
    ParticleMap::iterator mapIter;
    for (int i = 0; i < p.neighbors.size(); i++)
    {   mapIter = pmap.find(p.neighbors[i]);
        if (mapIter == pmap.end()) continue;
        
        expandSurfaces(mapIter->second); } }

/** placeOnSurfaceA()
 *  
 *  Add a particle to the surface---it should already have been removed from any
 *  other surface.
 */
inline void placeOnSurfaceA(Particle& p)
{   p.onSurface = onSurfaceA;
    surfaceA.push_back(&p); }

/** placeOnSurfaceB()
 *  
 *  Add a particle to the surface---it should already have been removed from any
 *  other surface.
 */
inline void placeOnSurfaceB(Particle& p)
{   p.onSurface = onSurfaceB;
    surfaceB.push_back(&p); }

/** hasDissolved(Particle* ptr)
 *  
 *  Indicate if the particle indicated by ptr has dissolved.
 */
inline bool hasDissolved(Particle* ptr)
{   return (ptr->state >= dissolnStates); }

