/** process.cpp
 *  29 Oct 2009--07 Jun 2010
 *  Minas Charalambides and Neal Davis
 *  
 */

#include <iostream> //  Used for debugging.
#include <list>     //  
#include <map>      //  
#include <cmath>    //  
#include <mpi.h>    //  
#include <string.h> //  
#include <cstdlib>  //  
#include <ctime>    //  Used to initialize random seed ranseed.
#include <vector>   //  

#include <boost/random.hpp>             //  These may require explicit inclusion of the Boost C++
#include <boost/random/uniform_01.hpp>  //  library locations when compiling at the command line.

#include "definitions.h"
#include "error.h"
#include "particle.h"
#include "process.h"
#include "statistics.h" // FIXME debug

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
static uint     sendSurfaceBuffer[MAX_ARRAY_SIZE];
static uint     sendChangesBufferSize;

static vector<Particle*>    internalChangesBuffer;
static vector<state_t>      internalStatesBuffer;

static id_t     recvChangesBuffer[MAX_ARRAY_SIZE];
static state_t  recvStatesBuffer[2 * MAX_ARRAY_SIZE];
static uint     recvSurfaceBuffer[MAX_ARRAY_SIZE];
static uint     recvChangesCount;

/*  process()
 *  
 *  Cycle through the iterations of the system, updating each particle both
 *  locally and remotely.
 */
void process()
{   calcProbs();
    findSurface();
    
    ranseed = (unsigned int) time(0);
    oldParticleCount= new uint[numStates];
    newParticleCount= new uint[numStates];
    myParticleCount = new uint[numStates];
    for (uint i = 0; i < numStates; i++)
    {   oldParticleCount[i] = 0;
        newParticleCount[i] = 0;
        myParticleCount[i]  = 0; }
    rates.resize(numStates, 0);
    
    for (uint t = 0; t <= tmax; t++)
    {   resetVariables();
        
        //  Apply rules to system particle-by-particle.
        for (list<Particle*>::iterator iter = surfaceA.begin(); iter != surfaceA.end(); iter++)
        {   updateState(*iter); }
        if (deposition)
        {   for (list<Particle*>::iterator iter = surfaceB.begin(); iter != surfaceB.end(); iter++)
            {   updateState(*iter); } }
        
        //  Update particle maps and surfaces locally and across all processors.
        exchangeInterNodeChanges();
        updateLocal();
        updateExternal();
        packSurface();
        
        if (!(t % outputInterval)) { collateStatistics(t); outputSurface(surfaceA, t); } //FIXME
        MPI_Barrier(MPI_COMM_WORLD);    //  Remain synchronous.
        
        if (verbose && !(t % outputInterval) && !rank)
        {   cout << "  " << 100 * (double) t / (double) tmax << "% complete.\n";
            cout.flush(); } }
    
    delete[] oldParticleCount;
    delete[] newParticleCount;
    delete[] myParticleCount; }

/*  hasDissolved(Particle* ptr)
 *  
 *  Determine if the particle indicated by ptr has dissolved.
 */
inline bool hasDissolved(Particle* ptr)
{   return (ptr->state >= dissolnStates); }

/*  findSurface()
 *  
 *  Determine the particles whose z-coordinate is less than SURFACE_CUTOFF from
 *  the global maximum z-coordinate and add them to the appropriate surface.
 */
void findSurface()
{   //  Determine local (and thence global) maximum z-coordinates.
    static coord_t  maxCoord,
                    globalMaxCoord;
    maxCoord = -1e8;
    for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
    {   if (iter->second.z > maxCoord)   maxCoord = iter->second.z; }
    MPI_Allreduce(&maxCoord, &globalMaxCoord, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Add particles within SURFACE_CUTOFF of the top to the appropriate surface.
    static uint numNN;
    for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
    {   numNN = accumulate(iter->second.countN.begin(),
                           iter->second.countN.begin() + dissolnStates, 0);
        
        if (globalMaxCoord - iter->second.z > SURFACE_CUTOFF || 
            iter->second.onSurface != notOnSurface || numNN >= maxNN) continue;
        
        if (!hasDissolved(&(iter->second)))
        //  XXX  Although it would be preferable to use something like
        //      placeOnSurfaceA(iter->second);,
        //  it overflows memory quickly and causes a core dump.
        {   surfaceA.push_back(&(iter->second));
            iter->second.onSurface = onSurfaceA; }
        else
        {   surfaceB.push_back(&(iter->second));
            iter->second.onSurface = onSurfaceB; } } }

/*  updateOnBoundary()
 *  
 *  Add the particle indicated by ptr to the sendChangesBuffer.
 */
void updateOnBoundary(Particle* ptr, state_t oldState, state_t newState)
{   //  Write the particle indicated by ptr for local updating as well.
    updateOffBoundary(ptr, oldState);
    
    static uint initialNCount;
    initialNCount = accumulate(ptr->initialN.begin(),
                               ptr->initialN.begin() + dissolnStates, 0);
    
    //  Write the particle's neighbors, which must be updated, to a buffer.
    memcpy(&(sendChangesBuffer[sendChangesBufferSize]),
           ptr->getNList(), initialNCount * sizeof(id_t));
    
    //  Write the particle's state changes, which must be updated for each
    //  neighbor, to a buffer.
    for (uint i = 0; i < initialNCount; i++)
    {   sendStatesBuffer[2 * sendChangesBufferSize + 2 * i    ] = oldState;
        sendStatesBuffer[2 * sendChangesBufferSize + 2 * i + 1] = newState; }
    
    //  Write the particle's surface to a buffer.
    for (uint i = 0; i < initialNCount; i++)
    {   sendSurfaceBuffer[sendChangesBufferSize + i] = ptr->onSurface; }
    
    sendChangesBufferSize += initialNCount; }

/*  updateOffBoundary()
 *
 *  Add the particle indicated by ptr to the internalChangesBuffer.
 */
inline void updateOffBoundary(Particle* ptr, state_t oldState)
{   //  Write the pointer to the particle to a buffer so the neighbors can be updated.
    internalChangesBuffer.push_back(ptr);
    
    //  Write the particle's state changes, which must be updated for each
    //  neighbor, to a buffer.
    internalStatesBuffer.push_back(oldState); }

/*  exchangeInterNodeChanges()
 *  
 *  Initiate the exchange of inter-node data by even processors first, then odd
 *  ones.
 */
void exchangeInterNodeChanges()
{   int rankLeft, rankRight;
    
    //  Let each process know how many ids to expect to receive.
    static uint *numToRecv;
    numToRecv   = new uint[size];
    MPI_Allgather(&sendChangesBufferSize, 1, MPI_UNSIGNED,
                  numToRecv,              1, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    //  Determine neighboring processors.
    rankRight = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
    rankLeft  = (rank == 0)        ? MPI_PROC_NULL : rank - 1;
    
    //  Create sendCount and recvCount, sendDispls and recvDispls buffers.
    static int *sendCount,
               *recvCount,
               *sendDispls,
               *recvDispls;
    sendCount = new int[size];
    recvCount = new int[size];
    sendDispls= new int[size];
    recvDispls= new int[size];
    int sendCt = 0,
        recvCt = 0;
    for (int i = 0; i < size; i++)
    {   sendCount[i]  = 0;
        sendDispls[i] = 0;
        if (i == rankRight || i == rankLeft)
        {   sendCount[i] = sendChangesBufferSize;
            sendDispls[i] = sendCt;
            sendCt += sendChangesBufferSize; } }
    for (int i = 0; i < size; i++)
    {   recvCount[i] = 0;
        recvDispls[i] = 0;
        if (i == rankRight || i == rankLeft)
        {   recvCount[i] = numToRecv[i];
            recvDispls[i] = recvCt;
            recvCt += numToRecv[i]; } }
    
    //  Send and receive particle neighbor ids.
    MPI_Alltoallv(sendChangesBuffer, sendCount, sendDispls, MPI_UNSIGNED,
                  recvChangesBuffer, recvCount, recvDispls, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
    //  Send and receive particle surface data.
    MPI_Alltoallv(sendSurfaceBuffer, sendCount, sendDispls, MPI_UNSIGNED,
                  recvSurfaceBuffer, recvCount, recvDispls, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
    //  Send and receive particle state data.
    for (int i = 0; i < size; i++)
    {   sendCount[i]  *= 2; recvCount[i]  *= 2;
        sendDispls[i] *= 2; recvDispls[i] *= 2; }
    MPI_Alltoallv(sendStatesBuffer, sendCount, sendDispls, MPI_UNSIGNED,
                  recvStatesBuffer, recvCount, recvDispls, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
    
    recvChangesCount = 0;
    for (int i = 0; i < size; i++)
    {   recvChangesCount += recvCount[i]; }
    
    delete[] numToRecv;
    delete[] sendCount;
    delete[] recvCount;
    delete[] sendDispls;
    delete[] recvDispls; }

/*  resetVariables()
 *
 *  Reset the values of variables unique to each iteration.
 */
void resetVariables()
{   sendChangesBufferSize = 0;
    internalChangesBuffer.clear();
    internalStatesBuffer.clear(); }

/*  placeOnSurfaceA()
 *  
 *  Add a particle to the surface and recursively affect the number of nearest
 *  neighbors accordingly.
 */
void placeOnSurfaceA(Particle& p)
{   if (!hasDissolved(&p) && p.onSurface != onSurfaceA)
    {   surfaceA.push_back(&p);
        p.onSurface = onSurfaceA; }
    
    //  Check the neighboring particles to see if any of them have less than
    //  a full complement of neighbors; if they do, add them to the surface.
    static uint initialN;
    initialN = accumulate(p.initialN.begin(), p.initialN.begin() + dissolnStates, 0);
    for (uint i = 0; i < initialN; i++)
    {   id_t *NList = p.getNList();
        ParticleMap::iterator iter = pmap.find(NList[i]);
        if (iter->second.onSurface != onSurfaceA && iter->second.state < dissolnStates)
        {   static uint countNN;
            countNN = accumulate(iter->second.countN.begin(),
                                 iter->second.countN.begin() + dissolnStates, 0);
            if (countNN < maxNN) placeOnSurfaceA(iter->second); } } }

/*  placeOnSurfaceB()
 *  
 *  Add a particle to the surface and recursively affect the number of nearest
 *  neighbors accordingly.
 */
void placeOnSurfaceB(Particle& p)
{   if (hasDissolved(&p) && p.onSurface != onSurfaceB)
    {   surfaceB.push_back(&p);
        p.onSurface = onSurfaceB; }
    
    //  Check the neighboring particles to see if any of them have less than
    //  a full complement of neighbors; if they do, add them to the surface.
    static uint initialN;
    initialN = accumulate(p.initialN.begin(), p.initialN.begin() + dissolnStates, 0);
    for (uint i = 0; i < initialN; i++)
    {   id_t *NList = p.getNList();
        ParticleMap::iterator iter = pmap.find(NList[i]);
        if (iter->second.onSurface != onSurfaceB && iter->second.state >= dissolnStates)
        {   static uint countNN;
            countNN = accumulate(iter->second.countN.begin(),
                                 iter->second.countN.begin() + dissolnStates, 0);
            if (countNN < maxNN) placeOnSurfaceB(iter->second); } } }

/*  updateLocal()
 *  
 *  Update the state of the neighbors of each local dissolved particle.
 */
void updateLocal()
{   static uint initialN;
    
    for (uint i = 0; i < internalChangesBuffer.size(); i++)
    {   initialN = accumulate(internalChangesBuffer[i]->initialN.begin(),
                              internalChangesBuffer[i]->initialN.begin() + dissolnStates, 0);
        for (uint n = 0; n < initialN; n++)
        {   ParticleMap::iterator iter = pmap.find(internalChangesBuffer[i]->getNList()[n]);
            if (iter != pmap.end())
            {   placeOnSurfaceA(iter->second);
                iter->second.countN[internalStatesBuffer[i]]--;
                iter->second.countN[iter->second.state]++; } } } }

/*  updateExternal()
 *
 *  Update the state of the neighbors of each remote dissolved particle.
 */
void updateExternal()
{   for (uint i = 0; i < recvChangesCount; i++)
    {   ParticleMap::iterator iter = pmap.find(recvChangesBuffer[i]);
        if (iter != pmap.end())
        {   placeOnSurfaceA(iter->second);
            iter->second.countN[recvStatesBuffer[2 * i    ]]--;
            iter->second.countN[recvStatesBuffer[2 * i + 1]]++; } } }

/*  uniformRand()
 *  
 *  Return a uniform random number (uses Boost library uniform_01).
 */
static inline double uniformRand()
{   static boost::mt19937 rng(ranseed);
    static boost::uniform_01<boost::mt19937> zeroOne(rng);
    return zeroOne(); }

/*  calcProbs()
 *  
 *  As there are only maxNN discrete states permissible for any given
 *  temperature, pre-calculate the probabilities for each arrangement.
 *  See for details:  Alekseev YV, Alekseev GY, Bityurin VA. (2002).  Atomic--
 *  Topological and Statistical Background for a Correct Dissolution Theory of
 *  Crystal Substances. Protection of Metals.  38(6):517-529.
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

/*  updateState(Particle* ptr)
 *  
 *  Calculate the probability that the particle indicated by ptr has its state
 *  altered by the currently applicable rule.
 */
void updateState(Particle* ptr)
{   static uint nn;
    nn = accumulate(ptr->countN.begin(), ptr->countN.begin() + dissolnStates, 0);
    if (!hasDissolved(ptr))
    {   /// Dissolution:
        if (uniformRand() < probA[ptr->state][nn])
        {   state_t oldState = ptr->state,
                    newState = rxn[ptr->state].newState;
            ptr->state = newState;
            if (ptr->onBoundary) updateOnBoundary(ptr, oldState, newState);
            else                 updateOffBoundary(ptr, oldState); } }
    else
    {   /// Deposition:
        if (uniformRand() < probB[ptr->state][nn])
        {   state_t oldState = ptr->state,
                    newState = rxn[ptr->state].newState;
            ptr->state = newState;
            if (ptr->onBoundary) updateOnBoundary(ptr, oldState, newState);
            else                 updateOffBoundary(ptr, oldState); } } }

/*  packSurface()
 *  
 *  Remove dissolved or deposited particles from the surface A or B.
 */
void packSurface()
{   static uint numNN;
    
    /// Dissolution:
    for (list<Particle*>::iterator iter = surfaceA.begin(); iter != surfaceA.end(); iter++)
    {   numNN = accumulate((*iter)->countN.begin(), (*iter)->countN.begin() + dissolnStates, 0);
        if (numNN == maxNN)
        {   if ((*iter)->onSurface == onSurfaceA)
            {   (*iter)->onSurface = notOnSurface;
                iter = surfaceA.erase(iter); } }
        else if (hasDissolved(*iter))
        {   placeOnSurfaceB(**iter);
            iter = surfaceA.erase(iter); } }
    
    /// Deposition:
    if (deposition)
    {   for (list<Particle*>::iterator iter = surfaceB.begin(); iter != surfaceB.end(); iter++)
        {   numNN = accumulate((*iter)->countN.begin(), (*iter)->countN.begin() + dissolnStates, 0);
            if (numNN == 0)
            {   if ((*iter)->onSurface == onSurfaceB)
                {   (*iter)->onSurface = notOnSurface;
                    iter = surfaceB.erase(iter); } }
            else if (!hasDissolved(*iter))
            {   placeOnSurfaceA(**iter);
                iter = surfaceB.erase(iter); } } } }

