//  Another thing I'd like to fix is the logick of this file and process_helpers.h and process.h.
//  TODO:  work on validating process()
#include <math.h>
#include <iostream>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "particle.h"
#include "process_helpers.h"
#include "statistics.h"

using namespace std;

extern bool verbose;

//extern uint         receiveCount; //*** extern?
extern ParticleMap  pmap;

extern uint         numStates,
                    maxNN;
extern double       T,
                    phi;
extern Reaction    *rxn;

extern int          outputInterval,
                    tmax;

extern vector<vector<double> >  prob;

extern uint        *oldCount;
extern uint        *count;
extern double      *rates;

/*  resetVariables()
 *
 *  Reset the values of variables unique to each iteration.
 */
inline void resetVariables()
{   interNodeChangesBufferSize = 0;
    receiveCount = 0;
    internalChangesBuffer.clear(); }

/*  myrand()
 *  
 *  Return a uniform random number. *** requires validation of ergodicity
 */
static inline double myrand()
{   return (double) rand() / (double) RAND_MAX; }

/*  updateState(ParticlePtr ptr)
 *  
 *  Calculate the probability that the particle indicated by ptr has its state
 *  altered by the currently applicable rule. *** this needs to take into account more possible combinations
 */
inline void updateState(ParticlePtr ptr)
{   if(myrand() < prob[ptr->state][accumulate(ptr->countN.begin(), ptr->countN.begin() + dissolnStates, 0)])
    {   //cout << rank << ": " << ptr->id << " from " << ptr->state << " to " << rxn[ptr->state].newState << endl; //***
        ptr->state = rxn[ptr->state].newState; } }

/*  calcProbs()
 *  
 *  As there are only maxNN discrete states permissible for any given
 *  temperature, pre-calculate the probabilities for each arrangement.
 */
inline void calcProbs()
{   //  Prepare array of probabilities indexed by state and nn.
    prob.resize(numStates);
    for (uint i = 0; i < prob.size(); i++)
    {   prob[i].resize(maxNN); }
    
    //  Calculate values for each state and nn condition.
    for (uint i = 0; i < numStates; i++)
    {   for (uint j = 0; j < maxNN; j++)
        {   prob[i][j] = rxn[i].prefactor
                       * exp(-rxn[i].alpha * j * rxn[i].E_s      / (k_B * T))
                       * exp(-rxn[i].alpha     * rxn[i].E_r      / (k_B * T))
                       * exp(-rxn[i].alpha * rxn[i].z * ec * phi / (k_B * T)); } } }

/*  packSurface()
 *  
 *  Remove dissolved particles from the surface. ***
 */
void packSurface()
{   uint count = 0;
    
    for(uint i = 0; i < surface.size(); i++)
    {   if(hasDissolved(surface[i])) // && accumulate(surface[i]->countN.begin(), surface[i]->countN.begin() + dissolnStates, 0) < 1
        {   count++;
            if(i == surface.size() - 1)
            {   surface.pop_back(); }
            else
            {   surface[i] = surface.back();
                surface.pop_back(); }
        }
    }
    //***
    if (verbose && count > 0) cout << "rank: " << rank << " packed " << count
                                   << " particles" << endl;
}

/*  process()
 *  
 *  Cycle through the iterations of the system, updating each particle both
 *  locally and remotely.
 */
void process()
{   calcProbs();
    
    findSurface();
    if (verbose) cout << rank << ":  Surface found with " << surface.size()
                      << " elements." << endl;
    
    for(uint iter = 0; iter < tmax; iter++)
    {   resetVariables();
        
        for(uint i = 0; i < surface.size(); i++)
        {   ParticlePtr ptr = surface[i];
            updateState(ptr);
            if (hasDissolved(ptr))
            {   cerr << rank << ": " << ptr->id << " dissolved;" << endl; //***
                if(ptr->onBoundary) processOnBoundaryDissolved(ptr);
                else processNotOnBoundaryDissolved(ptr); } }
        
        //***
        #ifdef PROC_DEBUG
        cerr<<"rank: "<<rank<<" now sending/receiving"<<endl;
        #endif
        
        exchangeInterNodeChanges();
        
        #ifdef PROC_DEBUG
        cerr<<"rank: "<<rank<<" now updating local"<<endl;
        #endif
        
        updateLocal();
        
        #ifdef PROC_DEBUG
        cerr<<"rank: "<<rank<<" now updating external"<<endl;
        #endif
        
        updateExternal();
        
        #ifdef PROC_DEBUG
        cerr<<"rank: "<<rank<<" now packing surface"<<endl;
        #endif
        
        packSurface();
        
        //if (!(iter % outputInterval)) collateStatistics(iter);
        
        MPI_Barrier(MPI_COMM_WORLD);    //  Unfortunately, we must remain synchronous.
        /*if (verbose && !rank)
        {   cout << 100 * (double) iter / (double) tmax << "% complete." << endl;
            cout.flush(); }*/ }
    
    //collectStatFiles();
}

