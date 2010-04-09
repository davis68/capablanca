//  Another thing I'd like to fix is the logick of this file and process_helpers.h and process.h.
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "particle.h"
#include "process.h"
#include "process_helpers.h"
#include "statistics.h"

using namespace std;

extern bool verbose;

extern uint         numStates,
                    maxNN;
extern double       T,
                    phi;
extern Reaction    *rxn;

extern int          outputInterval,
                    tmax;

extern vector<vector<double> >  prob;
extern ParticleMap              pmap;
extern vector<ParticlePtr>      surface;

extern vector<uint>     oldParticleCount;
extern vector<uint>     myParticleCount;
extern vector<double>   rates;

/*  myrand()
 *  
 *  Return a uniform random number.
 */
static inline double myrand()
{   return (double) rand() / (double) RAND_MAX; }

/*  updateState(ParticlePtr ptr)
 *  
 *  Calculate the probability that the particle indicated by ptr has its state
 *  altered by the currently applicable rule.
 */
inline void updateState(ParticlePtr ptr)
{   if (myrand() < prob[ptr->state][accumulate(ptr->countN.begin(), ptr->countN.begin() + dissolnStates, 0)])
    {   //***cerr << rank << "->" << ptr->id << " from " << ptr->state << " to " << rxn[ptr->state].newState << endl;
        ptr->state = rxn[ptr->state].newState;
        if (ptr->onBoundary) updateOnBoundary(ptr);
        else                 updateOffBoundary(ptr); } }

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
    
    for (uint i = 0; i < surface.size(); i++)
    {   if (hasDissolved(surface[i])) // && accumulate(surface[i]->countN.begin(), surface[i]->countN.begin() + dissolnStates, 0) < 1
        {   count++;
            if (i == surface.size() - 1)
            {   surface.pop_back(); }
            else
            {   surface[i] = surface.back();
                surface.pop_back(); } } } }

/*  process()
 *  
 *  Cycle through the iterations of the system, updating each particle both
 *  locally and remotely.
 */
void process()
{   calcProbs();
    findSurface();
    
    oldParticleCount.resize(numStates, 0);
    myParticleCount.resize(numStates, 0);
    rates.resize(numStates, 0);
    
    for (uint t = 0; t < tmax; t++)
    {   resetVariables();
        
        for (vector<ParticlePtr>::iterator iter = surface.begin(); iter != surface.end(); iter++)
        {   updateState(*iter); }
        
        exchangeInterNodeChanges();
        updateLocal();
        updateExternal();
        packSurface();
        
        if (!(t % outputInterval)) collateStatistics(t);
        MPI_Barrier(MPI_COMM_WORLD);    //  Unfortunately, we must remain synchronous.
        
        if (verbose && !(t % outputInterval) && !rank)
        {   cout << "  " << 100 * (double) t / (double) tmax << "% complete.\n";
            cout.flush(); } }
    
    if (!rank)
    {   cerr << "  Collecting & finalizing data output.\n";
        collectStatFiles(); } }

