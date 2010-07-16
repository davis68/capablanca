/** process.h
 *  29 Oct 2009--14 Jun 2010
 *  Minas Charalambides and Neal Davis
 *  
 */

#ifndef _PROCESS_H
#define	_PROCESS_H

#include <boost/random.hpp>             //  These may require explicit inclusion of the Boost C++
#include <boost/random/uniform_01.hpp>  //  library location when compiling.

#include "definitions.h"

//  External global variables from global.cpp
extern uint ranseed;

//  Function prototypes
bool hasDissolved(Particle*);

void findSurface();
void updateBoundaryParticle(Particle*, state_t, state_t);
void updateInternalParticle(Particle*, state_t, state_t);
void packSurface();

void updateLocalNeighbors();
void exchangeInterNodeChanges();
void updateExternalNeighbors();

void placeOnSurfaceA(Particle&);
void placeOnSurfaceB(Particle&);
void expandSurfaces(Particle&);
void expandSurfaceA(Particle&);
void expandSurfaceB(Particle&);
void removeFromSurfaceA(Particle);
void removeFromSurfaceB(Particle);

void calcProbs();
void resetVariables();
void transitionParticle(Particle*);
void setParticleStates();

void process();

/*  uniformRand()
 *  
 *  Return a uniform random number (uses Boost library uniform_01).
 */
static inline double uniformRand()
{   static boost::mt19937 rng(ranseed);
    static boost::uniform_01<boost::mt19937> zeroOne(rng);
    return zeroOne(); }

#endif	/* _PROCESS_H */

