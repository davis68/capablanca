/** particle.h
 *  29 Oct 2009--14 Jun 2010
 *  Neal Davis
 *  
 */

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <map>
#include <vector>

#include "definitions.h"

using namespace std;

//  External global variables from global.cpp
extern uint numStates,
            maxNN;

//  Global variables
const uint  notOnSurface = 0,
            onSurfaceA = 1,
            onSurfaceB = 2;

struct Particle
{   vector<uint>    countN;     //  Current number of neighbors.  One should increment before decrementing to prevent underflow.
    vector<id_t>    neighbors;  //  List of neighboring particle ids.
    coord_t x, y, z;            //  This particle's coordinates.
    id_t    id;                 //  This particle's id.
    state_t state;              //  Current state of this particle.
    bool    onBoundary;         //  Indicates that the particle may have neighbors on a neighboring processor.
    uint    onSurface;          //  Indicate which surface the particle is on, if any.
    
    Particle()
     : x(0), y(0), z(0), id(0), state(0), onBoundary(false)
    {   countN.resize(numStates, 0); }
    
    Particle(id_t _id, coord_t _x, coord_t _y, coord_t _z, state_t _st)
     : onBoundary(false)
    {   x = _x; y = _y; z = _z; id = _id; state = _st;
        countN.resize(numStates, 0); } };

typedef map<id_t, Particle>     ParticleMap;

#endif  /* _PARTICLE_H */

