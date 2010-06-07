/** particle.h
 *  29 Oct 2009--07 Jun 2010
 *  Neal Davis
 */

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "definitions.h"
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

using namespace std;

extern uint numStates,
            maxNN;

//FIXME:enum surface {notOnSurface, onSurfaceA, onSurfaceB};
const uint  notOnSurface = 0,           //  FIXME enum this?
            onSurfaceA = 1,
            onSurfaceB = 2;

struct Particle
{   vector<uint>    countN,     //  Current number of neighbors.
                    initialN;   //  Initial number of neighbors.
    vector<id_t>    neighbors;  //  List of neighboring particle ids.
    coord_t x, y, z;        //  This particle's coordinates.
    id_t    id;             //  This particle's id.
    state_t state;          //  Current state of this particle.
    bool    onBoundary;     //  Indicates that the particle may have neighbors on a neighboring processor.
    uint    onSurface;      //  Indicate whether the particle is already on a surface.
    
    Particle()
     : x(0), y(0), z(0), id(0), state(0), onBoundary(false)
    {   countN.resize(numStates, 0);
        initialN.resize(numStates, 0); }
    Particle(id_t _id, coord_t _x, coord_t _y, coord_t _z, state_t _st)
     : onBoundary(false)
    {   x = _x; y = _y; z = _z; id = _id; state = _st;
        countN.resize(numStates, 0);
        initialN.resize(numStates, 0); }
    
    /*  id_t* getNList()
     *  
     *  Returns an array of the ids of the neighboring particles.
     */
    id_t* getNList()
    {   id_t *nnList = new id_t[maxNN];
        uint k = 0;
        
        for (uint i = 0; i < neighbors.size(); i++)
        {   nnList[k++] = neighbors[i]; }
        
        return nnList; } };

typedef map<id_t, Particle>     ParticleMap;

#endif  /* _PARTICLE_H */

