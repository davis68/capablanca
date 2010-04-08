#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "definitions.h"
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

using namespace std;

extern uint numStates,
            dissolnStates,
            maxNN;

struct Particle
{   vector<uint>    countN,     //  Current number of neighbors.
                    initialN;   //  Initial number of neighbors.
    vector<vector<id_t> >   neighbors;  //  List of neighboring particle ids.
    id_t    id;         //  This particle's id.
    state_t state;      //  Current state of this particle.
    bool    onBoundary, //  Indicates that the particle may have neighbors on a neighboring processor.
            onSurface;  //  Indicates whether the particle is already on the surface.
    coord_t x, y, z;    //  This particle's coordinates.
    
    Particle()
     : onBoundary(false), id(0), x(0), y(0), z(0), state(0)
    {   countN.resize(numStates, 0);
        initialN.resize(numStates, 0);
        neighbors.resize(numStates);
        for (uint i = 0; i < neighbors.size(); i++)
        {   neighbors[i].resize(0); } }
    Particle(id_t _id, coord_t _x, coord_t _y, coord_t _z, state_t _st)
     : onBoundary(false)
    {   id = _id; x = _x; y = _y; z = _z; state = _st;
        countN.resize(numStates, 0);
        initialN.resize(numStates, 0);
        neighbors.resize(numStates);
        for (uint i = 0; i < neighbors.size(); i++)
        {   neighbors[i].resize(0); } }
    
    void output()
    {   cout << id << " [" << state << "] @ ("
             << x << ", " << y << ", " << z << ") with "
             << accumulate(countN.begin(), countN.begin() + dissolnStates, 0) << endl;
        cout.flush(); }
    void outputN()
    {   cout << id << ": ";
        for (uint i = 0; i < neighbors.size(); i++)
        {   for (uint j = 0; j < neighbors[i].size(); j++)
            {   cout << neighbors[i][j] << "; "; } }
        cout << endl; cout.flush(); }
    
    /*  id_t* getNList()
     *  
     *  Returns an array of the ids of the neighboring particles.
     */
    id_t* getNList()
    {   id_t *nnList = new id_t[maxNN];
        uint k = 0;
        
        for (uint i = 0; i < neighbors.size(); i++)
        {   for (uint j = 0; j < neighbors[i].size(); j++)
            {   nnList[k] = neighbors[i][j]; } }
        return nnList; }
};

typedef Particle*               ParticlePtr;
typedef vector<Particle>        ParticleList; //*** perhaps a better name?
typedef map<id_t, Particle>     ParticleMap;
typedef ParticleMap::iterator   pmapiter;     //*** perhaps a better name?

#endif  /* _PARTICLE_H */

