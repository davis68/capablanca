#ifndef _NEIGHBORS_H
#define	_NEIGHBORS_H

#include "particle.h"

void findMaxNN();
void checkBorder(Particle& p, const coord_t y);
double distSqrd(Particle& p1, Particle& p2);
double distSqrdCyc(Particle& p1, Particle& p2);

uint getNumCells(const coord_t range);
void createCells();

bool alreadyNeighbors(Particle&, Particle&);
void makeNeighbors(Particle& p1, Particle& p2);
void makeNeighborsAffectOnlyLeft(Particle& p1, Particle& p2);
void findLocalNeighbors();
void findRemoteNeighbors();

void fillMap();

bool loadNeighbors();
void outputNeighbors();

void calculateNeighbors();

#endif	/* _NEIGHBORS_H */

