#ifndef _PROCESS_H
#define	_PROCESS_H

#include "definitions.h"

inline bool hasDissolved(Particle* ptr);

void findSurface();
void updateOnBoundary(Particle* ptr, state_t oldState, state_t newState);
void updateOffBoundary(Particle* ptr, state_t oldState);

void recvChanges(const int src);
void recvStates(const int src);
void recvSurface(const int src);
void sendChanges(const int dest);
void sendStates(const int dest);
void sendSurface(const int dest);
void exchangeInterNodeChanges();

void resetVariables();
void placeOnSurfaceA(Particle& p);
void placeOnSurfaceB(Particle& p);
void updateLocal();
void updateExternal();

void calcProbs();
void updateState(Particle* ptr);
void packSurface();

void process();

#endif	/* _PROCESS_H */

