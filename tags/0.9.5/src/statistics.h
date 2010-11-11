#ifndef _STATISTICS_H
#define _STATISTICS_H

#include <list>
#include "definitions.h"

bool fileExists(char*);
void outputSurface(list<Particle*>, uint);
void outputToFile(uint, uint);
void collateStatistics(uint);
void collectStatFiles();

#endif  /* _STATISTICS_H */

