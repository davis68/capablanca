#ifndef _READXYZ_H
#define	_READXYZ_H

#include <vector>
#include "definitions.h"
#include "error.h"

unsigned int getProc(coord_t y);
void dump(std::vector<double> &list, int bucket);
void readXYZ();

#endif	/* _READXYZ_H */

