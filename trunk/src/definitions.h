/** process.h
 *  29 Oct 2009--15 Jul 2010
 */
#ifndef _DEFINITIONS_H
#define _DEFINITIONS_H

#include <vector>

using namespace std;

typedef unsigned int    uint;
typedef unsigned char   byte;

typedef uint            id_t;           //  the type used for particle ids
typedef double          coord_t;        //  the type used for particle coordinates
typedef uint            state_t;        //  the type used for a particle's current state

typedef vector<double>  Vector;

//  Input defaults.
#define PROGRAM_VERSION "0.8.0"         //  version of this code
#define CONF_FILENAME   "standard.conf" //  the default configuration file name
#define RULE_FILENAME   "rules.rs"      //  the default rule set file name
#define DATA_FILENAME   "particles.xyz" //  the default file name of generated data

//  Calculation.
#define MAX_ARRAY_SIZE  2000000         //  the maximum size of an array to be communicated between two CPUs FIXME

//  Nearest-neighbor detection.
#define NEIGHBOR_SQUARE_CUTOFF  1.100   //  determined by precision in XYZ file

//  System parameters.
const double    k_B = 1.0,              //  Boltzmann constant in reduced units
                ec  = 1.0;              //  electron charge in reduced units

struct Reaction
{   double  prefactor,
            E_s,
            E_r,
            alpha,
            z;
    state_t newState; };

#endif  /* _DEFINITIONS_H */

