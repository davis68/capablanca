#ifndef _DEFINITIONS_H
#define _DEFINITIONS_H

typedef unsigned int    uint;
typedef unsigned long   ulong;
typedef unsigned short  ushort;
typedef unsigned char   byte;

typedef uint    id_t;                   // the type used for particle ids
typedef double  coord_t;                // the type used for particle coordinates
typedef uint    state_t;                // the type used for a particle's current state

//  Input defaults.
#define PROGRAM_VERSION "0.5.0"         // version of this code
#define CONF_FILENAME   "standard.conf" // the file into which data is generated
#define RULE_FILENAME   "rules.rs"      // the file into which data is generated
#define DATA_FILENAME   "particles.xyz" // the file into which data is generated

//  Calculation.
#define MAX_ARRAY_SIZE      1500000     // the maximum size of an array to be communicated between two CPUs ***nonoptimal?
#define MAX_NEIGHBORS       12          // the maximum number of neighbors a particle can have
#define NUM_OF_THREADS      2           // the number of threads per node
#define DISSOLUTION_STATE   6           // the state which represents complete dissolution

//  Nearest-neighbor finding.
#define BORDER_CUTOFF   1.0
#define NEIGHBOR_SQUARE_CUTOFF  1.004   // determined by precision in XYZ file, 0.001 at best.
#define NUM_CELLS       128

//  Processing.
#define SURFACE_CUTOFF  1.0             // particles within that distance from the top are considered to be in the surface

//  Debugging.
//#define PROC_DEBUG
//#define STAT_DEBUG

const double    k_B = 1.0,
                ec = 1.0;

struct Reaction
{   double  prefactor,
            E_s,
            E_r,
            alpha,
            beta,
            z;
    state_t newState;
};

#endif  /* _DEFINITIONS_H */

