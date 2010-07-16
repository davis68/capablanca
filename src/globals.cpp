/** globals.cpp
 *  29 Oct 2009--16 Jun 2010
 */

#include <list>
#include <vector>

#include "definitions.h"
#include "particle.h"

//  MPI variables.
int     rank,           //  Init in main.cpp main()
        size;           //  Init in main.cpp main()

//  Input file names and initializers.
char   *progName,       //  Init in input.cpp parseInput()
       *progVers,       //  Init in input.cpp parseInput()
       *confFileName,   //  Init in input.cpp parseInput()
       *ruleFileName,   //  Init in input.cpp parseInput()
       *dataFileName;   //  Init in input.cpp parseInput()
bool    verbose;        //  Init in input.cpp parseInput()

//  Reaction and state variables.
uint    numStates,      //  Init in input.cpp loadRules()
        dissolnStates,  //  Init in input.cpp loadRules()
        maxNN;          //  Init in input.cpp loadRules()
double  T,              //  Init in input.cpp parseInput() or loadConfig()
        phi;            //  Init in input.cpp parseInput() or loadConfig()
Reaction   *rxnA,       //  Init in input.cpp loadRules() | Freed in main.cpp main()
           *rxnB;       //  Init in input.cpp loadRules() | Freed in main.cpp main()

//  Simulation parameters.
uint    outputInterval, //  Init in input.cpp parseInput() or loadConfig()
        tmax,           //  Init in input.cpp parseInput() or loadConfig()
        ranseed;        //  Init in process.cpp process()
bool    cyclical,       //  Init in input.cpp parseInput()
        deposition;     //  Init in input.cpp parseInput()

//  Calculation variables.
vector<vector<double> > probA;      //  Init in process.cpp calcProbs()
vector<vector<double> > probB;      //  Init in process.cpp calcProbs()
vector<Particle>        particles;  //  Init in readxyz.cpp readXYZ()
ParticleMap             pmap;       //  Init in neighbors.cpp fillMap()

list<Particle*>         surfaceA;   //  Init in process.cpp process()
list<Particle*>         surfaceB;   //  Init in process.cpp process()

//  Nearest-neighbor calculation variables.
coord_t myMinX =  1e8, myMinY =  1e8, myMinZ =  1e8,
        myMaxX = -1e8, myMaxY = -1e8, myMaxZ = -1e8,
        minY, maxY;

//  Output & statistics variables.
uint   *oldParticleCount;       //  Init in process.cpp process() | Freed in process.cpp process()
uint   *newParticleCount;       //  Init in process.cpp process() | Freed in process.cpp process()
uint   *myParticleCount;        //  Init in process.cpp process() | Freed in process.cpp process()
Vector  rates;                  //  Init in process.cpp process()
uint    initialTotalParticles;  //  Init in readxyz.cpp readXYZ()

