#include "definitions.h"
#include "particle.h"

//  Input file names and initializers.
char   *progName,       //  Init in input.cpp parseInput()
       *progVers,       //  Init in input.cpp parseInput()
       *confFileName,   //  Init in input.cpp parseInput()
       *ruleFileName,   //  Init in input.cpp parseInput()
       *dataFileName;   //  Init in input.cpp parseInput()
bool    verbose;        //  Init in input.cpp parseInput()

//  Reaction and state variables.
uint        numStates,      //  Init in input.cpp loadRules()
            dissolnStates,  //  Init in input.cpp loadRules()
            maxNN;          //  Init in input.cpp loadRules()
double      T,          //  Init in input.cpp parseInput() or loadConfig()
            phi;        //  Init in input.cpp parseInput() or loadConfig()
Reaction   *rxn;        //  Init in input.cpp loadRules()

//  Simulation parameters.
uint    outputInterval, //  Init in input.cpp parseInput() or loadConfig()
        tmax;           //  Init in input.cpp parseInput() or loadConfig()

//  MPI and OpenMP variables.
int     rank,           //  Init in main.cc main()
        size;           //  Init in main.cc main()

//  Calculation variables.
vector<vector<double> > prob;       //  Init in process.cc calcProbs()
ParticleList            particles;  //  Init in ***
ParticleMap             pmap;       //  Init in ***
vector<ParticlePtr>     surface;    //  Init in ***

//  Nearest-neighbor calculation variables.
coord_t myMinX = 0.0, myMinY = 0.0, myMinZ = 0.0,
        myMaxX = 0.0, myMaxY = 0.0, myMaxZ = 0.0,
        minY, maxY;

//  Output & statistics variables.***
vector<uint>    oldParticleCount;   //  Init in process.cc process()
vector<uint>    myParticleCount;    //  Init in process.cc process()
vector<double>  rates;              //  Init in process.cc process()
uint            initialTotalParticles;

