/** input.cpp
 *  27 Feb 2010--28 Jul 2010
 *  Neal Davis
 */

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "definitions.h"
#include "error.h"
#include "input.h"

using namespace std;

//  External global variables from global.cpp
extern char    *progName,
               *progVers,
               *progCL,
               *confFileName,
               *ruleFileName,
               *dataFileName;
extern bool     verbose;

extern uint     numStates,
                dissolnStates,
                maxNN;
extern double   T,
                phi;
extern Reaction *rxnA,
                *rxnB;

extern int      outputInterval,
                tmax,
                ranseed;
extern bool     deposition,
                ranseedspec,
                recalcNN;

extern int      rank;

void outputHelp()
{   cout << "capablanca is a state transition machine for simulating electrochemical" << endl
         << "reactions, especially as a function of nearest neighbors.  It is" << endl
         << "designed for use on a cluster supporting MPI, and it is recommended" << endl
         << "that the terminal support Unicode UTF-8 or higher." << endl << endl
         << "Usage:  " << progName << " [OPTION]..." << endl << endl
         << "  -c CONFIG.CFG     specify configuration file;" << endl
         << "  -d                suppress deposition;" << endl
         << "  -h                display this help message;" << endl
         << "  -E φ              specify potential φ;" << endl
         << "  -i Τ              specify T steps between output files;" << endl
         << "  -n                force recalculation of nearest neighbor files;" << endl
         << "  -p DATAFILE.XYZ   specify particle position file;" << endl
         << "  -R s              specify random seed s;" << endl
         << "  -r RULE_SET.RS    specify rule set file;" << endl
         << "  -T ϑ              specify temperature ϑ;" << endl
         << "  -t τ              specify τ time steps to run;" << endl
         << "  -V                use verbose output;" << endl
         << "  -v                display version information." << endl; }

void outputVersion()
{   cout << progName << " v" << progVers << endl; }

/*  loadRules()
 *  Loads the ruleset file specified in ruleFileName.
 */
void loadRules()
{   ifstream    ruleFile;
    ruleFile.open(ruleFileName, ifstream::in);
    if (!ruleFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", ruleFileName);
        error(err); }
    
    try
    {   //  Set number of rules and states.
        uint    numRules;
        ruleFile >> numRules >> numStates >> dissolnStates;
        rxnA = new Reaction[numStates];
        rxnB = new Reaction[numStates];
        
        //  Load parameters for each rule.
        char    tempStr[256];
        uint    state;
        double  tmp;
        for (uint i = 0; i < numRules; i++)
        {   //  Read and ignore comment explaining reaction.
            ruleFile.getline(tempStr, 256);
            ruleFile.getline(tempStr, 256);
            
            //  Input the appropriate parameters for this reaction.
            ruleFile >> state;
            ruleFile >> rxnA[state].prefactor >> tmp
                     >> rxnA[state].alpha >> rxnA[state].E_s
                     >> rxnA[state].E_r >> rxnA[state].z >> rxnA[state].newState;
            rxnB[rxnA[state].newState].prefactor = tmp;
            rxnB[rxnA[state].newState].alpha     = 1.0 - rxnA[state].alpha; //  This acts as beta.
            rxnB[rxnA[state].newState].E_s       = rxnA[state].E_s;
            rxnB[rxnA[state].newState].E_r       = rxnA[state].E_r;
            rxnB[rxnA[state].newState].z         = rxnA[state].z;
            rxnB[rxnA[state].newState].newState  = state; }
        ruleFile.close(); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Rule set file %s incorrectly formatted.", ruleFileName);
        error(err); } }

/*  loadConfig()
 *  Loads the configuration file specified in confFileName.
 */
 void loadConfig()
{   ifstream    confFile;
    confFile.open(confFileName, ifstream::in);
    if (!confFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", confFileName);
        error(err); }
    
    try
    {   //  Get temperature and other parameters from input file.
        double  tmpT;
        confFile >> tmpT;
        if (T == 0.0)   T = tmpT;
        
        double  tmpphi;
        confFile >> tmpphi;
        if (phi == 0.0) phi = tmpphi;
        
        uint    tmpOI;
        confFile >> tmpOI;
        if (outputInterval == 0)  outputInterval = tmpOI;
        
        uint    tmptmax;
        confFile >> tmptmax;
        if (tmax == 0)    tmax = tmptmax;
        
        confFile.close(); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Configuration file %s incorrectly formatted.", confFileName);
        error(err); } }

/*  parseInput(const int argc, char** argv)
 *  
 *  Interpret the command line input.
 */
void parseInput(const int argc, char** argv)
{   confFileName= new char[32];
    ruleFileName= new char[32];
    dataFileName= new char[32];
    verbose     = false;
    deposition  = true;
    ranseedspec = false;
    recalcNN    = false;
    
    bool    confFlag = false,
            ruleFlag = false,
            dataFlag = false;
    
    progName = new char[12];
    progVers = new char[8];
    progCL   = new char[256];
    strcpy(progName, "capablanca");
    strcpy(progVers, PROGRAM_VERSION);
    strcpy(progCL, argv[0]);
    for (int i = 1; i < argc; i++)
    {   strcat(progCL, " ");
        strcat(progCL, argv[i]); }
    
    //  Set variables which will be affected by the configuration file to zero
    //  as a flag that they should be overridden if they are still zero.
    T   = 0.0;
    phi = 0.0;
    outputInterval = 0;
    tmax= 0;
    
    char   *val;
    if (argc > 1)
    {   int i = 1;
        val = strtok(argv[i++], " ,.-");
        while (val != NULL)
        {   switch (val[0])
            {   case 'c':   //  Configuration file specified.
                    val = strtok(argv[i++], " ,");
                    confFileName = val;
                    confFlag = true;
                    break;
                    
                case 'd':   //  Suppress deposition.
                    deposition = false;
                    break;
                    
                case 'E':   //  Potential specified.
                    val = strtok(argv[i++], " ,-");
                    phi = atof(val);
                    break;
                    
                case 'h':   //  Help.
                    if (!rank) outputHelp();
                    cout.flush();
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                    break;
                    
                case 'i':   //  Output interval specified.
                    val = strtok(argv[i++], " ,.-");
                    outputInterval = atoi(val);
                    break;
                    
                case 'n':   //  Force recalculation of nearest neighbor list.
                    recalcNN = true;
                    break;
                    
                case 'p':   //  Particle position file specified.
                    val = strtok(argv[i++], " ,");
                    dataFileName = val;
                    dataFlag = true;
                    break;
                    
                case 'R':   //  Random seed specified.
                    val = strtok(argv[i++], " ,");
                    ranseed = atoi(val);
                    ranseedspec = true;
                    break;
                    
                case 'r':   //  Rule set file specified.
                    val = strtok(argv[i++], " ,");
                    ruleFileName = val;
                    ruleFlag = true;
                    break;
                    
                case 'T':   //  Temperature specified.
                    val = strtok(argv[i++], " ,-");
                    T   = atof(val);
                    break;
                    
                case 't':   //  Number of time steps specified.
                    val = strtok(argv[i++], " ,.-");
                    tmax= atoi(val);
                    break;
                    
                case 'V':   //  Verbose output desired.
                    verbose = true;
                    break;
                    
                case 'v':   //  Version.
                    if (!rank) outputVersion();
                    cout.flush();
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                    break;
                    
                default:    //  Unknown parameter encountered.
                    char err[64];
                    sprintf(err, "⚠ Unknown input parameter %s", val);
                    if (!rank) outputHelp();
                    cout.flush();
                    error(err);
                    break; }
            val = strtok(argv[i++], " ,.-"); } }
    
    if (!confFlag)  strcpy(confFileName, CONF_FILENAME);
    if (!ruleFlag)  strcpy(ruleFileName, RULE_FILENAME);
    if (!dataFlag)  strcpy(dataFileName, DATA_FILENAME); }

