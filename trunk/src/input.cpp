#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "definitions.h"
#include "error.h"
#include "input.h"
#include "mpi.h"

using namespace std;

extern char    *progName,
               *progVers,
               *confFileName,
               *ruleFileName,
               *dataFileName;
extern bool     verbose;

extern uint     numStates,
                dissolnStates,
                maxNN;
extern double   T,
                phi;
extern Reaction *rxn;

extern int      outputInterval,
                tmax;

extern int      rank;

void outputHelp()
{   cout << "capablanca is a state transition machine for simulating electrochemical" << endl
         << "reactions, especially as a function of nearest neighbors.  It is" << endl
         << "designed for use on a cluster supporting MPI, and it is recommended" << endl
         << "that the terminal support Unicode UTF-8 or higher." << endl << endl
         << "Usage:  " << progName << " [OPTION]..." << endl << endl
         << "  -c CONFIG.CFG     specify configuration file;" << endl
         << "  -h                display this help message;" << endl
         << "  -E φ              specify potential φ;" << endl
         << "  -i Τ              specify Τ steps between output files;" << endl
         << "  -p DATAFILE.XYZ   specify particle position file;" << endl
         << "  -r RULE_SET.RS    specify rule set file;" << endl
         << "  -T ϑ              specify temperature ϑ;" << endl
         << "  -t τ              specify τ time steps to run;" << endl
         << "  -V                use verbose output;" << endl
         << "  -v                display version information." << endl;
}

void outputVersion()
{   cout << progName << " v" << progVers << endl; }

/*  loadRules()
 *  Loads the ruleset file specified in ruleFileName.
 */
void loadRules()
{   uint        numRules;
    
    ifstream    ruleFile;
    ruleFile.open(ruleFileName, ifstream::in);
    if (!ruleFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", ruleFileName);
        error(err); }
    
    //  Set number of rules and states.
    ruleFile >> numRules >> numStates >> dissolnStates;
    rxn = new Reaction[numStates];
    
    //  Load parameters for each rule.
    char    tempStr[32];
    uint    state;
    for (uint i = 0; i < numRules; i++)
    {   //  Read and ignore comment explaining reaction.
        ruleFile >> tempStr;
        //  Input the appropriate parameters for this reaction.
        ruleFile >> state;
        ruleFile >> rxn[state].prefactor >> rxn[state].alpha >> rxn[state].E_s
                 >> rxn[state].E_r >> rxn[state].z >> rxn[state].newState;
        rxn[state].beta = 1 - rxn[state].alpha; }
    ruleFile.close(); }

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
    
    //  Get temperature and other parameters from input file.
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
    
    confFile.close();
}

void parseInput(const int argc, char** argv)
{   confFileName= new char[32];
    ruleFileName= new char[32];
    dataFileName= new char[32];
    verbose     = false;
    
    bool    confFlag = false,
            ruleFlag = false,
            dataFlag = false;
    
    progName = new char[32];
    progVers = new char[8];
    strcpy(progName, argv[0]);
    strcpy(progVers, PROGRAM_VERSION);
    
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
                    
                case 'p':   //  Particle position file specified.
                    val = strtok(argv[i++], " ,");
                    dataFileName = val;
                    dataFlag = true;
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
                    sprintf(err, "⚠ Unknown input parameter %d", val[0]);
                    error(err);
                    break;
            }
            val = strtok(argv[i++], " ,.-");
        }
    }
    
    if (!confFlag)  strcpy(confFileName, CONF_FILENAME);
    if (!ruleFlag)  strcpy(ruleFileName, RULE_FILENAME);
    if (!dataFlag)  strcpy(dataFileName, DATA_FILENAME);
}
