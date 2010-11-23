/** neighbors.cpp
 *  29 Oct 2009--11 Nov 2010
 *  Neal Davis, Minas Charalambides, and Simon Jenkins
 *  
 *  This file uses MPI buffered send mode, which precludes any other portion
 *  of the program from using this mode simultaneously.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <cstdio>
#include <cstring>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "neighbors.h"
#include "particle.h"
#include "statistics.h" //  Contains fileExists()

#define  SENT_PARTICLE_SIZE 5

using namespace std;

typedef vector<Particle> Cell;

//  External global variables from global.cpp
extern char *dataFileName;
extern bool verbose,
            recalcNN;

extern uint numStates,
            dissolnStates,
            maxNN;

extern vector<Particle>    particles;
extern map<id_t, Particle> pmap;

extern int rank,
           size;

extern coord_t myMinX, myMinY, myMinZ,
               myMaxX, myMaxY, myMaxZ;
extern coord_t xRange, yRange, zRange;

//  Global variables
uint numCellsX, numCellsY, numCellsZ;
vector<vector<vector<Cell> > > cells;

/** void findMaxNN()
 *  
 *  Find the expected number of maximum nearest neighbors based on the geometry.
 *  The algorithm looks first for four mutual neighbors, indicating FCC.  If
 *  that is not found, then three mutual neighbors indicates BCC.  Otherwise,
 *  SC is assumed.
 */
void findMaxNN()
{   maxNN = 0;
    //  Determine global maximum z-coordinate.
    coord_t  maxZ;
    MPI_Allreduce(&myMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if (!rank)
    {   Particle i1, i2, i3;
        bool    twoNN,
                threeNN,
                univTwoNN   = false,
                univThreeNN = false;
        vector<Particle>::iterator  iter,
                                    i1Iter, i2Iter,  i3Iter,
                                    i2Start, i3Start;
        
        //  Cycle through until successful identification of the geometry.
        i1Iter = particles.begin();
        while (!maxNN)
        {   //  Get a base particle for the triangulation.
            i1 = *(i1Iter++);
            if (i1Iter == particles.end()) break;
            if (i1.z < maxZ - 2 * NEIGHBOR_SQUARE_CUTOFF) continue; //  Only look near the top to save processing.
            
            //  Look for two mutual neighbors of the base particle.
            i2Start = i1Iter;   //  Set the starting point for searching for neighbors.
            twoNN = false;
            while (!twoNN)
            {   //  Find a neighbor of i1.
                for (iter = i2Start; iter != particles.end(); iter++)
                {   if (distSqrd(i1, *iter) <= NEIGHBOR_SQUARE_CUTOFF)
                    {   i2          = *iter;
                        i2Iter      = iter;
                        i2Start     = i2Iter;
                        twoNN       = true;
                        univTwoNN   = true;
                        break; } }
                if (iter == particles.end()) break;     //  No neighbor of i1 found.
                
                i3Start = i2Iter + 1;   //  Set the starting point for searching for neighbors.
                threeNN = false;
                while (!threeNN)
                {   //  Find a mutual neighbor of i1 and i2.
                    for (iter = i3Start; iter != particles.end(); iter++)
                    {   if (distSqrd(i1, *iter) <= NEIGHBOR_SQUARE_CUTOFF &&
                            distSqrd(i2, *iter) <= NEIGHBOR_SQUARE_CUTOFF)
                        {   i3          = *iter;
                            i3Iter      = iter;
                            i3Start     = i3Iter;
                            threeNN     = true;
                            univThreeNN = true;
                            break; } }
                    if (iter == particles.end()) break; //  No mutual neighbor of i1 and i2 found.
                    
                    //  Try for a third mutual neighbor; if found, this is FCC.
                    for (iter = i3Iter + 1; iter != particles.end(); iter++)
                    {   if (distSqrd(i1, *iter) <= NEIGHBOR_SQUARE_CUTOFF &&
                            distSqrd(i2, *iter) <= NEIGHBOR_SQUARE_CUTOFF &&
                            distSqrd(i3, *iter) <= NEIGHBOR_SQUARE_CUTOFF)
                        {   maxNN = 12;                     //  FCC
                            break; } } } } }
            if (!maxNN)
            {   if (univThreeNN)  maxNN = 8;      //  BCC
                else              maxNN = 6; } }  //  SC
    
    MPI_Bcast(&maxNN, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD); }

/** void checkBorder()
 *  
 *  Find whether the particle lies within a cutoff radius of the border.
 */
inline void checkBorder(Particle& p, const coord_t y)
{   if(abs(y - p.y) <= NEIGHBOR_SQUARE_CUTOFF) p.onBoundary = true; }

/** double distSqrd()
 *  
 *  Return square of Euclidean distance between two particles.
 */
inline double distSqrd(Particle& p1, Particle& p2)
{   double x, y, z;
    x = p2.x - p1.x;
    y = p2.y - p1.y;
    z = p2.z - p1.z;
    return (x * x + y * y + z * z); }

/** uint getNumCells()
 *  
 *  Find the appropriate number of cells in a given dimension.
 */
inline uint getNumCells(const coord_t range)
{   static double inv3 = 1.0 / (6.0 * NEIGHBOR_SQUARE_CUTOFF);
    uint numCells;
    
    numCells = floor(range * inv3);
    if (numCells <= 0) numCells = 1;
    
    return numCells; }

/*  void createCells()
 *  
 *  Divide the volume allotted to this processor into cells and assign
 *  each particle to the appropriate cell.  The result is held in cells.
 */
void createCells()
{   //  Initialize cells.
    xRange = myMaxX - myMinX;
    yRange = myMaxY - myMinY;
    zRange = myMaxZ - myMinZ;
    numCellsX = getNumCells(xRange);
    numCellsY = getNumCells(yRange);
    numCellsZ = getNumCells(zRange);
    
    cells.resize(numCellsX);
    for (uint i = 0; i < numCellsX; i++)
    {   cells[i].resize(numCellsY);
        for (uint j = 0; j < numCellsY; j++)
        {   cells[i][j].resize(numCellsZ); } }
    
    //  Split particles into cells.
    while(particles.size() > 0)
    {   Particle& p = particles.back();
        coord_t x, y, z;
        uint    k, l, m;
        
        //  Normalize x, y, z into local units.
        x = p.x - myMinX; x /= xRange;
        y = p.y - myMinY; y /= yRange;
        z = p.z - myMinZ; z /= zRange;
        
        //  Determine the cell to which to assign p.
        k = x * numCellsX;  if(k == numCellsX)  k--;
        l = y * numCellsY;  if(l == numCellsY)  l--;
        m = z * numCellsZ;  if(m == numCellsZ)  m--;
        
        //  If p is on the boundary, its information will have to be tagged.
        checkBorder(p, myMinY);
        checkBorder(p, myMaxY);
        
        cells[k][l][m].push_back(p);
        particles.pop_back(); } }

/*  void makeNeighbors(Particle&, Particle&)
 *  
 *  Make two particles have each other's states & ids in their neighbor list.
 */
inline void makeNeighbors(Particle& p1, Particle& p2)
{   if(p1.neighbors.size() > maxNN)
    {   char err[64];
        sprintf(err, "⚠ %d neighbors detected for particle %d at (%f, %f, %f).", p1.neighbors.size(), p1.id, p1.x, p1.y, p1.z);
        warning(err); }
    if(p2.neighbors.size() > maxNN)
    {   char err[64];
        sprintf(err, "⚠ %d neighbors detected for particle %d at (%f, %f, %f).", p2.neighbors.size(), p2.id, p2.x, p2.y, p2.z);
        warning(err); }
    p1.neighbors.push_back(p2.id);
    p1.countN[p2.state]++;
    p2.neighbors.push_back(p1.id);
    p2.countN[p1.state]++; }

/** void makeNeighborsAffectOnlyLeft()
 *  
 *  Make a particle have another's state & id in its neighbor list.
 */
inline void makeNeighborsAffectOnlyLeft(Particle& p1, Particle& p2)
{   if(p1.neighbors.size() > maxNN)
    {   char err[64];
        sprintf(err, "⚠ %d neighbors detected for particle %d at (%f, %f, %f).", p1.neighbors.size(), p1.id, p1.x, p1.y, p1.z);
        warning(err); }
    p1.neighbors.push_back(p2.id);
    p1.countN[p2.state]++; }

/** void findIntraCellNeighbors()
 *  
 *  Test each particle within a cell to see if they are neighbors.
 */
inline void findIntraCellNeighbors(Cell& cell1)
{   for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   for(Cell::iterator iter2 = iter1 + 1; iter2 != cell1.end(); iter2++)
        {   if(distSqrd(*iter1, *iter2) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(*iter1, *iter2); } } } }

/** void findInterCellNeighbors()
 *  
 *  Test each particle within two cells to see if they are neighbors.
 */
inline void findInterCellNeighbors(Cell& cell1, Cell& cell2)
{   for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   for(Cell::iterator iter2 = cell2.begin(); iter2 != cell2.end(); iter2++)
        {   if(distSqrd(*iter1, *iter2) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(*iter1, *iter2); } } } }

/** void findLocalNeighbors()
 *  
 *  Test all particles on this processor to see if they are neighbors.
 */
void findLocalNeighbors()
{   //  Look for new nearest neighbors.
    for(uint i = 0; i < numCellsX; i++)
    {   for(uint j = 0; j < numCellsY; j++)
        {   for(uint k = 0; k < numCellsZ; k++)
            {   //  Find local neighbors.
                findIntraCellNeighbors(cells[i][j][k]);
                
                //  Find adjacent neighbors.
                if(i > 0) findInterCellNeighbors(cells[i][j][k], cells[i - 1][j][k]);
                if(j > 0) findInterCellNeighbors(cells[i][j][k], cells[i][j - 1][k]);
                if(k > 0) findInterCellNeighbors(cells[i][j][k], cells[i][j][k - 1]);
                
                if(i > 1) findInterCellNeighbors(cells[i][j][k], cells[i - 2][j][k]);
                if(j > 1) findInterCellNeighbors(cells[i][j][k], cells[i][j - 2][k]);
                if(k > 1) findInterCellNeighbors(cells[i][j][k], cells[i][j][k - 2]);
                
                //  Find diagonal neighbors.
                if(i > 0 && j > 0) findInterCellNeighbors(cells[i][j][k], cells[i - 1][j - 1][k]);
                if(i > 0 && k > 0) findInterCellNeighbors(cells[i][j][k], cells[i - 1][j][k - 1]);
                if(j > 0 && k > 0) findInterCellNeighbors(cells[i][j][k], cells[i][j - 1][k - 1]); } } } }

/** void getYSlice()
 *  
 *  Pack a slice of a cell into a buffer for messaging.
 */
inline void getYSlice(coord_t* buffer, int& bufferSize, const uint j)
{   bufferSize = 0;
    for (uint i = 0; i < numCellsX; i++)
    {   for (uint k = 0; k < numCellsZ; k++)
        {   for (uint l = 0; l < cells[i][j][k].size(); l++)
            {   buffer[bufferSize++] = (coord_t) cells[i][j][k][l].id;
                buffer[bufferSize++] =           cells[i][j][k][l].x;
                buffer[bufferSize++] =           cells[i][j][k][l].y;
                buffer[bufferSize++] =           cells[i][j][k][l].z;
                buffer[bufferSize++] = (coord_t) cells[i][j][k][l].state; } } } }

/** void checkNeighbors()
 *  
 *  Find out if the sent data describes a particle neighboring any in this cell.
 */
inline void checkNeighbors(Cell& cell1, coord_t* data)
{   Particle other(data[0], data[1], data[2], data[3], data[4]);
    
    for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   if(distSqrd(*iter1, other) <= NEIGHBOR_SQUARE_CUTOFF)
        {   makeNeighborsAffectOnlyLeft(*iter1, other); } } }

/** void findRemoteNeighbors()
 *  
 *  Find neighbors on the borders between processors using MPI buffered send mode.
 */
void findRemoteNeighbors()
{   int         mpi_buffer_size = MAX_ARRAY_SIZE * 2 * sizeof(coord_t);
    byte        mpi_buffer[mpi_buffer_size];
    MPI_Status  status;
    MPI_Request mpr1;
    
    MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);
    
    //  Exchange data from end slices of each processor and check the slices for
    //  nearest-neighbor pairs.
    coord_t buffer[MAX_ARRAY_SIZE];
    int     bufferSize;
    
    if (rank > 0)
    {   getYSlice(buffer, bufferSize, 0);
        MPI_Ibsend(buffer, bufferSize, MPI_DOUBLE, rank - 1, bufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    if (rank < size - 1)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint l = 0; l < (uint) status.MPI_TAG; l++)
        {   for(uint i = 0; i < numCellsX; i++)
            {   for(uint k = 0; k < numCellsZ; k++)
                {   checkNeighbors(cells[i][numCellsY - 1][k], &(buffer[l * SENT_PARTICLE_SIZE])); } } } }
    if (rank < size - 1)
    {   getYSlice(buffer, bufferSize, numCellsY - 1);
        MPI_Ibsend(buffer, bufferSize, MPI_DOUBLE, rank + 1, bufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    if (rank > 0)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint l = 0; l < (uint) status.MPI_TAG; l++)
        {   for(uint i = 0; i < numCellsX; i++)
            {   for(uint k = 0; k < numCellsZ; k++)
                {   checkNeighbors(cells[i][numCellsY - 1][k], &(buffer[l * SENT_PARTICLE_SIZE])); } } } }
    
    MPI_Buffer_detach(mpi_buffer, &mpi_buffer_size); }

/** void fillMap()
 *  
 *  Map each particle into the calculation map so it can be indexed by ID.
 */
inline void fillMap()
{   pmap.clear();
    for(uint x = 0; x < numCellsX; x++)
    {   for(uint y = 0; y < numCellsY; y++)
        {   for(uint z = 0; z < numCellsZ; z++)
            {   for(Cell::iterator iter = cells[x][y][z].begin(); iter != cells[x][y][z].end(); iter++)
                {   Particle& p = *iter;
                    p.onSurface = notOnSurface;
                    pmap.insert(pair<id_t, Particle>(p.id, p)); }
                cells[x][y][z].clear(); } } } }

/** bool loadNeighbors()
 *  
 *  Load the previously-generated list of nearest neighbors.  Note that if there
 *  is a problem with the existence of one of the files but not all, undefined
 *  behavior may result.  Also, the nearest-neighbor files generated are
 *  specific to the number of processors running the program, and must be
 *  recalculated, especially if you previously ran with more processors than
 *  you are now using.
 */
bool loadNeighbors()
{   if (recalcNN) return false;
    
    fstream inDataFile;
    char    inDataFileName[36];
    sprintf(inDataFileName, "%s.P=%d.nn", dataFileName, rank);
    if (!fileExists(inDataFileName))
    {   char err[64];
        if (!rank) sprintf(err, "⚠ Unable to find neighbor files for %s; calculating instead.", dataFileName);
        warning(err);
        return false; }
    
    inDataFile.open(inDataFileName, fstream::in);
    if (!inDataFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to load file %s", inDataFileName);
        error(err); }
    
    try
    {   ParticleMap::iterator mapIter;
        char line[256];
        char *word;
        
        //  Map each particle into the calculation map so it can be indexed by ID.
        for (vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
        {   pmap.insert(pair<id_t, Particle>((*iter).id, *iter)); }
        
        while (inDataFile.getline(line, 256, '\n'))
        {   //  Input the particle id and the ids of its neighbors.
            word = strtok(line, " \t");
            mapIter = pmap.find(atoi(word));
            if (mapIter == pmap.end()) throw 1;
            
            while ((word = strtok(NULL, " \t")) != NULL)
            {   mapIter->second.neighbors.push_back(atoi(word)); }
            
            inDataFile.getline(line, 256, '\n');
            word = strtok(line, " \t");
            uint i = 0;
            do
            {   if (i > numStates) throw 2;
                mapIter->second.countN[i++] = atoi(word); }
            while ((word = strtok(NULL, " \t")) != NULL); }
        
        inDataFile.close(); }
    
    catch (int errType)
    {   char err[64];
        if (errType == 1 && !rank) sprintf(err, "⚠ Output file %s incorrectly formatted (extra newline at end?);\n  calculating neighbors instead.", inDataFileName);
        else if (errType == 2 && !rank) sprintf(err, "⚠ Output file %s incorrectly formatted (too many particle states);\n  calculating neighbors instead.", inDataFileName);
        warning(err);
        return false; }
    
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to load output file %s.", inDataFileName);
        error(err); }
    
    return true; }

/** void outputNeighbors()
 *  
 *  Output the calculated list of nearest neighbors.
 */
void outputNeighbors()
{   fstream outDataFile;
    char    outDataFileName[36];
    sprintf(outDataFileName, "%s.P=%d.nn", dataFileName, rank);
    outDataFile.open(outDataFileName, fstream::out);
    if (!outDataFile)
    {   char err[64];
        sprintf(err, "⚠ Unable to create file %s", outDataFileName);
        error(err); }
    
    try
    {   for (ParticleMap::iterator iter = pmap.begin(); iter != pmap.end(); iter++)
        {   //  Output the particle id and the ids of its neighbors.
            outDataFile << iter->second.id << "\t";
            for (vector<id_t>::iterator iter1 = iter->second.neighbors.begin(); iter1 != iter->second.neighbors.end(); iter1++)
            {   outDataFile << *iter1 << "\t"; }
            outDataFile << endl;
            
            //  Output the states of the nearest neighbors.
            for (vector<uint>::iterator iter2 = iter->second.countN.begin(); iter2 != iter->second.countN.end(); iter2++)
            {   outDataFile << *iter2 << "\t"; }
            outDataFile << endl; }
        outDataFile.close(); }
    
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to write to output file %s.", outDataFileName);
        error(err); } }

/** void calculateNeighbors()
 *  
 *  Calculate the number of nearest neighbors of each particle.
 */
void calculateNeighbors()
{   findMaxNN();            if (verbose && !rank) cout << "  Cell geometry identified as " 
                                                       << (maxNN == 12 ? "face-centered cubic." : (maxNN == 8 ? "body-centered cubic." : "simple cubic."))
                                                       << endl;    cout.flush();
    
    if (!loadNeighbors())
    {   createCells();          if (verbose && !rank) cout << "  Cell partitioning complete." << endl;  cout.flush();
        
        findLocalNeighbors();   if (verbose && !rank) cout << "  Local neighbors found." << endl;       cout.flush();
        findRemoteNeighbors();  if (verbose && !rank) cout << "  Remote neighbors found." << endl;      cout.flush();
        
        fillMap();
        outputNeighbors(); } }

