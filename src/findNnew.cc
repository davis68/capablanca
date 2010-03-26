#include <cassert>
#include <iostream>
#include <map>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <cstdio>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "neighbors.h"
#include "particle.h"

#include <fstream>
#include <sstream>
#include <string>

#define  SENT_PARTICLE_SIZE 5

using namespace std;

typedef vector<Particle>    Cell;

extern bool verbose;

extern uint numStates,
            dissolnStates;

extern vector<Particle>    particles;
extern map<id_t, Particle> pmap;

extern int rank,
           size;

extern coord_t myMinX, myMinY, myMinZ,
               myMaxX, myMaxY, myMaxZ,
               minY, maxY;

uint numCellsX, numCellsY, numCellsZ;
vector<vector<vector<Cell> > > cells;

/*  void findMaxNN()
 *  
 *  Find the expected number of maximum nearest neighbors. ***
 */
void findMaxNN()
{   maxNN = MAX_NEIGHBORS; }

/*  void checkBorder(Particle&, const coord_t)
 *  
 *  Find whether the particle lies within the cutoff radius of the border.
 */
inline void checkBorder(Particle& p, const coord_t y)
{   if(abs(y - p.y) <= BORDER_CUTOFF) p.onBoundary = true; }

/*  double distSqrd(Particle& p1, Particle& p2)
 *  
 *  Return square of Euclidean distance between two particles.
 */
inline double distSqrd(Particle& p1, Particle& p2)
{   double x, y, z;
    x = p2.x - p1.x;
    y = p2.y - p1.y;
    z = p2.z - p1.z;
    return (x * x + y * y + z * z); }

/*  uint getNumCells(const coord_t)
 *  
 *  Find the appropriate number of cells in a given dimension. *** optimize
 */
inline uint getNumCells(const coord_t range)
{   /*double  inv3 = 1.0f / 3.0f;
    if (range != floor(range))
        return ceil(inv3 * range);
    else
        return inv3 * floor(range);*/ //***
    return NUM_CELLS; }

/*  void createCells()
 *  
 *  Divide the volume allotted to this processor into cells and assign
 *  each particle to the appropriate cell.  The result is held in cells.
 */
inline void createCells()
{   //  Initialize cells.
    coord_t xRange = myMaxX - myMinX;
    coord_t yRange = myMaxY - myMinY;
    coord_t zRange = myMaxZ - myMinZ;
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
        if (l == 0)                         checkBorder(p, myMinY);
        else if (l == numCellsY - 1)        checkBorder(p, myMaxY);
        
        cells[k][l][m].push_back(p);
        particles.pop_back(); } }

/*  void makeNeighbors(Particle&, Particle&)
 *  
 *  Make two particles have each other's states & ids in their neighbor list.
 */
inline void makeNeighbors(Particle& p1, Particle& p2)
{   int numNN = accumulate(p1.countN.begin(), p1.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particle %d: %d.", p1.id, numNN);
        warning(err); }
    numNN = accumulate(p2.countN.begin(), p2.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particles %d: %d.", p2.id, numNN);
        warning(err); }
    p1.neighbors[p2.state].push_back(p2.id);
    p1.countN[p2.state]++;
    p2.neighbors[p1.state].push_back(p1.id);
    p2.countN[p1.state]++; }

/*  void makeNeighborsAffectOnlyLeft(Particle&, Particle&)
 *  
 *  Make a particle have another's state & id in its neighbor list.
 */
inline void makeNeighborsAffectOnlyLeft(Particle& p1, Particle& p2)
{   int numNN = accumulate(p1.countN.begin(), p1.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particle %d: %d.", p1.id, numNN);
        warning(err); }
    p1.neighbors[p2.state].push_back(p2.id);
    p1.countN[p2.state]++; }

/*  void findIntraCellNeighbors(Cell&, Cell&)
 *  
 *  Test each particle within a cell to see if they are neighbors.
 */
inline void findIntraCellNeighbors(Cell& cell1)
{   for(uint i = 0; i < cell1.size(); i++)
    {   for(uint j = i + 1; j < cell1.size(); j++)
        {   if(distSqrd(cell1[i], cell1[j]) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(cell1[i], cell1[j]); } } } }

/*  void findInterCellNeighbors(Cell&, Cell&)
 *  
 *  Test each particle within two cells to see if they are neighbors.
 */
inline void findInterCellNeighbors(Cell& cell1, Cell& cell2)
{   for(uint i = 0; i < cell1.size(); i++)
    {   for(uint j = 0; j < cell2.size(); j++)
        {   if(distSqrd(cell1[i], cell2[j]) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(cell1[i], cell2[j]); } } } }

/*  void findLocalNeighbors()
 *  
 *  Test all particles on this processor to see if they are neighbors.
 */
inline void findLocalNeighbors()
{   int currProg, maxProg;
    for(uint i = 0; i < numCellsX; i++)
    {   for(uint j = 0; j < numCellsY; j++)
        {   for(uint k = 0; k < numCellsZ; k++)
            {   //  Find local neighbors.
                findIntraCellNeighbors(cells[i][j][k]);
                //  Find adjacent neighbors.
                if(i > 0) findInterCellNeighbors(cells[i][j][k],
                                                 cells[i - 1][j][k]);
                if(j > 0) findInterCellNeighbors(cells[i][j][k],
                                                 cells[i][j - 1][k]);
                if(k > 0) findInterCellNeighbors(cells[i][j][k],
                                                 cells[i][j][k - 1]);
                //  Find diagonal neighbors.
                if(i > 0 && j > 0) findInterCellNeighbors(cells[i][j][k],
                                                          cells[i - 1][j - 1][k]);
                if(i > 0 && k > 0) findInterCellNeighbors(cells[i][j][k],
                                                          cells[i - 1][j][k - 1]);
                if(j > 0 && k > 0) findInterCellNeighbors(cells[i][j][k],
                                                          cells[i][j - 1][k - 1]); } } } }

/*  void getSlice(coord_t*, int&, const uint)
 *  
 *  Pack a slice of a cell into a buffer for messaging.
 */
inline void getSlice(coord_t* buffer, int& bsize, const uint y)
{   bsize = 0;
    for (uint x = 0; x < numCellsX; x++)
    {   for (uint z = 0; z < numCellsZ; z++)
        {   for (uint i = 0; i < cells[x][y][z].size(); i++)
            {   buffer[bsize++] = (coord_t) cells[x][y][z][i].id;
                buffer[bsize++] =           cells[x][y][z][i].x;
                buffer[bsize++] =           cells[x][y][z][i].y;
                buffer[bsize++] =           cells[x][y][z][i].z;
                buffer[bsize++] = (coord_t) cells[x][y][z][i].state; } } } }

/*  void checkNeighbors(Cell&, coord_t*)
 *  
 *  Find out if the sent data describes a particle neighboring any in this cell.
 */
inline void checkNeighbors(Cell& cell1, coord_t* data)
{   Particle other(data[0], data[1], data[2], data[3], data[4]);
    for(uint i = 0; i < cell1.size(); i++)
    {   if(distSqrd(cell1[i], other) <= NEIGHBOR_SQUARE_CUTOFF)
        {   makeNeighborsAffectOnlyLeft(cell1[i], other); } } }

/*  void findRemoteNeighbors()
 *  
 *  Find neighbors on the borders between processors.
 */
inline void findRemoteNeighbors()
{   int     mpi_buffer_size = MAX_ARRAY_SIZE * 2 * sizeof(coord_t);
    byte    mpi_buffer[mpi_buffer_size];
    coord_t buffer[MAX_ARRAY_SIZE];
    int     bsize;
    MPI_Status  status;
    MPI_Request mpr1;
    
    MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);
    //*** clean up following with MPI_PROC_NULL
    if(rank > 0)
    {   getSlice(buffer, bsize, 0);
        MPI_Ibsend(buffer, bsize, MPI_DOUBLE, rank - 1,
                  bsize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    if(rank < size - 1)
    {   getSlice(buffer, bsize, numCellsY - 1);
        MPI_Ibsend(buffer, bsize, MPI_DOUBLE, rank + 1,
                  bsize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    
    if(rank > 0)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank - 1,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint i = 0; i < status.MPI_TAG; i++)
        {   for(uint x = 0; x < numCellsX; x++)
            {   for(uint z = 0; z < numCellsZ; z++)
                {   checkNeighbors(cells[x][0][z],
                                   &(buffer[i * SENT_PARTICLE_SIZE])); } } } }
    if(rank < size - 1)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank + 1,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint i = 0; i < status.MPI_TAG; i++)
        {   for(uint x = 0; x < numCellsX; x++)
            {   for(uint z = 0; z < numCellsZ; z++)
                {   checkNeighbors(cells[x][numCellsY - 1][z],
                                   &(buffer[i * SENT_PARTICLE_SIZE])); } } } }
    
    MPI_Buffer_detach(mpi_buffer, &mpi_buffer_size); }

/*  void fillMap()
 *  
 *  Map each particle into the calculation map so it can be indexed by ID.
 */
inline void fillMap()
{   for(uint x = 0; x < numCellsX; x++)
    {   for(uint y = 0; y < numCellsY; y++)
        {   for(uint z = 0; z < numCellsZ; z++)
            {   for(uint i = 0; i < cells[x][y][z].size(); i++)
                {   Particle& p = cells[x][y][z].back();
                    p.initialN  = p.countN;
                    pmap.insert(pair<id_t, Particle>(p.id, p));
                    cells[x][y][z].pop_back(); } } } } }

/*  void calculateNeighbors()
 *  
 *  Calculate the number of nearest neighbors of each particle.
 */
void calculateNeighbors()
{   createCells();          if (verbose && rank == 1) cout << "  Cell division complete." << endl;
    findLocalNeighbors();   if (verbose && rank == 1) cout << "  Local neighbors found." << endl;
    findRemoteNeighbors();  if (verbose && rank == 1) cout << "  Remote neighbors found." << endl;
    fillMap();
    
    for (uint i = 0; i < particles.size(); i++)
    {   if (accumulate(particles[i].countN.begin(), particles[i].countN.begin() + dissolnStates, 0) == 0) cerr << particles[i].id << "!" << endl; } //***
     }

