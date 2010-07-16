/** neighbors.cpp
 *  29 Oct 2009--14 Jul 2010
 *  Neal Davis, Minas Charalambides, and Simon Jenkins
 *  
 *  This file uses MPI buffered send mode, which precludes any other portion
 *  of the program from using this mode simultaneously.
 */

#include <algorithm>
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

#define  SENT_PARTICLE_SIZE 5

using namespace std;

typedef vector<Particle> Cell;

//  External global variables from global.cpp
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

//  Global variables
uint numCellsX, numCellsY, numCellsZ;
vector<vector<vector<Cell> > > cells;

/*  void findMaxNN()
 *  
 *  Find the expected number of maximum nearest neighbors.
 */
void findMaxNN()
{   if (!rank)
    {   Particle i1, i2;
        i1 = particles[0];
        for (vector<Particle>::iterator iter = particles.begin() + 1; iter != particles.end(); iter++)
        {   if (distSqrd(i1, *iter) <= NEIGHBOR_SQUARE_CUTOFF)
            {   i2 = *iter; break; } }
        
        //  Determine orientation type by the distance in each dimension.
        double  dx = abs(i1.x - i2.x),
                dy = abs(i1.y - i2.y),
                dz = abs(i1.z - i2.z),
                d  = max(dx, max(dy, dz));
        if (abs(d - 0.707) < 1e-2)      maxNN = 12;     //  FCC
        else if (abs(d - 0.577) < 1e-2) maxNN = 8;      //  BCC
        else                            maxNN = 6; }    //  SC
    MPI_Bcast(&maxNN, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD); }

/*  void checkBorder(Particle&, const coord_t)
 *  
 *  Find whether the particle lies within a cutoff radius of the border.
 */
inline void checkBorder(Particle& p, const coord_t y)
{   if(abs(y - p.y) <= NEIGHBOR_SQUARE_CUTOFF) p.onBoundary = true; }

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
 *  Find the appropriate number of cells in a given dimension.
 */
inline uint getNumCells(const coord_t range)
{   static double inv3 = 1.0 / 3.0;
    static uint numCells;
    
    numCells = floor(range * inv3);
    if (numCells <= 0) numCells = 1;
    
    return numCells; }

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
        checkBorder(p, myMinY);
        checkBorder(p, myMaxY);
        
        cells[k][l][m].push_back(p);
        particles.pop_back(); } }

/*  bool alreadyNeighbors()
 *  
 *  Indicate whether two particles are already neighbors to each other.
 */
inline bool alreadyNeighbors(Particle& p1, Particle& p2)
{   for (vector<id_t>::iterator iter = p1.neighbors.begin(); iter != p1.neighbors.end(); iter++)
    {   if (find(p2.neighbors.begin(), p2.neighbors.end(), *iter) != p2.neighbors.end()) return true; }
    
    return false; }

/*  void makeNeighbors(Particle&, Particle&)
 *  
 *  Make two particles have each other's states & ids in their neighbor list.
 */
inline void makeNeighbors(Particle& p1, Particle& p2)
{   if (alreadyNeighbors(p1, p2))   return;
    
    static uint numNN;
    numNN = accumulate(p1.countN.begin(), p1.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particle %d: %d.", p1.id, numNN);
        warning(err); }
    numNN = accumulate(p2.countN.begin(), p2.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particles %d: %d.", p2.id, numNN);
        warning(err); }
    p1.neighbors.push_back(p2.id);
    p1.countN[p2.state]++;
    p2.neighbors.push_back(p1.id);
    p2.countN[p1.state]++; }

/*  void makeNeighborsAffectOnlyLeft(Particle&, Particle&)
 *  
 *  Make a particle have another's state & id in its neighbor list.
 */
inline void makeNeighborsAffectOnlyLeft(Particle& p1, Particle& p2)
{   if (alreadyNeighbors(p1, p2))   return;
    
    static uint numNN;
    numNN = accumulate(p1.countN.begin(), p1.countN.begin() + dissolnStates, 0);
    if(numNN > maxNN)
    {   char err[64];
        sprintf(err, "Too many neighbors detected for particle %d: %d.", p1.id, numNN);
        warning(err); }
    p1.neighbors.push_back(p2.id);
    p1.countN[p2.state]++; }

/*  void findIntraCellNeighbors(Cell&, Cell&)
 *  
 *  Test each particle within a cell to see if they are neighbors.
 */
inline void findIntraCellNeighbors(Cell& cell1)
{   for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   for(Cell::iterator iter2 = iter1 + 1; iter2 != cell1.end(); iter2++)
        {   if(distSqrd(*iter1, *iter2) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(*iter1, *iter2); } } } }

/*  void findInterCellNeighbors(Cell&, Cell&)
 *  
 *  Test each particle within two cells to see if they are neighbors.
 */
inline void findInterCellNeighbors(Cell& cell1, Cell& cell2)
{   for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   for(Cell::iterator iter2 = cell2.begin(); iter2 != cell2.end(); iter2++)
        {   if(distSqrd(*iter1, *iter2) <= NEIGHBOR_SQUARE_CUTOFF)
            {   makeNeighbors(*iter1, *iter2); } } } }

/*  void findLocalNeighbors()
 *  
 *  Test all particles on this processor to see if they are neighbors.
 */
inline void findLocalNeighbors()
{   for(uint i = 0; i < numCellsX; i++)
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
inline void getSlice(coord_t* buffer, int& bufferSize, const uint y)
{   bufferSize = 0;
    for (uint x = 0; x < numCellsX; x++)
    {   for (uint z = 0; z < numCellsZ; z++)
        {   for (Cell::iterator iter = cells[x][y][z].begin(); iter != cells[x][y][z].end(); iter++)
            {   if ((*iter).onBoundary)
                {   buffer[bufferSize++] = (coord_t) (*iter).id;
                    buffer[bufferSize++] =           (*iter).x;
                    buffer[bufferSize++] =           (*iter).y;
                    buffer[bufferSize++] =           (*iter).z;
                    buffer[bufferSize++] = (coord_t) (*iter).state; } } } } }


/*  void getBSlice(coord_t*, int&, const uint) FIXME
 *  
 *  Pack a slice of a cell into a buffer for messaging.
 */
inline void getBSlice(coord_t* buffer, int& bufferSize, const uint y)
{   bufferSize = 0;
    for (uint x = 0; x < numCellsX; x++)
    {   for (uint z = 0; z < numCellsZ; z++)
        {   for (uint i = 0; i < cells[x][y][z].size(); i++)
            {   buffer[bufferSize++] = (coord_t) cells[x][y][z][i].id;
                buffer[bufferSize++] =           cells[x][y][z][i].x;
                buffer[bufferSize++] =           cells[x][y][z][i].y;
                buffer[bufferSize++] =           cells[x][y][z][i].z;
                buffer[bufferSize++] = (coord_t) cells[x][y][z][i].state; } } } }

/*  void checkNeighbors(Cell&, coord_t*)
 *  
 *  Find out if the sent data describes a particle neighboring any in this cell.
 */
inline void checkNeighbors(Cell& cell1, coord_t* data)
{   Particle other(data[0], data[1], data[2], data[3], data[4]);
    for(Cell::iterator iter1 = cell1.begin(); iter1 != cell1.end(); iter1++)
    {   if(distSqrd(*iter1, other) <= NEIGHBOR_SQUARE_CUTOFF)
        {   makeNeighborsAffectOnlyLeft(*iter1, other); } } }

/*  void findRemoteNeighbors()
 *  
 *  Find neighbors on the borders between processors.
 */
inline void findRemoteNeighbors()
{   static coord_t  sendBuffer[MAX_ARRAY_SIZE],
                    recvBuffer[MAX_ARRAY_SIZE];
    int sendBufferSize = 0,
        recvBufferSize = 0;
    MPI_Status  status;
    
    //  Determine neighboring processors.
    int rankLeft, rankRight;
    rankRight = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
    rankLeft  = (rank == 0)        ? MPI_PROC_NULL : rank - 1;
    
    //  Even processors send to left, then odd processors send to left.
    getSlice(sendBuffer, sendBufferSize, 0);
    if (rank % 2 == 0)
    {   MPI_Ssend(sendBuffer, sendBufferSize, MPI_DOUBLE, rankLeft,
                  sendBufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD); }
    else
    {   MPI_Recv(recvBuffer + recvBufferSize, MAX_ARRAY_SIZE, MPI_DOUBLE, rankRight,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        recvBufferSize += status.MPI_TAG; }
    if (rank % 2 != 0)
    {   MPI_Ssend(sendBuffer, sendBufferSize, MPI_DOUBLE, rankLeft,
                  sendBufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD); }
    else
    {   MPI_Recv(recvBuffer + recvBufferSize, MAX_ARRAY_SIZE, MPI_DOUBLE, rankRight,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        recvBufferSize += status.MPI_TAG; }
    
    //  Even processors send to right, then odd processors send to right.
    getSlice(sendBuffer, sendBufferSize, numCellsY - 1);
    if (rank % 2 == 0)
    {   MPI_Ssend(sendBuffer, sendBufferSize, MPI_DOUBLE, rankRight,
                  sendBufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD); }
    else
    {   MPI_Recv(recvBuffer + recvBufferSize, MAX_ARRAY_SIZE, MPI_DOUBLE, rankLeft,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        recvBufferSize += status.MPI_TAG; }
    if (rank % 2 != 0)
    {   MPI_Ssend(sendBuffer, sendBufferSize, MPI_DOUBLE, rankRight,
                  sendBufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD); }
    else
    {   MPI_Recv(recvBuffer + recvBufferSize, MAX_ARRAY_SIZE, MPI_DOUBLE, rankLeft,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        recvBufferSize += status.MPI_TAG; }
    
    if (rank > 0)
    {   for(uint i = 0; i < (uint) status.MPI_TAG; i++)
        {   for(uint x = 0; x < numCellsX; x++)
            {   for(uint z = 0; z < numCellsZ; z++)
                {   checkNeighbors(cells[x][0][z],
                                   &(sendBuffer[i * SENT_PARTICLE_SIZE])); } } } }
    if(rank < size - 1)
    {   for(uint i = 0; i < (uint) status.MPI_TAG; i++)
        {   for(uint x = 0; x < numCellsX; x++)
            {   for(uint z = 0; z < numCellsZ; z++)
                {   checkNeighbors(cells[x][numCellsY - 1][z],
                                   &(sendBuffer[i * SENT_PARTICLE_SIZE])); } } } } }

/*  void findBRemoteNeighbors() FIXME
 *  
 *  Find neighbors on the borders between processors using MPI buffered send mode.
 */
inline void findBRemoteNeighbors()
{   int     mpi_buffer_size = MAX_ARRAY_SIZE * 2 * sizeof(coord_t);
    byte    mpi_buffer[mpi_buffer_size];
    coord_t buffer[MAX_ARRAY_SIZE];
    int     bufferSize;
    MPI_Status  status;
    MPI_Request mpr1;
    
    MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);
    
    if(rank > 0)
    {   getBSlice(buffer, bufferSize, 0);
        MPI_Ibsend(buffer, bufferSize, MPI_DOUBLE, rank - 1,
                  bufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    if(rank < size - 1)
    {   getBSlice(buffer, bufferSize, numCellsY - 1);
        MPI_Ibsend(buffer, bufferSize, MPI_DOUBLE, rank + 1,
                  bufferSize / SENT_PARTICLE_SIZE, MPI_COMM_WORLD, &mpr1); }
    
    if(rank > 0)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank - 1,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint i = 0; i < (uint) status.MPI_TAG; i++)
        {   for(uint x = 0; x < numCellsX; x++)
            {   for(uint z = 0; z < numCellsZ; z++)
                {   checkNeighbors(cells[x][0][z],
                                   &(buffer[i * SENT_PARTICLE_SIZE])); } } } }
    if(rank < size - 1)
    {   MPI_Recv(buffer, MAX_ARRAY_SIZE, MPI_DOUBLE, rank + 1,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for(uint i = 0; i < (uint) status.MPI_TAG; i++)
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
            {   for(Cell::iterator iter = cells[x][y][z].begin(); iter != cells[x][y][z].end(); iter++)
                {   Particle& p = *iter;
                    p.onSurface = notOnSurface;
                    pmap.insert(pair<id_t, Particle>(p.id, p)); }
                cells[x][y][z].clear(); } } } }

/*  void calculateNeighbors()
 *  
 *  Calculate the number of nearest neighbors of each particle.
 */
void calculateNeighbors()
{   findMaxNN();
    createCells();          if (verbose && !rank) cout << "  Cell partitioning complete." << endl;  cout.flush();
    findLocalNeighbors();   if (verbose && !rank) cout << "  Local neighbors found." << endl;       cout.flush();
    findBRemoteNeighbors();  if (verbose && !rank) cout << "  Remote neighbors found." << endl;      cout.flush();
    fillMap(); }

