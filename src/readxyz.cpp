/** readxyz.cpp
 *  30 Nov 2009--02 Jun 2010
 *  Bill Tuohy and Neal Davis
 */

#include <algorithm>
#include <cassert> //FIXME
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "particle.h"
#include "readxyz.h"

#define  TEMP_PARTICLE_SIZE 4

using namespace std;

//  External global variables from global.cpp
extern char *dataFileName;
extern bool verbose;

extern int  rank,
            size;

extern vector<Particle> particles;

extern coord_t myMinX, myMinY, myMinZ,
               myMaxX, myMaxY, myMaxZ, minY, maxY;

extern uint initialTotalParticles;

//  Global variables
class TempParticle
{   public:
        coord_t x, y, z, state;
        
        TempParticle(coord_t _x, coord_t _y, coord_t _z, state_t _state)
        {   x = _x;
            y = _y;
            z = _z;
            state = (coord_t) _state; } };

int interval,
    idCount;

/*  getProc(coord_t y)
 *  
 *  Assign the node number for a particle based on its y coordinate.  Normalize
 *  to Y range and scale to number of nodes, convert to int.  Check for exactly
 *  equality to the upper limit so the integer isn't bigger than highest rank.
 */
unsigned int getProc(coord_t y)
{   coord_t range   = maxY - minY,
            procFloat = ((y - minY) / range) * size;
    int     procNum   = (int) procFloat;
    if (procNum == size) procNum--;
    
    return procNum; }

/*  int assignID()
 *  
 *  Generate a unique ID for a given particle based on the rank.
 */
int assignID()
{   //TODO:static int idCount = 0;
    return (rank * interval + idCount++); }

/*  void setMinMax()
 *  
 *  Find the minima and maxima in each dimension.
 */
inline void setMinMax (Particle *p1)
{   if (p1->x < myMinX) myMinX = p1->x;
    if (p1->x > myMaxX) myMaxX = p1->x;
    if (p1->y < myMinY) myMinY = p1->y;
    if (p1->y > myMaxY) myMaxY = p1->y;
    if (p1->z < myMinZ) myMinZ = p1->z;
    if (p1->z > myMaxZ) myMaxZ = p1->z; }

/*  void readXYZ()
 *  
 *  Reads XYZ file, builds particle struct for each particle, and separates the
 *  particles into buckets for each processor to own.
 *  
 *  The first line of XYZ file is an integer with the number of particles in the
 *  file.  The rest of the lines contain 4 tab-separated fields,
 *      id, x, y, z.
 *  
 *  This function reads the first line to get the number of particles, then
 *  loops that many times to read lines and split the fields.  Particles are
 *  constructed with the fields from each line and then pushed onto a vector.
 */
void readXYZ()
{   coord_t x, y, z;
    state_t state;
    int     thisProc;    //  Rank to which a particle will be assigned.
    interval = idCount = 0;
    
    vector<TempParticle*> tempParticles;
    double  localMinY = 1e8, localMaxY = -1e8;
    
    // 2-D array to hold raw particle data per destination node
    vector<vector<coord_t> > sendBuckets(size);
    vector<vector<coord_t> > recvBuckets(size);
    
    FILE *xyzFile;
    xyzFile = fopen(dataFileName, "r");
    if (xyzFile == NULL)
    {   char err[64];
        sprintf(err, "⚠ Unable to open file %s", dataFileName);
        error(err);
        return; }
    
    try
    {   int cnt = fscanf(xyzFile, "%u", &initialTotalParticles);
        
        //  Determine file size.
        long int curr_fptr;             //  Test file ptr finding EOL, EOF data.
        curr_fptr = ftell(xyzFile);     //  Note that this is currently EOL 1.
        fseek(xyzFile, 0, SEEK_END);
        int fileBytes   = ftell(xyzFile) - curr_fptr,       //  Total bytes of particle data in file.
            bytesPerNode= fileBytes / size,                 //  Bytes of particle data per proc.
            fileStartPtr= bytesPerNode * rank + curr_fptr,  //  Start of data for this proc.
            fileEndPtr  = fileStartPtr + bytesPerNode - 1;  //  End of data for this proc.
        
        //  Align start/end ptrs to EOLs.
        fseek(xyzFile, fileStartPtr, SEEK_SET);
        int c = fgetc(xyzFile);
        while(c != '\n' && c != EOF)
        {   c = fgetc(xyzFile); }
        fileStartPtr = ftell(xyzFile) - 1;
        
        fseek(xyzFile, fileEndPtr, SEEK_SET);
        c = fgetc(xyzFile); //FIXME:  while((c = fgetc(xyzFile)) != '\n' && c != EOF) { }
        while(c != '\n' && c != EOF)
        {   c = fgetc(xyzFile); }
        fileEndPtr = ftell(xyzFile) - 1;
        if (rank == size - 1)                       //  Force last proc to EOF.
        {   fseek(xyzFile, -1, SEEK_END);
            fileEndPtr = ftell(xyzFile); }
        fseek(xyzFile, fileStartPtr, SEEK_SET);    //  Seek to start ptr.
        
        //  Get interval as estimated number of particles plus one thousand.
        int bytesInFile = fileEndPtr - fileStartPtr,
            numParticlesInFile = bytesInFile / (34 * sizeof(char)); //*** validate
        interval = numParticlesInFile / size + 1000;
        
        while(ftell(xyzFile) != fileEndPtr)
        {   cnt = fscanf(xyzFile, "%u %lf %lf %lf", &state, &x, &y, &z);
            
            //  Assign to localMinY, localMaxY.
            if(y < localMinY)   localMinY = y;
            if(y > localMaxY)   localMaxY = y;
            
            tempParticles.push_back(new TempParticle(x, y, z, state)); }
        fclose(xyzFile); }
    catch (...)
    {   char err[64];
        sprintf(err, "⚠ Unable to read data file %s.", dataFileName);
        error(err); }
    
    //  Find global min and max Y-coordinates.
    MPI_Allreduce(&localMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Read and distribute from vector of id, x, y, z data.
    for(uint i = 0; i < tempParticles.size(); i++)
    {   thisProc = getProc(tempParticles[i]->y);
        sendBuckets[thisProc].push_back(tempParticles[i]->x);
        sendBuckets[thisProc].push_back(tempParticles[i]->y);
        sendBuckets[thisProc].push_back(tempParticles[i]->z);
        sendBuckets[thisProc].push_back(tempParticles[i]->state); }
    
    //  Track particles retained locally, requiring only one call to size().
    uint kept = sendBuckets[rank].size() / TEMP_PARTICLE_SIZE;
    
    //  Send particle data to each node from the buckets created by readXYZ.
    //  These portions were left as partially blocking to not overflow the
    //  receive buffer.
    //    Send and receive number of particles to be transferred.
    MPI_Request mpr1;
    uint        sendSize, recvSize;
    for (int dest = 0; dest < size; dest++)
    {   if (dest != rank)
        {   sendSize = sendBuckets[dest].size();
            MPI_Isend(&sendSize, 1, MPI_INT, dest, 1, MPI_COMM_WORLD, &mpr1); } }
    for (int src = 0; src < size; src++)
    {   if (src != rank)
        {   MPI_Recv(&recvSize, 1, MPI_INT, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recvBuckets[src].resize(recvSize); } }
    
    //    Receive particle data and handle it.
    for (int dest = 0; dest < size; dest++)
    {   if (dest != rank)
        {   MPI_Isend(&(sendBuckets[dest][0]), sendBuckets[dest].size(),
                      MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &mpr1); } }
    for (int src = 0; src < size; src++)
    {   if (src != rank)
        {   MPI_Recv(&(recvBuckets[src][0]), recvBuckets[src].size(),
                     MPI_DOUBLE, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); } }
    
    for (int src = 0; src < size; src++)
    {   if (src != rank)
        {   for (uint i = 0; i < recvBuckets[src].size(); i += TEMP_PARTICLE_SIZE)
            {   Particle p1 ;
                p1.x    =           recvBuckets[src][i + 0];
                p1.y    =           recvBuckets[src][i + 1];
                p1.z    =           recvBuckets[src][i + 2];
                p1.state= (state_t) recvBuckets[src][i + 3];
                particles.push_back(p1);
                setMinMax(&p1); } } }
    
    //  Copy kept data back in from sendBucket.
    for (uint i = 0; i < kept * TEMP_PARTICLE_SIZE; i += TEMP_PARTICLE_SIZE)
    {   Particle p1 ;
        p1.x    =           sendBuckets[rank][i + 0];
        p1.y    =           sendBuckets[rank][i + 1];
        p1.z    =           sendBuckets[rank][i + 2];
        p1.state= (coord_t) sendBuckets[rank][i + 3];
        particles.push_back(p1);
        setMinMax(&p1); }
    
    //  Validate that the number of particles loaded is identical to the
    //  original number of particles in the XYZ file.
    uint localParticleCount = particles.size(),
         totalParticleCount;
    MPI_Allreduce(&localParticleCount, &totalParticleCount, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (verbose && !rank && (initialTotalParticles != totalParticleCount))
        cout << "⚠ " << initialTotalParticles << " particles expected; "
             << totalParticleCount << " found.\n";
    //assert(totalParticleCount == initialTotalParticles); ***
    
    //  Assign particle IDs.
    int *particleCount; particleCount = new int[size];
    MPI_Allgather(&localParticleCount, 1, MPI_INT,
                  particleCount, 1, MPI_INT, MPI_COMM_WORLD);
    int basis = accumulate(particleCount, particleCount + rank, 0);
    for (uint i = 0; i < particles.size(); i++)
    {   particles[i].id = i + basis; } }

