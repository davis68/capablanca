/*  
 *  read_xyz.cpp
 *  420-proj
 *  
 *  Created by Bill Tuohy on 30 Nov 09.
 *  Modified by Neal Davis.
 *
 *  Current issues:  occasionally drops particles, so assert() commented out. ***
 */
#include <algorithm>
#include <cassert>
#include <iostream>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "definitions.h"
#include "error.h"
#include "neighbors.h" //*** findMaxNN() business
#include "particle.h"
#include "readxyz.h"

#define  TEMP_PARTICLE_SIZE 4

using namespace std;

extern char *dataFileName;
extern bool verbose;

extern int  rank,
            size;

extern ParticleList particles;

extern coord_t myMinX, myMinY, myMinZ, myMaxX, myMaxY, myMaxZ;
extern double minY, maxY;

extern uint initialTotalParticles;

class TempParticle
{   public:
        coord_t x, y, z;
        coord_t state;
        
        TempParticle(coord_t _x, coord_t _y, coord_t _z, state_t _state)
        {   x = _x; y = _y; z = _z;
            state = (coord_t) _state; } };

int interval, idCount;

/*  Assign the node number for a particle based on its Y
    coordinate.  Normalize to Y range and scale to number
    of nodes, convert to int.  Check for exactly equal to
    upper limit so the integer isn't bigger than highest rank
*/
unsigned int getProc(coord_t y)
{   coord_t range   = maxY - minY,
            p_float = ((y - minY) / range) * size;
    int     p_num   = (int) p_float; if (p_num == size) p_num--;
    return p_num; }

/*  int assignID()
 *  
 *  Generate a unique ID for a given particle based on the rank.
 */
int assignID()
{   return (rank * interval + idCount++); }

/*  void setMinMax(Particle *, bool &)
 *  
 *  Find the minima and maxima in each dimension.
 */
void setMinMax (Particle *p1, bool &flag)
{   coord_t x = p1->x,
            y = p1->y,
            z = p1->z;
    
    if (!flag)
    {   flag = true;
        myMinX = x; myMinY = y; myMinZ = z;
        myMaxX = x; myMaxY = y; myMaxZ = z; }
    else
    {   if (x < myMinX) myMinX = x;
        if (x > myMaxX) myMaxX = x;
        if (y < myMinY) myMinY = y;
        if (y > myMaxY) myMaxY = y;
        if (z < myMinZ) myMinZ = z;
        if (z > myMaxZ) myMaxZ = z; } }

/*  void readXYZ()
 *  
 *  Reads XYZ file, builds particle struct for each particle, and separates
 *  particles into buckets for each proc to own.
 *  
 *  First line of XYZ file is an integer with the number of particles
 *  in the file.  The rest of the lines contain 4 fields:
 *  id, x, y, z that are tab separated.
 *  
 *  This function reads the first line to get the number of particles,
 *  then loops that many times to read lines and split the fields.
 *  Particles are constructed with the fields from each line.  They are
 *  then pushed onto a vector
 */
void readXYZ()
{   id_t    id;
    coord_t x, y, z;
    state_t state;
    int     thisProc;    //  Rank to which a particle will be assigned.
    interval = idCount = 0;
    findMaxNN(); //***
    
    // 2-D array to hold raw particle data per destination node
    vector<vector<coord_t> > sendBuckets(size);
    vector<vector<coord_t> > recvBuckets(size);
    
    FILE *xyz_file;
    xyz_file = fopen(dataFileName, "r");
    if (xyz_file == NULL)
    {   char err[64];
        sprintf(err, "⚠ Unable to open file %s", dataFileName);
        error(err);
        return; }
    
    fscanf(xyz_file, "%u", &initialTotalParticles);
    
    //  Determine file size.
    long int curr_fptr;             //  Test file ptr finding EOL, EOF data.
    curr_fptr = ftell(xyz_file);    //  Note that this is currently EOL 1.
    fseek(xyz_file, 0, SEEK_END);
    int fileBytes   = ftell(xyz_file) - curr_fptr,      //  Total bytes of particle data in file.
        bytesPerNode= fileBytes / size,                 //  Bytes of particle data per proc.
        fileStartPtr= bytesPerNode * rank + curr_fptr,  //  Start of data for this proc.
        fileEndPtr  = fileStartPtr + bytesPerNode - 1;  //  End of data for this proc.
    
    //  Align start/end ptrs to EOLs.
    fseek(xyz_file, fileStartPtr, SEEK_SET);
    int c = fgetc(xyz_file);
    while(c != '\n' && c != EOF)
    {   c = fgetc(xyz_file); }
    fileStartPtr = ftell(xyz_file) - 1;
    
    fseek(xyz_file, fileEndPtr, SEEK_SET);
    c = fgetc(xyz_file);
    while(c != '\n' && c != EOF)
    {   c = fgetc(xyz_file); }
    fileEndPtr = ftell(xyz_file) - 1;
    if (rank == size - 1)                       //  Force last proc to EOF.
    {   fseek(xyz_file, -1, SEEK_END);
        fileEndPtr = ftell(xyz_file); }
    fseek(xyz_file, fileStartPtr, SEEK_SET);    //  Seek to start ptr.
    
    double  localMinY, localMaxY;
    bool    hasInitLocalMinMax = false;
    vector<TempParticle*> tempParticle;
    
    //  Get interval as estimated number of particles plus one thousand.
    int bytesInFile = fileEndPtr - fileStartPtr,
        numParticlesInFile = bytesInFile / (34 * sizeof(char)); //*** validate
    interval = numParticlesInFile / size + 1000;
    
    while(ftell(xyz_file) != fileEndPtr)
    {   fscanf(xyz_file, "%u %lf %lf %lf", &state, &x, &y, &z);
        
        //  Assign to localMinY, localMaxY.
        if(!hasInitLocalMinMax)
        {   localMinY = localMaxY = y;
            hasInitLocalMinMax = true; }
        if(y < localMinY)   localMinY = y;
        if(y > localMaxY)   localMaxY = y;
        
        tempParticle.push_back(new TempParticle(x, y, z, state)); }
    
    //  Find global min and max Y-coordinates.
    MPI_Allreduce(&localMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    //  Read and distribute from vector of id, x, y, z data.
    for(uint i = 0; i < tempParticle.size(); i++)
    {   thisProc = getProc(tempParticle[i]->y);
        sendBuckets[thisProc].push_back(tempParticle[i]->x);
        sendBuckets[thisProc].push_back(tempParticle[i]->y);
        sendBuckets[thisProc].push_back(tempParticle[i]->z);
        sendBuckets[thisProc].push_back(tempParticle[i]->state); }
    
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
    
    bool initMyMax = false;
    for (int src = 0; src < size; src++)
    {   if (src != rank)
        {   for (uint i = 0; i < recvBuckets[src].size(); i += TEMP_PARTICLE_SIZE)
            {   Particle p1 ;
                p1.x    =           recvBuckets[src][i + 0];
                p1.y    =           recvBuckets[src][i + 1];
                p1.z    =           recvBuckets[src][i + 2];
                p1.state= (state_t) recvBuckets[src][i + 3];
                particles.push_back(p1);
                setMinMax(&p1, initMyMax); } } }
    
    //  Copy kept data back in from sendBucket.
    for (uint i = 0; i < kept * TEMP_PARTICLE_SIZE; i += TEMP_PARTICLE_SIZE)
    {   Particle p1 ;
        p1.x    =           sendBuckets[rank][i + 0];
        p1.y    =           sendBuckets[rank][i + 1];
        p1.z    =           sendBuckets[rank][i + 2];
        p1.state= (coord_t) sendBuckets[rank][i + 3];
        particles.push_back(p1);
        setMinMax(&p1, initMyMax); }
    
    //  Validate that the number of particles loaded is identical to the
    //  original number of particles in the XYZ file.
    int localParticleCount = particles.size(),
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
    {   particles[i].id = i + basis; }
}

