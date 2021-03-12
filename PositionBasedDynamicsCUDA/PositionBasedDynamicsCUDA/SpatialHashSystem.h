#ifndef _SPATIAL_HASH_SYS_H_
#define _SPATIAL_HASH_SYS_H_
#include "cuda_runtime.h"
#include"PBD_Basic.cuh"
#include <assert.h>
#include "helper_math.h"

class SpatialHashSystem
{
public:
	SpatialHashSystem(string filename, uint3 gridSize, HardwareType ht);
	SpatialHashSystem(BufferVector3f posBuffer, BufferInt triIdxBuffer, HardwareType ht);
	//SpatialHashSystem(BufferVector3f posBuffer, BufferInt triIdxBuffer, uint numTriangles, uint3 gridSize, float oneCellSize, HardwareType ht);
	~SpatialHashSystem() {};

	void SetGridCenter(float3 worldOrigin);
	void SetGridSize(uint3 gridSize);
	void SetDivision(float3 cellSize);  // change voxel size / cell size 
	void InitSH();
	void UpdateSH(float dt);  // reserved for dynamic update for accelerating PBD collision detection
	void FindNeighbors(BufferInt& neighbors,  // output: a list of triangle IDs (int)
		uint targetTriID);   // input: a triangle ID 
// Test Method
	void SetTargetTriId(uint targetId) { m_targetId = targetId; }

protected: // functions
	// I/O
	void readMeshFromTxt(string filename);
	void splitString(const std::string& src, std::vector<std::string>& v, const std::string& split);
	std::vector<std::string> splitString(const std::string& src, const std::string& split);

	// void initSpatialHash(uint1 numTriangles, uint3 gridSize, HardwareType ht);
	// void initSpatialHash(uint1 numTriangles, uint3 gridSize, float cellSize, HardwareType ht);
	// void initTriangleIndex();  // initialize m_hTriangleIndex
	void evalTriCentroids();  // evaluate the centroid for each triangle
	void calcHashCPU();   // calculate the hash for each triangle in the grid basd on the grid cells
	void calcHashGPU(uint* gridTriHash,	 // output
		uint* gridtriIndex,	 // output
		float* pos,			 // input: positions
		int    numTriangles);   // calculate the hash for each triangle in the grid basd on the grid cells
	void sortTrianglesCPU();
	void sortTrianglesGPU(uint* dGridParticleHash,   // will change
		uint* dGridParticleIndex,  // will change
		uint numParticles);  // input
	void findCellStartCPU();


	// TODO: Reserved for Optimization 
	void reorderDataAndFindCellStart(uint* cellStart,
		uint* cellEnd,
		float* sortedPos,
		float* sortedVel,
		uint* gridParticleHash,
		uint* gridParticleIndex,
		float* oldPos,
		float* oldVel,
		uint   numParticles,
		uint   numCells);

	int3 calcGridPosCPU(float3 triPos);
	uint calcGridHashCPU(int3 gridPos);

protected: // data
	// Grid configuration
	HardwareType m_ht;
	uint3 m_gridSize;
	uint m_numTriangles;
	uint m_numGridCells;
	bool m_initialized = false;
	float3 m_worldOrigin;
	float3 m_girdStartPos;
	float3 m_cellSize;  // voxel size

	// I/O
	bool m_loadSuccess = false;
	std::vector<std::string> dataStringList;  // intermediate vector storing strings from file

	// CPU data
	BufferVector3f m_hPosBuffer;   // triangles positions (read from txt for now)
	BufferInt m_hTriIdxBuffer;   // triangle indices (read from txt for now)
	BufferVector3f m_hTriCentroidBuffer;  // triangle centriod (used for determining the grid pos)
	// grid data for sorting method
	BufferInt m_hTriangleHash;  // grid hash value for each triangle (aka cell id)
	BufferInt m_hTriangleId; // particle index for each triangle (aka triangle id)
	BufferInt m_hCellStart;     // index of start of each cell in sorted list
	BufferInt m_hCellEnd;		// index of end of cell

	// Test var
	uint m_targetId;

};

#endif

