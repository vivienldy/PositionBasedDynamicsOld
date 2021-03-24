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
	SpatialHashSystem(BufferVector3f& posBuffer, BufferInt& triIdxBuffer, HardwareType ht);
	SpatialHashSystem(BufferVector3f& posBuffer, BufferInt& triIdxBuffer, HardwareType ht, float3 gridCenter, uint3 gridSize, float3 cellSize);
	~SpatialHashSystem() {};

	void SetGridCenter(float3 worldOrigin);
	void SetGridSize(uint3 gridSize);
	void SetDivision(float3 cellSize);  // change voxel size / cell size 
	void InitSH();
	void UpdateSH(BufferVector3f& prdPBuffer);  // reserved for dynamic update for accelerating PBD collision detection
	void FindNeighbors(BufferInt& neighbors,  // output: a list of triangle IDs (int)
		int targetTriID);   // input: a triangle ID 
	void SetTimer(Timer* timer) { this->m_shsTimer = timer; }
	// Test Method
	void SetTargetTriId(int targetId) { m_targetId = targetId; }

protected: // functions
	void evalTriCentroids();  // evaluate the centroid for each triangle
	void calcHashCPU();   // calculate the hash for each triangle in the grid basd on the grid cells
	void calcHashGPU(int* gridTriHash,	 // output
		int* gridtriIndex,	 // output
		float* pos,			 // input: positions
		int    numTriangles);   // calculate the hash for each triangle in the grid basd on the grid cells
	void sortTrianglesCPU();
	void sortTrianglesGPU(int* dGridParticleHash,   // will change
		int* dGridParticleIndex,  // will change
		int numParticles);  // input
	void findCellStartCPU();


	// TODO: Reserved for Optimization 
	void reorderDataAndFindCellStart(int* cellStart,
		int* cellEnd,
		float* sortedPos,
		float* sortedVel,
		int* gridParticleHash,
		int* gridParticleIndex,
		float* oldPos,
		float* oldVel,
		int   numParticles,
		int   numCells);

	int3 calcGridPosCPU(float3 triPos);
	int calcGridHashCPU(int3 gridPos);

	bool outsideGrid(int3 gridPos);

	void checkBuffersSizes();

protected: // data
	// Grid configuration
	HardwareType m_ht;
	uint3 m_gridSize;
	uint m_numTriangles;
	int m_numGridCells;
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

	// bool for optimizing when the entire object is outside the grid
	bool m_hAllOutOfGrid = true;

	// Test var
	int m_targetId;

	// Timer 
	Timer* m_shsTimer;
	// bool m_timerStatus;  // 1 - on; 0 - off
};

#endif