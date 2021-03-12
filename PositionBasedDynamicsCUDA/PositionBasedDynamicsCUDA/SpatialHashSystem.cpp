#include "SpatialHashSystem.h"

SpatialHashSystem::SpatialHashSystem(string filename, uint3 gridSize, HardwareType ht) :
	m_gridSize(gridSize),
	m_ht(ht)
{
	// init m_numTriangles, m_loadSuccess, m_hPosBuffer, m_hTriIdxBuffer
	readMeshFromTxt(filename);
	assert(m_loadSuccess);

	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f); // should be initialized by the user
	m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;
	//cout << "m_numGridCells: " << m_numGridCells << endl;

	m_cellSize = make_float3(2.0f, 2.0f, 2.0f);

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hCellStart.m_Data.resize(m_numTriangles, -1);
	m_hCellEnd.m_Data.resize(m_numTriangles, -1);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));
	m_targetId = 645;
	evalTriCentroids();

}

SpatialHashSystem::SpatialHashSystem(BufferVector3f& posBuffer, BufferInt& triIdxBuffer, HardwareType ht) :
	m_ht(ht)
{
	// init m_numTriangles, m_loadSuccess, m_hPosBuffer, m_hTriIdxBuffer
	m_hPosBuffer = posBuffer;
	m_hTriIdxBuffer = triIdxBuffer;
	m_numTriangles = m_hTriIdxBuffer.GetSize() / 3;
	// assert(m_loadSuccess);

	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f); // should be initialized by the user
	m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;
	//cout << "m_numGridCells: " << m_numGridCells << endl;

	m_cellSize = make_float3(2.0f, 2.0f, 2.0f);

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hCellStart.m_Data.resize(m_numTriangles, -1);
	m_hCellEnd.m_Data.resize(m_numTriangles, -1);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));
	m_targetId = 0;
	evalTriCentroids();

}

void SpatialHashSystem::SetGridCenter(float3 worldOrigin)
{
	m_worldOrigin = worldOrigin;
}

void SpatialHashSystem::SetGridSize(uint3 gridSize)
{
	m_gridSize = gridSize;
	// TODO: update the SH grid
}

void SpatialHashSystem::SetDivision(float3 cellSize)
{
	m_cellSize = cellSize;
	// TODO: update the SH grid
}

void SpatialHashSystem::InitSH()
{
	float shiftX = (m_gridSize.x * m_cellSize.x) / 2;
	float shiftY = (m_gridSize.y * m_cellSize.y) / 2;
	float shiftZ = (m_gridSize.z * m_cellSize.z) / 2;
	m_girdStartPos = m_worldOrigin - make_float3(shiftX, shiftY, shiftZ);
	//printf("grid start pos is %f, %f, %f\n", m_girdStartPos.x, m_girdStartPos.y, m_girdStartPos.z);
	m_initialized = true;
}

/*
// this constructor will be used ONLY if the mesh comes from a file
SpatialHashSystem::SpatialHashSystem(string filename, uint numTriangles, uint3 gridSize, float3 gridStart, float oneCellSize, HardwareType ht) :
	m_numTriangles(numTriangles),
	m_gridSize(gridSize),
	m_ht(ht)
{
	m_cellSize.x = m_cellSize.y = m_cellSize.z = oneCellSize;
	readMeshFromTxt(filename);
	assert(m_loadSuccess);

	// m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f);
	m_worldOrigin = gridStart;
	m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;
	cout << "m_numGridCells: " << m_numGridCells << endl;

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hCellStart.m_Data.resize(m_numTriangles, -1);
	m_hCellEnd.m_Data.resize(m_numTriangles, -1);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));

	evalTriCentroids();

	m_initialized = true;
}

SpatialHashSystem::SpatialHashSystem(BufferVector3f posBuffer, BufferInt triIdxBuffer, uint numTriangles, uint3 gridSize, HardwareType ht) :
	m_numTriangles(numTriangles),
	m_gridSize(gridSize),
	m_ht(ht)
{
	m_hPosBuffer = posBuffer;
	m_hTriIdxBuffer = triIdxBuffer;
	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f);
	m_numGridCells = m_gridSize.x*m_gridSize.y*m_gridSize.z;

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hCellStart.m_Data.resize(m_numTriangles, -1);
	m_hCellEnd.m_Data.resize(m_numTriangles, -1);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));

	evalTriCentroids();

	m_initialized = true;
}

SpatialHashSystem::SpatialHashSystem(BufferVector3f posBuffer, BufferInt triIdxBuffer, uint numTriangles, uint3 gridSize, float oneCellSize, HardwareType ht) :
	m_numTriangles(numTriangles),
	m_gridSize(gridSize),
	m_ht(ht)
{
	m_cellSize.x = m_cellSize.y = m_cellSize.z = oneCellSize;
	m_hPosBuffer = posBuffer;
	m_hTriIdxBuffer = triIdxBuffer;
	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f);
	m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hCellStart.m_Data.resize(m_numTriangles, -1);
	m_hCellEnd.m_Data.resize(m_numTriangles, -1);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));

	evalTriCentroids();

	m_initialized = true;
}
*/

void SpatialHashSystem::UpdateSH(float dt)
{
	assert(m_initialized);
	switch (m_ht)
	{
	case CPU:
		calcHashCPU();
		sortTrianglesCPU();
		findCellStartCPU();
		break;
	case GPU:
		break;
	default:
		break;
	}
}

bool SpatialHashSystem::outsideGrid(int3 gridPos)
{
	if (gridPos.x < 0 || gridPos.y < 0 || gridPos.z < 0 ||
		gridPos.x >= m_gridSize.x || gridPos.y >= m_gridSize.y || gridPos.z >= m_gridSize.z)
		return true;
	else 
		return false;
}

void SpatialHashSystem::FindNeighbors(BufferInt& neighbors,  // output: a list of triangle IDs (int)
	uint targetTriID)   // input: a triangle ID 
{
	neighbors.m_Data.clear();
	float3 targetTriCentroid = m_hTriCentroidBuffer.m_Data[targetTriID];
	// get address in grid
	int3 gridPos = calcGridPosCPU(targetTriCentroid);
	// printf("\t gridPos: %d %d %d\n", gridPos.x, gridPos.y, gridPos.z);
	uint cellStart;
	uint cellEnd;
	for (int z = -1; z <= 1; z++)
	{
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				int3 neighborPos = gridPos + make_int3(x, y, z);
				// check neighbor pos valid
				if (outsideGrid(neighborPos))
					continue;
				//printf("\t valid neighbor pos : (%d, %d, %d)", neighborPos.x, neighborPos.y, neighborPos.z);
				uint hashId = calcGridHashCPU(neighborPos);
				cellStart = m_hCellStart.m_Data[hashId];
				if (cellStart != -1)
				{
					//printf("\t hashId: %d\n", hashId);
					cellEnd = m_hCellEnd.m_Data[hashId];
					// cell range of a hash: [cellStart, cellEnd]
					for (uint i = cellStart; i <= cellEnd; ++i)
					{
						uint neighborId = m_hTriangleId.m_Data[i];
						if (neighborId != targetTriID)
							neighbors.m_Data.push_back(neighborId);
					}
				}
			}
		}
	}
}

// I/O
void SpatialHashSystem::readMeshFromTxt(string filename)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	m_hPosBuffer.m_Data.clear();
	m_hTriIdxBuffer.m_Data.clear();
	m_hTriangleId.m_Data.clear();
	m_loadSuccess = false;

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = splitString(buffer, firstSplit)[1];
		dataStringList = splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				// std::cout << dataStringList[i] << std::endl;
				float3 vertexPos = make_float3(std::stof(dataStringList[i]),
					std::stof(dataStringList[i + 1]),
					std::stof(dataStringList[i + 2]));
				// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				m_hPosBuffer.m_Data.push_back(vertexPos);
			}
			//std::cout << "m_hPosBuffer: " << m_hPosBuffer.GetSize() << endl;
			assert(m_hPosBuffer.GetSize() == dataStringList.size() / 3);
			/*printf("0 (%f, %f, %f),615 (%f, %f, %f),1569 (%f, %f, %f)\n", m_hPosBuffer.m_Data[0].x, m_hPosBuffer.m_Data[0].y, m_hPosBuffer.m_Data[0].z,
				m_hPosBuffer.m_Data[615].x, m_hPosBuffer.m_Data[615].y, m_hPosBuffer.m_Data[615].z, m_hPosBuffer.m_Data[1569].x, m_hPosBuffer.m_Data[1569].y,
				m_hPosBuffer.m_Data[1569].z);*/
		}
		else  // second line of the input file: vertices tet
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				m_hTriIdxBuffer.m_Data.push_back(std::stoi(dataStringList[i]));
				m_hTriIdxBuffer.m_Data.push_back(std::stoi(dataStringList[i + 1]));
				m_hTriIdxBuffer.m_Data.push_back(std::stoi(dataStringList[i + 2]));
				m_hTriangleId.m_Data.push_back((i / 3));
			}
			//std::cout << "m_hTriangleId: " << m_hTriangleId.GetSize() << endl;
			m_numTriangles = m_hTriangleId.GetSize();
			assert(m_hTriangleId.GetSize() == dataStringList.size() / 3);
		}
		++lineNum;
	}
	in.close();
	m_loadSuccess = true;
	//printf("Read File Done!\n");
}

void SpatialHashSystem::splitString(const std::string& src, std::vector<std::string>& v, const std::string& split)
{
	std::string::size_type pos1, pos2;
	pos2 = src.find(split);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(src.substr(pos1, pos2 - pos1));
		pos1 = pos2 + split.size();
		pos2 = src.find(split, pos1);
	}
	if (pos1 != src.length())
		v.push_back(src.substr(pos1));
}

std::vector<std::string> SpatialHashSystem::splitString(const std::string& src, const std::string& split)
{
	std::vector<std::string> _ret = std::vector<std::string>();
	splitString(src, _ret, split);
	return _ret;
}

// evaluate the centroid for each triangle
void SpatialHashSystem::evalTriCentroids()
{
	auto triCentroidBuffer = &(m_hTriCentroidBuffer);
	auto posBuffer = &(m_hPosBuffer);
	auto triIdxBuffer = &(m_hTriIdxBuffer);
	float3 currCentroid;
	for (uint triId = 0; triId < m_numTriangles; ++triId)
	{
		int triStart = triId * 3;
		float3 p1 = posBuffer->m_Data[triIdxBuffer->m_Data[triStart]];
		//if (triId == m_targetId) printf("p1 %f, %f, %f\n", p1.x, p1.y, p1.z);
		float3 p2 = posBuffer->m_Data[triIdxBuffer->m_Data[triStart + 1]];
		//if (triId == m_targetId) printf("p2 %f, %f, %f\n", p2.x, p2.y, p2.z);
		float3 p3 = posBuffer->m_Data[triIdxBuffer->m_Data[triStart + 2]];
		//if (triId == m_targetId)  printf("p3 %f, %f, %f\n", p3.x, p3.y, p3.z);
		currCentroid.x = (p1.x + p2.x + p3.x) / 3;
		currCentroid.y = (p1.y + p2.y + p3.y) / 3;
		currCentroid.z = (p1.z + p2.z + p3.z) / 3;

		triCentroidBuffer->m_Data[triId] = currCentroid;
		//printf("triCentroidBuffer %f, %f, %f\n", m_hTriCentroidBuffer->m_Data[triId].x, m_hTriCentroidBuffer->m_Data[triId].y, m_hTriCentroidBuffer->m_Data[triId].z);
	}
}

void SpatialHashSystem::calcHashCPU()
{
	auto triCentroidBuffer = &(m_hTriCentroidBuffer);
	auto triangleHash = &(m_hTriangleHash);
	auto triangleId = &(m_hTriangleId);
	for (uint triId = 0; triId < m_numTriangles; ++triId)
	{
		float3 triCentroid = triCentroidBuffer->m_Data[triId];
		// get pos in grid
		int3 gridPos = calcGridPosCPU(triCentroid);
		uint hashId = calcGridHashCPU(gridPos);
		// store grid hash and particle index
		triangleHash->m_Data[triId] = hashId;
		triangleId->m_Data[triId] = triId;  // has been initialized if read mesh from txt
		/*if (triId == m_targetId)
		{
			printf("triCentroidBuffer[triId]: (%f, %f, %f)\n", triCentroid.x, triCentroid.y, triCentroid.z);
			printf("gridPos: %d, %d, %d\n", gridPos.x, gridPos.y, gridPos.z);
			printf("triangleHash[triId]:triangleId[triId] -- ( %d, %d )\n", hashId, triangleId->m_Data[triId]);
		}*/
	}
}

void SpatialHashSystem::sortTrianglesCPU()
{
	thrust::sort_by_key(m_hTriangleHash.m_Data.begin(),
		m_hTriangleHash.m_Data.end(),
		m_hTriangleId.m_Data.begin());
	/*printf("\n");
	for (uint triId = 0; triId < m_numTriangles; ++triId)
	{
		printf("Hash:triangleId -- ( %d, %d )\n", m_hTriangleHash.m_Data[triId], m_hTriangleId.m_Data[triId]);
	}*/
}

// cell range of a hash: [cellStart, cellEnd]
void SpatialHashSystem::findCellStartCPU()
{
	uint index;
	uint hashId;
	uint hashIdPrev;
	uint hashIdNext;
	for (index = 0; index < m_numTriangles - 1; ++index)
	{
		hashId = m_hTriangleHash.m_Data[index];
		hashIdNext = m_hTriangleHash.m_Data[index + 1];
		// printf("hashId %d, hashIdNext %d\n", hashId, hashIdNext);
		if (index == 0)
			m_hCellStart.m_Data[hashId] = index;
		else if (hashId != hashIdNext)
		{
			m_hCellStart.m_Data[hashIdNext] = index + 1;
			m_hCellEnd.m_Data[hashId] = index;
		}
	}
	hashId = m_hTriangleHash.m_Data[index];
	hashIdPrev = m_hTriangleHash.m_Data[index - 1];
	if ((index == m_numTriangles - 1) && (hashId != hashIdPrev))
		m_hCellStart.m_Data[hashId] = index;
	m_hCellEnd.m_Data[hashId] = index;

	// printf("start: %d, end: %d,", m_hCellStart.m_Data[95], m_hCellEnd.m_Data[95]);
}

int3 SpatialHashSystem::calcGridPosCPU(float3 triPos)
{
	int3 gridPos;
	gridPos.x = floor((triPos.x - m_girdStartPos.x) / m_cellSize.x);
	gridPos.y = floor((triPos.y - m_girdStartPos.y) / m_cellSize.y);
	gridPos.z = floor((triPos.z - m_girdStartPos.z) / m_cellSize.z);
	// printf("gridPos (%d, %d, %d) --\n", gridPos.x, gridPos.y, gridPos.z);
	return gridPos;
}

uint SpatialHashSystem::calcGridHashCPU(int3 gridPos)
{
	//gridPos.x = gridPos.x & (m_gridSize.x - 1);  // wrap grid, assumes size is power of 2
	//gridPos.y = gridPos.y & (m_gridSize.y - 1);
	//gridPos.z = gridPos.z & (m_gridSize.z - 1);
	uint hashId = ((gridPos.z * m_gridSize.y) * m_gridSize.x) + (gridPos.y * m_gridSize.x) + gridPos.x;
	// printf("gridPos (%d, %d, %d) -- currHashId: %d\n", gridPos.x, gridPos.y, gridPos.z, hashId);
	return hashId;
}