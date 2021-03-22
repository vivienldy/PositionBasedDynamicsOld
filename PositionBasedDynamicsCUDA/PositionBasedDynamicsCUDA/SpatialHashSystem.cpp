#include "SpatialHashSystem.h"

SpatialHashSystem::SpatialHashSystem(string filename, uint3 gridSize, HardwareType ht) :
	m_gridSize(gridSize),
	m_ht(ht)
{
	PBD_DEBUG;
	// init m_numTriangles, m_loadSuccess, m_hPosBuffer, m_hTriIdxBuffer
	//readMeshFromTxt(filename);
	//assert(m_loadSuccess);
	m_loadSuccess = IO::ReadClothFromTxt(filename, m_hPosBuffer, m_hTriIdxBuffer, m_hTriangleId, &m_numTriangles);
	assert(m_loadSuccess);
	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f); // should be initialized by the user

	m_cellSize = make_float3(2.0f, 2.0f, 2.0f);   // default value

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));
	m_targetId = 645;
	evalTriCentroids();
	// these are initialized after setting the grid info: m_numGridCells, m_hCellStart, m_hCellEnd
}

SpatialHashSystem::SpatialHashSystem(BufferVector3f& posBuffer, BufferInt& triIdxBuffer, HardwareType ht) :
	m_ht(ht)
{
	// init m_numTriangles, m_loadSuccess, m_hPosBuffer, m_hTriIdxBuffer
	m_hPosBuffer = posBuffer;
	m_hTriIdxBuffer = triIdxBuffer;
	m_numTriangles = m_hTriIdxBuffer.GetSize() / 3;

	m_worldOrigin = make_float3(0.0f, 0.0f, 0.0f); // should be initialized by the user

	m_cellSize = make_float3(2.0f, 2.0f, 2.0f);

	m_hTriangleHash.m_Data.resize(m_numTriangles, 0);
	m_hTriangleId.m_Data.resize(m_numTriangles, 0);
	m_hTriCentroidBuffer.m_Data.resize(m_numTriangles, make_float3(0.0f, 0.0f, 0.0f));
	m_targetId = 0;
	evalTriCentroids();
	// these are initialized after setting the grid info: m_numGridCells, m_hCellStart, m_hCellEnd
}

void SpatialHashSystem::SetGridCenter(float3 worldOrigin)
{
	m_worldOrigin = worldOrigin;
}

void SpatialHashSystem::SetGridSize(uint3 gridSize)
{
	m_gridSize = gridSize;
	m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;
	m_hCellStart.m_Data.resize(m_numGridCells, -1);
	m_hCellEnd.m_Data.resize(m_numGridCells, -1);
	// TODO: update the SH grid
}

void SpatialHashSystem::SetDivision(float3 cellSize)
{
	m_cellSize = cellSize;
	// TODO: update the SH grid
}

void SpatialHashSystem::InitSH()
{
	PBD_DEBUG;
	float shiftX = (m_gridSize.x * m_cellSize.x) / 2;
	float shiftY = (m_gridSize.y * m_cellSize.y) / 2;
	float shiftZ = (m_gridSize.z * m_cellSize.z) / 2;
	m_girdStartPos = m_worldOrigin - make_float3(shiftX, shiftY, shiftZ);
	// printf("grid start pos is %f, %f, %f\n", m_girdStartPos.x, m_girdStartPos.y, m_girdStartPos.z);
	m_initialized = true;
	UpdateSH(m_hPosBuffer);
}

void SpatialHashSystem::UpdateSH(BufferVector3f& prdPBuffer)
{
	m_shsTimer->Tick();
	assert(m_initialized);
	switch (m_ht)
	{
	case CPU:
		m_hPosBuffer = prdPBuffer;
		evalTriCentroids();
		calcHashCPU();
		sortTrianglesCPU();
		findCellStartCPU();
		break;
	case GPU:
		break;
	default:
		break;
	}
	m_shsTimer->Tock();
	PBD_DEBUGTIME(m_shsTimer->GetFuncTime());
	// printf("all out of is %d", m_hAllOutOfGrid);
}

bool SpatialHashSystem::outsideGrid(int3 gridPos)
{
	if (gridPos.x < 0 || gridPos.y < 0 || gridPos.z < 0 ||
		gridPos.x >= m_gridSize.x || gridPos.y >= m_gridSize.y || gridPos.z >= m_gridSize.z)
		return true;
	else
		return false;
}

static std::set<int> maskSet = { 17,21,22,23,26,34,41,47,52,55,57,58,87,92,103,106,107,134,136,137,152,153,162,164,168,170,171,174,192,214,217,228,233,248,251,252,259,274,280,285,307,308,310,317,318,319,320,321,323,373,384,411,412,414,415,419,420,421,422,424,425,432,447,448,454,455,456,466,478,479,488,494,525,541,542,544,596,599,601,606,609,610,611,631,643,644,646,650,658,662,664,666,683,728,729,730,732,742,753,758,759,761,763,787,809,827,836,869,871,874,875,877,878,885,887,902,904,905,908,909,913,914,915,919,921,929,930,936,937,938,939,979,1028,1115,1118,1119,1121,1127,1161,1169,1170,1171,1172,1174,1177,1187,1258,1297,1298,1299,1300,1301,1302,1303,1308,1309,1311,1313,1317,1318,1319,1335,1338,1340,1341,1342,1343,1345,1346,1347,1348,1349,1350,1351,1352,1354,1358,1359,1361,1362,1363,1365,1367,1370,1372,1373,1374,1393,1400,1401,1453,1468,1471,1475,1477,1479,1482,1483,1484,1485,1487,1510,1511,1512,1514,1515,1516,1518,1519,1520,1522,1556,1557,1560,1569,1600,1601,1606,1609,1611,1612,1613,1614,1632,1641,1644,1645,1648,1668,1670,1671,1689,1691,1692,1702,1707,1763,1787,1788,1789,1790,1791,1796,1798,1799,1813,1828,1835,1836,1837,1839,1842,1845,1846,1848,1849,1850,1851,1852,1853,1859,1861,1864,1896,1917,1976,1977,1978,1979,1981,1992,2014,2016,2017,2019,2021,2024,2030,2041,2043,2073,2076,2078,2089,2091,2092,2098,2099,2103,2106,2167,2168,2170,2291,2292,2294,2297,2303,2304,2330,2331,2333,2334,2335,2337,2338,2339,2340,2344,2345,2348,2359,2384,2386,2388,2425,2426,2428,2442,2445,2446,2448,2449,2450,2460,2461,2462,2464,2467,2468,2469,2483,2538,2562,2563,2566,2638,2639,2640,2642,2648,2649,2651,2653,2655,2662,2717,2718,2719,2731,2738,2752,2759,2760,2762,2771,2773,2800,2820,2823,2917,2919,2920,2921,2922,2928,2929,2930,2931,2937,2939,2940,2941,2942,2943,2944,2946,2947,2949,2950,2951,2952,2953,2955,2958,2978,2983,2986,2987,2988,2989,2990,2991,2996,2997,2998,3071,3072,3075,3078,3084,3085,3086,3087,3090,3097,3098,3099,3100,3101,3102,3103,3106,3118,3120,3122,3124,3125,3126,3127,3128,3129,3130,3131,3148,3149,3155,3157,3158,3159,3162,3165,3166,3188,3197,3199,3200,3201,3204,3205,3229,3233,3235,3236,3237,3238,3239,3240,3241,3242,3243,3261,3263,3264,3284,3286,3287,3303,3304,3305,3307,3308,3310,3331,3354,3364,3372,3387,3398,3407,3412,3428,3429,3430,3478,3504,3505,3508,3509,3511,3512,3516,3518,3519,3525,3526,3527,3548,3549,3555,3558,3560,3585,3586,3595,3641,3664,3665,3668,3706,3707,3721,3732,3752,3769,3775,3782,3783,3793,3795,3803,3805,3806,3807,3823,3824,3834,3836,3837,3838,3839,3840,3842,3845,3846,3848,3873,4056,4057,4058,4085,4086,4092,4093,4094,4097,4098,4099,4100,4101,4102,4104,4106,4127,4128,4129,4144,4221,4222,4223,4238,4239,4240,4244,4245,4246,4247,4248,4254,4255,4256,4492,4493,4494,4626,4627,4628,4738,4740,4747,4753,4758,4760,4770,4772,4780,4781,4827,4828,4831,4835,4838,4842,4844,4846,4847,4850,4851,4852,4859,4861,4863,4866,4867,4873,4875,4876,4878,4879,4922,4925,4928,4929,4931,4936,4938,4940,4942,4944,4963,5013,5022,5026,5030,5031,5047,5059,5063,5072,5073,5074,5077,5084,5090,5096,5098,5135,5136,5139,5141,5179,5181,5184,5186,5187,5190,5192,5199,5200,5212,5240,5247,5254,5259,5295,5305,5307,5308,5309,5312,5315,5317,5322,5323,5353,5355,5362,5363,5366,5371,5374,5406,5410,5411,5421,5466,5467,5472,5478,5489,5490,5499,5500,5503,5504,5505,5506,5507,5509,5515,5516,5517,5518,5519,5523,5524,5528,5531,5532,5533,5534,5536,5537,5538,5539,5540,5542,5547,5548,5551,5577,5578,5579,5584,5589,5590,5591,5598,5600,5603,5604,5606,5609,5638,5665,5680,5682,5683,5691,5740,5748,5750,5756,5765,5769,5770,5771,5776,5780,5781,5784,5786,5790,5798,5799,5801,5802,5819,5821,5823,5824,5844,5851,5854,5855,5856,5860,5861,5865,5867,5869,5871,5876,5877,5880,5882,5884,5888,5892,5896,5905,5908,5913,5917,5920,5987,6005,6028,6040,6042,6043,6049,6057,6058,6061,6066,6081,6082,6083,6114,6117,6124,6127,6128,6129,6130,6134,6143,6144,6162,6163,6166,6169,6172,6207,6213,6216,6251,6258,6261,6263,6264,6265,6270,6283,6285,6289,6300,6301,6304,6307,6318,6326,6327,6328,6330,6333,6335,6344,6350,6355,6361,6362,6376,6377,6411,6415,6434,6456,6457,6459,6461,6467,6470,6485,6495,6502,6504,6509,6511,6515,6517,6520,6521,6524,6525,6527,6528,6530,6535,6536,6539,6540,6542,6552,6574,6586,6598,6606,6661,6662,6663,6664,6665,6671,6673,6675,6677,6682,6699,6704,6706,6707,6708,6710,6713,6714,6715,6720,6724,6727,6729,6731,6732,6736,6738,6739,6740,6797,6801,6819,6828,6830,6838,6839,6841,6842,6843,6854,6855,6856,6857,6881,6907,6910,6918,6923,6925,6936,6939,6940,6941,6945,6951,6952,6990,7002,7003,7004,7008,7009,7017,7024,7025,7028,7029,7030,7032,7034,7035,7036,7073,7101,7108,7111,7124,7125,7126,7127,7128,7131,7132,7137,7148,7151,7165,7183,7184,7185,7191,7192,7202,7203,7208,7220,7221,7225,7229,7257,7264,7266,7270,7273,7276,7277,7278,7281,7284,7295,7305,7313,7319,7321,7322,7326,7351,7352,7353,7369,7378,7389,7392,7401,7416,7417,7418,7422,7425,7438,7441,7457,7459,7466,7472,7477,7486,7492,7498,7504,7510,7511,7513,7515,7519 };

void SpatialHashSystem::FindNeighbors(BufferInt& neighbors,  // output: a list of triangle IDs (int)
	int targetTriID)   // input: a triangle ID 
{
	neighbors.m_Data.clear();
	if (m_hAllOutOfGrid)  // the entire object is out of the grid
		return;
	float3 targetTriCentroid = m_hTriCentroidBuffer.m_Data[targetTriID];
	// get address in grid
	int3 gridPos = calcGridPosCPU(targetTriCentroid);
	// printf("\t gridPos: %d %d %d\n", gridPos.x, gridPos.y, gridPos.z);
	
	// deal with triangles that are out of the box
	if (outsideGrid(gridPos))
		return;  // current triangle is outside the SH domain
	
	int cellStart;
	int cellEnd;
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
				int hashId = calcGridHashCPU(neighborPos);
				if (hashId == -1)
					continue;
				cellStart = m_hCellStart.m_Data[hashId];
				if (cellStart != -1)
				{
					//printf("\t hashId: %d\n", hashId);
					cellEnd = m_hCellEnd.m_Data[hashId];
					// cell range of a hash: [cellStart, cellEnd]
					for (int i = cellStart; i <= cellEnd; ++i)
					{
						// std::cout << "cellStart: "  << i << std::endl;
						int neighborId = m_hTriangleId.m_Data[i];

						// For Mask Debug
						if (neighborId != targetTriID && maskSet.count(neighborId) != 0)
							neighbors.m_Data.push_back(neighborId);	
						
						/*if (neighborId != targetTriID)
							neighbors.m_Data.push_back(neighborId);*/
					}
				}
			}
		}
	}
}

// I/O
/*
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
			printf("0 (%f, %f, %f),615 (%f, %f, %f),1569 (%f, %f, %f)\n", m_hPosBuffer.m_Data[0].x, m_hPosBuffer.m_Data[0].y, m_hPosBuffer.m_Data[0].z,
				m_hPosBuffer.m_Data[615].x, m_hPosBuffer.m_Data[615].y, m_hPosBuffer.m_Data[615].z, m_hPosBuffer.m_Data[1569].x, m_hPosBuffer.m_Data[1569].y,
				m_hPosBuffer.m_Data[1569].z);
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
*/

// evaluate the centroid for each triangle
void SpatialHashSystem::evalTriCentroids()
{
	auto triCentroidBuffer = &(m_hTriCentroidBuffer);
	auto posBuffer = &(m_hPosBuffer);
	auto triIdxBuffer = &(m_hTriIdxBuffer);
	float3 currCentroid;
	for (int triId = 0; triId < m_numTriangles; ++triId)
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
	for (int triId = 0; triId < m_numTriangles; ++triId)
	{
		float3 triCentroid = triCentroidBuffer->m_Data[triId];
		// get pos in grid
		int3 gridPos = calcGridPosCPU(triCentroid);
		int hashId = calcGridHashCPU(gridPos);
		if (hashId != -1)
			m_hAllOutOfGrid = false;
		// store grid hash and particle index
		triangleHash->m_Data[triId] = hashId;
		triangleId->m_Data[triId] = triId;  // has been initialized if read mesh from txt
		/*if (triId == m_targetId)
		{
			printf("triCentroidBuffer[triId]: (%f, %f, %f)\n", triCentroid.x, triCentroid.y, triCentroid.z);
			printf("gridPos: %d, %d, %d\n", gridPos.x, gridPos.y, gridPos.z);
			printf("triangleHash[triId]:triangleId[triId] -- ( %d, %d )\n", hashId, triangleId->m_Data[triId]);
		}*/
		// printf("triangleHash[triId]:triangleId[triId] -- ( %d, %d )\n", hashId, triangleId->m_Data[triId]);
	}
}

void SpatialHashSystem::sortTrianglesCPU()
{
	thrust::sort_by_key(m_hTriangleHash.m_Data.begin(),
		m_hTriangleHash.m_Data.end(),
		m_hTriangleId.m_Data.begin());
	/*printf("\n");
	for (int triId = 0; triId < m_numTriangles; ++triId)
	{
		printf("Hash:triangleId -- ( %d, %d )\n", m_hTriangleHash.m_Data[triId], m_hTriangleId.m_Data[triId]);
	}*/
}

void SpatialHashSystem::checkBuffersSizes()
{
	printf("m_hPosBuffer size: %d\n", m_hPosBuffer.GetSize());
	printf("m_hTriIdxBuffer size: %d\n", m_hTriIdxBuffer.GetSize());
	printf("m_hTriCentroidBuffer size: %d\n", m_hTriCentroidBuffer.GetSize());
	printf("m_hTriangleHash size: %d\n", m_hTriangleHash.GetSize());
	printf("m_hTriangleId size: %d\n", m_hTriangleId.GetSize());
	printf("m_hCellStart size: %d\n", m_hCellStart.GetSize());
	printf("m_hCellEnd size: %d\n", m_hCellEnd.GetSize());
}
// cell range of a hash: [cellStart, cellEnd]
void SpatialHashSystem::findCellStartCPU()
{
	int index = 0;
	int hashId = 0;
	int hashIdPrev = 0;
	int hashIdNext = 0;
	//checkBuffersSizes();
	//printf("num of triangles: %d\n", m_numTriangles);
	for (index = 0; index < m_numTriangles - 1; ++index)
	{
		hashId = m_hTriangleHash.m_Data[index];
		hashIdNext = m_hTriangleHash.m_Data[index + 1];

		//printf("index: %d\n", index);
		/*printf("hashId: %d\n", hashId);
		printf(" hashIdNext: %d\n", hashIdNext);*/

		if (index == 0 && hashId != -1)
			m_hCellStart.m_Data[hashId] = index;
		else if (hashId != -1 && hashIdNext != -1 && hashId != hashIdNext)
		{
			//printf("\tNew hashid started by %d\n", (index + 1));
			//printf("\thashIdNext: %d\n", hashIdNext);
			m_hCellStart.m_Data[hashIdNext] = index + 1;
			m_hCellEnd.m_Data[hashId] = index;
		}
		else
		{
			; // should do nothing
		}
	}
	// printf("start: %d, end: %d,", m_hCellStart.m_Data[95], m_hCellEnd.m_Data[95]);
	hashId = m_hTriangleHash.m_Data[index];
	hashIdPrev = m_hTriangleHash.m_Data[index - 1];
	if ((index == m_numTriangles - 1) && (hashId != hashIdPrev) && hashId != -1)
		m_hCellStart.m_Data[hashId] = index;
	if(hashId != -1)
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

int SpatialHashSystem::calcGridHashCPU(int3 gridPos)
{
	//gridPos.x = gridPos.x & (m_gridSize.x - 1);  // wrap grid, assumes size is power of 2
	//gridPos.y = gridPos.y & (m_gridSize.y - 1);
	//gridPos.z = gridPos.z & (m_gridSize.z - 1);
	int hashId;
	if (outsideGrid(gridPos))
		hashId = -1;
	else
		hashId = ((gridPos.z * m_gridSize.y) * m_gridSize.x) + (gridPos.y * m_gridSize.x) + gridPos.x;
	// printf("gridPos (%d, %d, %d) -- currHashId: %d\n", gridPos.x, gridPos.y, gridPos.z, hashId);
	return hashId;
}
