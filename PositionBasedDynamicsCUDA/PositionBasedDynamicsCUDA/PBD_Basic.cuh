#pragma once
#ifndef PBD_BASIC_H
#define PBD_BASIC_H

#include <iostream>
#include <fstream>
#include<ostream>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <assert.h>
#include <algorithm>
#include <vector_types.h>
#include <vector_functions.h>
#include <chrono>
#include "thrust/sort.h"
#include"timer.h"
#include <device_launch_parameters.h>
#include "Buffer.cuh"
#include "Utility.h"

using namespace std;

#define __USE_REPULSION 0
#define __DEBUG 1

#ifdef  __DEBUG
#define PBD_DEBUGTIME(t)  cout <<  __FUNCTION__ << "   " << t << "ms" << endl;
#define PBD_DEBUG  cout <<  __FUNCTION__  << endl;
#else
#define PBD_DEBUG
#endif //  __DEBUG

#define KERNEL_FUNC __device__ __host__ 
#define FORBIDBITS 64
#define MAXFORBID 18446744073709551615ul

enum ConstraintType
{
	DISTANCE = 0x001,
	BENDING = 0x010,
	ANCHOR = 0x100
};


enum SolverType
{
	JACOBI,
	GAUSSSEIDEL,
	NUMSOLVER
};

enum HardwareType
{
	GPU,
	CPU,
	NUMHARDWARE
};

enum Color
{
	RED,
	GREEN,
	BLUE,
	PINK,
	NUMCOLOR
};

inline std::ofstream& operator<<(std::ofstream& ofs, float3& val)
{
	ofs << val.x << "," << val.y << "," << val.z;
	return ofs;
}

inline std::ofstream& operator<<(std::ofstream& ofs, int2& val)
{
	ofs << val.x << "," << val.y;
	return ofs;
}

class CollisionSolver;

struct P2P
{
	BufferInt indices;
	BufferInt2 startNumList;  // store start and number of neighbor points
};

struct Topology
{
	BufferInt indices;
	BufferVector3f posBuffer;
	BufferInt2 primList;
	void Save(std::ofstream& ofs);
};

struct ConstraintPBD
{
	Topology topol;
	// determine the rest pos structure of constrs
	BufferVector3f prdPBuffer;
	BufferVector3f restPosBuffer;
	BufferVector3f restAngleBuffer;
	//BufferFloat thickness;
	BufferFloat restLengthBuffer;
	BufferFloat stiffnessBuffer;
	BufferInt constraintType;
	HardwareType ht;  // inialized in PBDObject

	// Util Maps
	std::map<int, std::vector<int>> Point2PrimsMap;
	P2P Prim2PrimsMap;

	// edge coloring
	BufferInt color; // reserved for edge coloring
	BufferInt prdColor; // reserved for edge coloring
	BufferInt sortedPrimId; // reserved for edge coloring
	BufferInt2 colorWorksets;
	int detectConflict;
	int* dFlag; // detectConflicts used for GPU allocation 

	// init three constraints
	void InitDistanceConstr(BufferVector3f& meshPosBuffer, float stiffness, int resY, int resX);
	void InitBendingConstr();
	void InitAnchorConstr(BufferVector3f& meshPosBuffer, float stiffness, int resY);

	void InitDistanceIndices(int resY, int resX);
	void InitDistanceInfo(BufferVector3f& meshPosBuffer, float stiffness);

	// init two maps for edge coloring
	void GenePoint2PrimsMap(Topology topol);
	void GenePrim2PrimsMap(Topology topol);

	// TODO: Edge Coloring 
	void AssignColorsCPU();
	void ResolveConflictsCPU();
	void EdgeColoring(int iterations);
	void EdgeColoringCPU(int iterations);
	void EdgeColoringGPU(int iterations);
	void SortEdgesColors();
	void EvalWorksets();
	void Save(std::ofstream& ofs);
};

class PBDObject
{
public:
	PBDObject()
	{
	}
	PBDObject(float dampingRate, float3 gravity, HardwareType ht) :
		dampingRate(dampingRate), gravity(gravity), ht(ht)
	{
	}
	PBDObject(float dampingRate, float3 gravity, int resX, int resY, float sizeX, float sizeY, HardwareType ht) :
		dampingRate(dampingRate), gravity(gravity), resX(resX), resY(resY), sizeX(sizeX), sizeY(sizeY), ht(ht)
	{

	}
	~PBDObject()
	{
		if (ht == GPU)
			freeGPUBuffers();
	}

	void Init();
	void Init(string topolFileName, string distConstrFileName);
	void ContinueSimInit(string meshTopolPath, string constrPath, HardwareType hardwareType);
	void SetConstrOption(uint ct, float* stiffnessSetting);
	void SetTimer(Timer* timer) { this->m_pbdObjTimer = timer; }

	void InitGPUBuffers();
	void Save(string path);
	void Save(std::ofstream& ofs);
	void SaveMeshTopol(string path);
	void SaveConstraint(string path);
	void Read(string path);
	/// <summary>
	/// file name: frame data
	/// frame:0
	/// funtion: acvect
	/// header: 3 2 1
	/// </summary>

	Topology meshTopol;  // opengl will draw this topol
	ConstraintPBD constrPBDBuffer;
	uint ct;
	BufferVector3f velBuffer;
	BufferFloat massBuffer;
	float dampingRate;
	float3 gravity;
	HardwareType ht;
	int resX, resY;
	float sizeX, sizeY;
	float* stiffnessSetting;  // [DISTANCE, BENDING, ANCHOR]

private:
	// Timer
	Timer* m_pbdObjTimer;
	//bool m_timerStatus;  // 1 - on; 0 - off
	void initMeshTopol();
	void initMassVel();
	void initPosition(float2 cord);
	void initMeshTopolIndices();
	void initConstr();

	void freeGPUBuffers();

	void groundTruthTest();

};

class SolverPBD
{
public:
	SolverPBD()
	{

	}
	~SolverPBD()
	{

	}
	void SetTarget(PBDObject* pbdObj)
	{
		this->m_pbdObj = pbdObj;
		this->m_ht = pbdObj->ht;
	}
	void Advect(float dt);
	void ProjectConstraint(SolverType st, int iterations);
	void ProjectConstraintWithColli(SolverType st, int iterations, CollisionSolver* colliSolver, BufferVector3f& fixedBuffer, BufferVector3f& vFixedBufferint, BufferVector3f& fFixedBuffer, int debug);
	void ProjectConstraintWithColli(
		SolverType st, int iterations, CollisionSolver* colliSolver,
		float deltaTime,
		BufferVector3f& fixedBuffer, BufferVector3f& vFixedBufferint, BufferVector3f& fFixedBuffer, int debug);
	void Integration(float dt);
	void SetTimer(Timer* timer) { this->m_pbdSolverTimer = timer; }

private:
	// Timer
	Timer* m_pbdSolverTimer;
	//bool m_timerStatus;  // 1 - on; 0 - off
	PBDObject* m_pbdObj;
	HardwareType m_ht;

	float3 m_sphereCenter = make_float3(0.0f, 0.0f, 0.0f);
	float m_sphereRadius = 1.0f;
	float m_groundHeight = 0.0f;

	void advectCPU(float dt);
	void advectGPU(float dt);
	void projectConstraintCPU(SolverType st, int iterations);
	void projectConstraintWithColliCPU(SolverType st, int iterations, CollisionSolver* colliSolver, BufferVector3f& fixedBuffer,
		BufferVector3f& vFixedBufferint, BufferVector3f& fFixedBuffer, int debug);
	void projectConstraintWithColliCPU(
		SolverType st, int iterations, CollisionSolver* colliSolver,
		float deltaTime,
		BufferVector3f& fixedBuffer,
		BufferVector3f& vFixedBufferint, BufferVector3f& fFixedBuffer, int debug);
	void projectConstraintGPU(SolverType st, int iterations);
	void integrationCPU(float dt);
	void integrationGPU(float dt);
	void ColliWithShpGrd();
	bool ColliderSphere(float3 pointPos, float3 sphereOrigin, float r);
	bool CollideGround(float3 pointPos, float groundHeight);
	float3 GenerateMoveVectorSphere(float3 sphereOrigin, float sphereRadius, float3  p);

};

namespace IO
{
	template <typename S>
	inline void SaveBuffer(Buffer<S>& buffer, std::ofstream& ofs)
	{
		buffer.Save(ofs);
	}

	template <typename S>
	inline void SaveBuffer(Buffer<S>& buffer, std::string path)
	{
		std::ofstream ofs(path);
		if (!ofs.is_open())
			return;
		buffer.Save(ofs);
	}

	template <typename T>
	inline void SaveData(T& data, string name, std::ofstream& ofs)
	{
		if (!ofs)
			return;
		auto type = typeid(T).name();
		ofs << name << "|";
		ofs.write(type, strlen(type));
		ofs << "|";
		if (0 == strcmp(type, "float * __ptr64"))
		{
			ofs << endl;
		}
		else
		{
			ofs << data << endl;
		}
		//cout << "typeName: " << typeid(T).name() << endl;
	}

	inline void SaveToplogy(Topology top, std::string path)
	{
		std::ofstream ofs(path);
		if (!ofs.is_open())
			return;
		//for (int i = 0; i < top.indices.GetSize(); ++i)
		//{
		//	printf("%d-", top.indices.m_Data[i]);
		//}
		top.indices.SetName("Indices");
		top.posBuffer.SetName("P");
		top.primList.SetName("primList");
		SaveBuffer<int>(top.indices, ofs);
		SaveBuffer<float3>(top.posBuffer, ofs);
		SaveBuffer<int2>(top.primList, ofs);

		ofs.flush();
		ofs.close();
	}


	static std::set<int> maskSet = { 17,21,22,23,26,34,41,47,52,55,57,58,87,92,103,106,107,134,136,137,152,153,162,164,168,170,171,174,192,214,217,228,233,248,251,252,259,274,280,285,307,308,310,317,318,319,320,321,323,373,384,411,412,414,415,419,420,421,422,424,425,432,447,448,454,455,456,466,478,479,488,494,525,541,542,544,596,599,601,606,609,610,611,631,643,644,646,650,658,662,664,666,683,728,729,730,732,742,753,758,759,761,763,787,809,827,836,869,871,874,875,877,878,885,887,902,904,905,908,909,913,914,915,919,921,929,930,936,937,938,939,979,1028,1115,1118,1119,1121,1127,1161,1169,1170,1171,1172,1174,1177,1187,1258,1297,1298,1299,1300,1301,1302,1303,1308,1309,1311,1313,1317,1318,1319,1335,1338,1340,1341,1342,1343,1345,1346,1347,1348,1349,1350,1351,1352,1354,1358,1359,1361,1362,1363,1365,1367,1370,1372,1373,1374,1393,1400,1401,1453,1468,1471,1475,1477,1479,1482,1483,1484,1485,1487,1510,1511,1512,1514,1515,1516,1518,1519,1520,1522,1556,1557,1560,1569,1600,1601,1606,1609,1611,1612,1613,1614,1632,1641,1644,1645,1648,1668,1670,1671,1689,1691,1692,1702,1707,1763,1787,1788,1789,1790,1791,1796,1798,1799,1813,1828,1835,1836,1837,1839,1842,1845,1846,1848,1849,1850,1851,1852,1853,1859,1861,1864,1896,1917,1976,1977,1978,1979,1981,1992,2014,2016,2017,2019,2021,2024,2030,2041,2043,2073,2076,2078,2089,2091,2092,2098,2099,2103,2106,2167,2168,2170,2291,2292,2294,2297,2303,2304,2330,2331,2333,2334,2335,2337,2338,2339,2340,2344,2345,2348,2359,2384,2386,2388,2425,2426,2428,2442,2445,2446,2448,2449,2450,2460,2461,2462,2464,2467,2468,2469,2483,2538,2562,2563,2566,2638,2639,2640,2642,2648,2649,2651,2653,2655,2662,2717,2718,2719,2731,2738,2752,2759,2760,2762,2771,2773,2800,2820,2823,2917,2919,2920,2921,2922,2928,2929,2930,2931,2937,2939,2940,2941,2942,2943,2944,2946,2947,2949,2950,2951,2952,2953,2955,2958,2978,2983,2986,2987,2988,2989,2990,2991,2996,2997,2998,3071,3072,3075,3078,3084,3085,3086,3087,3090,3097,3098,3099,3100,3101,3102,3103,3106,3118,3120,3122,3124,3125,3126,3127,3128,3129,3130,3131,3148,3149,3155,3157,3158,3159,3162,3165,3166,3188,3197,3199,3200,3201,3204,3205,3229,3233,3235,3236,3237,3238,3239,3240,3241,3242,3243,3261,3263,3264,3284,3286,3287,3303,3304,3305,3307,3308,3310,3331,3354,3364,3372,3387,3398,3407,3412,3428,3429,3430,3478,3504,3505,3508,3509,3511,3512,3516,3518,3519,3525,3526,3527,3548,3549,3555,3558,3560,3585,3586,3595,3641,3664,3665,3668,3706,3707,3721,3732,3752,3769,3775,3782,3783,3793,3795,3803,3805,3806,3807,3823,3824,3834,3836,3837,3838,3839,3840,3842,3845,3846,3848,3873,4056,4057,4058,4085,4086,4092,4093,4094,4097,4098,4099,4100,4101,4102,4104,4106,4127,4128,4129,4144,4221,4222,4223,4238,4239,4240,4244,4245,4246,4247,4248,4254,4255,4256,4492,4493,4494,4626,4627,4628,4738,4740,4747,4753,4758,4760,4770,4772,4780,4781,4827,4828,4831,4835,4838,4842,4844,4846,4847,4850,4851,4852,4859,4861,4863,4866,4867,4873,4875,4876,4878,4879,4922,4925,4928,4929,4931,4936,4938,4940,4942,4944,4963,5013,5022,5026,5030,5031,5047,5059,5063,5072,5073,5074,5077,5084,5090,5096,5098,5135,5136,5139,5141,5179,5181,5184,5186,5187,5190,5192,5199,5200,5212,5240,5247,5254,5259,5295,5305,5307,5308,5309,5312,5315,5317,5322,5323,5353,5355,5362,5363,5366,5371,5374,5406,5410,5411,5421,5466,5467,5472,5478,5489,5490,5499,5500,5503,5504,5505,5506,5507,5509,5515,5516,5517,5518,5519,5523,5524,5528,5531,5532,5533,5534,5536,5537,5538,5539,5540,5542,5547,5548,5551,5577,5578,5579,5584,5589,5590,5591,5598,5600,5603,5604,5606,5609,5638,5665,5680,5682,5683,5691,5740,5748,5750,5756,5765,5769,5770,5771,5776,5780,5781,5784,5786,5790,5798,5799,5801,5802,5819,5821,5823,5824,5844,5851,5854,5855,5856,5860,5861,5865,5867,5869,5871,5876,5877,5880,5882,5884,5888,5892,5896,5905,5908,5913,5917,5920,5987,6005,6028,6040,6042,6043,6049,6057,6058,6061,6066,6081,6082,6083,6114,6117,6124,6127,6128,6129,6130,6134,6143,6144,6162,6163,6166,6169,6172,6207,6213,6216,6251,6258,6261,6263,6264,6265,6270,6283,6285,6289,6300,6301,6304,6307,6318,6326,6327,6328,6330,6333,6335,6344,6350,6355,6361,6362,6376,6377,6411,6415,6434,6456,6457,6459,6461,6467,6470,6485,6495,6502,6504,6509,6511,6515,6517,6520,6521,6524,6525,6527,6528,6530,6535,6536,6539,6540,6542,6552,6574,6586,6598,6606,6661,6662,6663,6664,6665,6671,6673,6675,6677,6682,6699,6704,6706,6707,6708,6710,6713,6714,6715,6720,6724,6727,6729,6731,6732,6736,6738,6739,6740,6797,6801,6819,6828,6830,6838,6839,6841,6842,6843,6854,6855,6856,6857,6881,6907,6910,6918,6923,6925,6936,6939,6940,6941,6945,6951,6952,6990,7002,7003,7004,7008,7009,7017,7024,7025,7028,7029,7030,7032,7034,7035,7036,7073,7101,7108,7111,7124,7125,7126,7127,7128,7131,7132,7137,7148,7151,7165,7183,7184,7185,7191,7192,7202,7203,7208,7220,7221,7225,7229,7257,7264,7266,7270,7273,7276,7277,7278,7281,7284,7295,7305,7313,7319,7321,7322,7326,7351,7352,7353,7369,7378,7389,7392,7401,7416,7417,7418,7422,7425,7438,7441,7457,7459,7466,7472,7477,7486,7492,7498,7504,7510,7511,7513,7515,7519 };

	inline float Distance(float3 p1, float3 p2)
	{
		return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
	}

	inline bool ReadTopolFromCache(string path, PBDObject* pbdObj)
	{
		std::fstream in;
		in.open(path, std::ios::in);
		if (!in.is_open())
		{
			printf("Error opening the file\n");
			exit(1);
		}

		auto meshPosBuffer = &(pbdObj->meshTopol.posBuffer);
		auto meshIndices = &(pbdObj->meshTopol.indices);
		auto meshPrimList = &(pbdObj->meshTopol.primList);
		auto massBuffer = &(pbdObj->massBuffer);
		auto velBuffer = &(pbdObj->velBuffer);

		bool m_loadSuccess = false;
		std::vector<std::string> dataStringList;


		std::string buffer;
		while (std::getline(in, buffer))
		{
			std::vector<std::string> dataMean = UT::splitString(buffer, "|");
			if (dataMean.size() != 3)
				continue;
			std::string name = dataMean[0];
			std::vector<std::string> dataBuffer = UT::splitString(dataMean[2], ";");
			if (name == "Indices")
			{
				for (int i = 0; i < dataBuffer.size(); i++)
				{
					meshIndices->m_Data.push_back(std::stoi(dataBuffer[i]));
				}
			}
			else if (name == "P")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					float3 vertexPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					meshPosBuffer->m_Data.push_back(vertexPos);
				}
			}
			else if (name == "primList")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					int2 prim = make_int2(std::stoi(data[0]), std::stoi(data[1]));
					meshPrimList->m_Data.push_back(prim);
					/*
					bool ignore = false;
					for (int i = 0; i < prim.y; ++i)
					{
						int currId = meshIndices->m_Data[prim.x + i];
						if (maskSet.count(currId)==0)
						{
							ignore = true;
							break;
						}
					}
					if(!ignore)
						meshPrimList->m_Data.push_back(prim);
					*/
				}
				std::cout << "meshPrimList:" << meshPrimList->GetSize() << std::endl;
			}
			else if (name == "velBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					float3 velPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					velBuffer->m_Data.push_back(velPos);
				}
			}
			else if (name == "massBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					massBuffer->m_Data.push_back(std::stof(dataBuffer[i]));
				}
			}
			else if (name == "dampingRate")
			{
				pbdObj->dampingRate = std::stof(dataBuffer[0]);
			}
			else if (name == "gravity")
			{
				std::vector<std::string> data = UT::splitString(dataBuffer[0], ",");
				float3 g = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
				pbdObj->gravity = g;
			}
		}
		in.close();
		m_loadSuccess = true;
		//printf("Read Topol File Done!\n");
		return m_loadSuccess;
	}

	inline bool ReadConstraintFromCache(string path, PBDObject* pbdObj)
	{
		std::fstream in;
		in.open(path, std::ios::in);
		if (!in.is_open())
		{
			printf("Error opening the file\n");
			exit(1);
		}

		auto constrPosBuffer = &(pbdObj->constrPBDBuffer.topol.posBuffer);
		auto constrIndices = &(pbdObj->constrPBDBuffer.topol.indices);
		auto constrPrimList = &(pbdObj->constrPBDBuffer.topol.primList);
		auto constrPrdpBuffer = &(pbdObj->constrPBDBuffer.prdPBuffer);
		//auto restPosBuffer = &(pbdObj->constrPBDBuffer.restPosBuffer); // now is empty
		auto restLengthBuffer = &(pbdObj->constrPBDBuffer.restLengthBuffer);
		auto stiffnessBuffer = &(pbdObj->constrPBDBuffer.stiffnessBuffer);
		auto constraintType = &(pbdObj->constrPBDBuffer.constraintType);

		bool m_loadSuccess = false;
		std::vector<std::string> dataStringList;


		std::string buffer;
		while (std::getline(in, buffer))
		{
			std::vector<std::string> dataMean = UT::splitString(buffer, "|");
			if (dataMean.size() != 3)
				continue;
			std::string name = dataMean[0];
			std::vector<std::string> dataBuffer = UT::splitString(dataMean[2], ";");
			if (name == "Indices")
			{
				for (int i = 0; i < dataBuffer.size(); i++)
				{
					constrIndices->m_Data.push_back(std::stoi(dataBuffer[i]));
				}
			}
			else if (name == "P")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					float3 vertexPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					constrPosBuffer->m_Data.push_back(vertexPos);
				}
			}
			else if (name == "primList")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					int2 prim = make_int2(std::stoi(data[0]), std::stoi(data[1]));
					constrPrimList->m_Data.push_back(prim);
				}
			}
			else if (name == "prdPBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					float3 prdP = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					constrPrdpBuffer->m_Data.push_back(prdP);
				}
			}
			else if (name == "restPosBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = UT::splitString(dataBuffer[i], ",");
					float3 restPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					constrPrdpBuffer->m_Data.push_back(restPos);
				}
			}
			else if (name == "restLengthBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					float l = std::stof(dataBuffer[i]);
					restLengthBuffer->m_Data.push_back(l);
				}
			}
			else if (name == "stiffnessBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					float s = std::stof(dataBuffer[i]);
					stiffnessBuffer->m_Data.push_back(s);
				}
			}
			else if (name == "constraintType")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					float c = std::stoi(dataBuffer[i]);
					constraintType->m_Data.push_back(c);
				}
			}
		}
		in.close();
		m_loadSuccess = true;
		//printf("Read Topol File Done!\n");
		return m_loadSuccess;
	}

	inline bool ReadTopolFromTxt(string filename, PBDObject* pbdObj)
	{
		//printf("Read Topol From Txt\n");
		std::fstream in;
		in.open(filename, std::ios::in);
		if (!in.is_open()) {
			printf("Error opening the file\n");
			exit(1);
		}

		bool m_loadSuccess = false;
		std::vector<std::string> dataStringList;

		auto meshPosBuffer = &(pbdObj->meshTopol.posBuffer);
		auto meshIndices = &(pbdObj->meshTopol.indices);
		auto meshPrimList = &(pbdObj->meshTopol.primList);
		auto massBuffer = &(pbdObj->massBuffer);
		auto velBuffer = &(pbdObj->velBuffer);

		meshIndices->SetName("Indices");
		meshPosBuffer->SetName("P");
		meshPrimList->SetName("primList");

		float3 initVel = make_float3(0.0f, 0.0f, 0.0f);


		std::string buffer;
		int lineNum = 0;
		while (std::getline(in, buffer)) {
			// string to char *
			const std::string firstSplit = ":";
			const std::string secondSplit = ",";
			std::string dataLine = UT::splitString(buffer, firstSplit)[1];
			dataStringList = UT::splitString(dataLine, secondSplit);
			if (lineNum == 0)  // first line of the input file: position
			{
				//assert(dataStringList.size() % 3 == 0);
				for (uint i = 0; i < dataStringList.size(); i += 3)
				{
					// std::cout << dataStringList[i] << std::endl;
					float3 vertexPos = make_float3(std::stof(dataStringList[i]),
						std::stof(dataStringList[i + 1]),
						std::stof(dataStringList[i + 2]));
					// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
					meshPosBuffer->m_Data.push_back(vertexPos);
					massBuffer->m_Data.push_back(1.0);
					velBuffer->m_Data.push_back(initVel);
				}
				//std::cout << "m_hPosBuffer: " << m_hPosBuffer.GetSize() << endl;
				//assert(meshPosBuffer->GetSize() == dataStringList.size() / 3);
			}
			else  // second line of the input file: vertices tet
			{
				//assert(dataStringList.size() % 3 == 0);
				for (uint i = 0; i < dataStringList.size(); i += 3)
				{
					meshIndices->m_Data.push_back(std::stoi(dataStringList[i]));
					meshIndices->m_Data.push_back(std::stoi(dataStringList[i + 1]));
					meshIndices->m_Data.push_back(std::stoi(dataStringList[i + 2]));
					meshPrimList->m_Data.push_back(make_int2(i, 3));
				}
				//std::cout << "m_hTriangleId: " << m_hTriangleId.GetSize() << endl;
			}
			++lineNum;
		}
		in.close();
		m_loadSuccess = true;
		//printf("Read Topol File Done!\n");
		return m_loadSuccess;
	}

	// TODO: more info for initializing a distance constr
	// read file for 3 sheets of cloth distance constraint
	inline bool ReadDistConstrFromTxt(string filename, PBDObject* pbdObj)
	{
		//printf("Read Dist Constr From Txt\n");
		std::fstream in;
		in.open(filename, std::ios::in);
		if (!in.is_open()) {
			printf("Error opening the file\n");
			exit(1);
		}

		bool m_loadSuccess = false;
		std::vector<std::string> dataStringList;

		auto constrPBDBuffer = &(pbdObj->constrPBDBuffer);
		auto topolIndicies = &(constrPBDBuffer->topol.indices);
		auto topolPrimBuffer = &(constrPBDBuffer->topol.primList);
		auto sortedPrimId = &(constrPBDBuffer->sortedPrimId);
		auto stiffnessBuffer = &(constrPBDBuffer->stiffnessBuffer);
		auto restLengthBuffer = &(constrPBDBuffer->restLengthBuffer);
		auto constraintType = &(constrPBDBuffer->constraintType);
		auto color = &(constrPBDBuffer->color);
		auto prdColor = &(constrPBDBuffer->prdColor);
		auto meshPosBuffer = &(pbdObj->meshTopol.posBuffer);
		auto stiffness = pbdObj->stiffnessSetting[0];

		constrPBDBuffer->ht = pbdObj->ht;
		constrPBDBuffer->topol.posBuffer = pbdObj->meshTopol.posBuffer;
		constrPBDBuffer->prdPBuffer = pbdObj->meshTopol.posBuffer;

		topolPrimBuffer->SetName("primList");
		topolIndicies->SetName("Indices");
		color->SetName("color");
		prdColor->SetName("prdcolor");
		sortedPrimId->SetName("sortedPrimId");

		std::string buffer;
		int lineNum = 0;
		while (std::getline(in, buffer))
		{
			if (lineNum == 0)  // first line of the input file: position
			{
				++lineNum;
				continue; // don't need the first line
			}
			// string to char *
			const std::string firstSplit = ":";
			const std::string secondSplit = ",";
			std::string dataLine = UT::splitString(buffer, firstSplit)[1];
			dataStringList = UT::splitString(dataLine, secondSplit);

			// assert(dataStringList.size() % 2 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 2)
			{
				int i0 = std::stoi(dataStringList[i]);
				int i1 = std::stoi(dataStringList[i + 1]);
				// init primList & sortedPrimID 
				topolIndicies->m_Data.push_back(i0);
				topolIndicies->m_Data.push_back(i1);
				topolPrimBuffer->m_Data.push_back(make_int2(i, 2));
				sortedPrimId->m_Data.push_back(i);
				// init stiffness
				stiffnessBuffer->m_Data.push_back(stiffness);
				// init rest length
				float d = IO::Distance(meshPosBuffer->m_Data[i0], meshPosBuffer->m_Data[i1]);
				restLengthBuffer->m_Data.push_back(d);
				// init contraint type
				constraintType->m_Data.push_back(DISTANCE);
				// init color
				color->m_Data.push_back(-1);
				prdColor->m_Data.push_back(-1);
			}
		}
		in.close();
		//printf("Read Constr File Done!\n");
		constrPBDBuffer->GenePoint2PrimsMap(constrPBDBuffer->topol);
		constrPBDBuffer->GenePrim2PrimsMap(constrPBDBuffer->topol);
		if (pbdObj->ht == GPU)
		{
			pbdObj->InitGPUBuffers();  // changed this to a public method
			constrPBDBuffer->EdgeColoring(20000);
		}
		m_loadSuccess = true;
		return m_loadSuccess;
	}

	// read file for SHS files
	inline bool ReadClothFromTxt(string filename, BufferVector3f& posBuffer, BufferInt& idxBuffer, BufferInt& triIdBuffer, uint* numTriangles)
	{
		std::fstream in;
		in.open(filename, std::ios::in);
		if (!in.is_open()) {
			printf("Error opening the file\n");
			exit(1);
		}

		bool m_loadSuccess = false;
		std::vector<std::string> dataStringList;

		posBuffer.m_Data.clear();
		idxBuffer.m_Data.clear();
		triIdBuffer.m_Data.clear();

		std::string buffer;
		int lineNum = 0;
		while (std::getline(in, buffer)) {
			// string to char *
			const std::string firstSplit = ":";
			const std::string secondSplit = ",";
			std::string dataLine = UT::splitString(buffer, firstSplit)[1];
			dataStringList = UT::splitString(dataLine, secondSplit);
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
					posBuffer.m_Data.push_back(vertexPos);
				}
				//std::cout << "m_hPosBuffer: " << m_hPosBuffer.GetSize() << endl;
				assert(posBuffer.GetSize() == dataStringList.size() / 3);
				/*printf("0 (%f, %f, %f),615 (%f, %f, %f),1569 (%f, %f, %f)\n", m_hPosBuffer.m_Data[0].x, m_hPosBuffer.m_Data[0].y, m_hPosBuffer.m_Data[0].z,
					m_hPosBuffer.m_Data[615].x, m_hPosBuffer.m_Data[615].y, m_hPosBuffer.m_Data[615].z, m_hPosBuffer.m_Data[1569].x, m_hPosBuffer.m_Data[1569].y,
					m_hPosBuffer.m_Data[1569].z);*/
			}
			else  // second line of the input file: vertices tet
			{
				assert(dataStringList.size() % 3 == 0);
				for (uint i = 0; i < dataStringList.size(); i += 3)
				{
					idxBuffer.m_Data.push_back(std::stoi(dataStringList[i]));
					idxBuffer.m_Data.push_back(std::stoi(dataStringList[i + 1]));
					idxBuffer.m_Data.push_back(std::stoi(dataStringList[i + 2]));
					triIdBuffer.m_Data.push_back((i / 3));
				}
				//std::cout << "m_hTriangleId: " << m_hTriangleId.GetSize() << endl;
				*numTriangles = triIdBuffer.GetSize();
				assert(triIdBuffer.GetSize() == dataStringList.size() / 3);
			}
			++lineNum;
		}
		in.close();
		m_loadSuccess = true;
		//printf("Read Cloth File Done!\n");
		return m_loadSuccess;
	}
}
#endif