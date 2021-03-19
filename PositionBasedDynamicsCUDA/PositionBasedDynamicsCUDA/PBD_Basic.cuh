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


using namespace std;

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
	void Integration(float dt);
	void SetTimer(Timer* timer) { this->m_pbdSolverTimer = timer; }

private:
	// Timer
	Timer* m_pbdSolverTimer;
	//bool m_timerStatus;  // 1 - on; 0 - off
	PBDObject* m_pbdObj;
	HardwareType m_ht;
	void advectCPU(float dt);
	void advectGPU(float dt);
	void projectConstraintCPU(SolverType st, int iterations);
	void projectConstraintGPU(SolverType st, int iterations);
	void integrationCPU(float dt);
	void integrationGPU(float dt);
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
			ofs  << data << endl; 
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

	inline void splitString(const std::string& src, std::vector<std::string>& v, const std::string& split)
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

	inline std::vector<std::string> splitString(const std::string& src, const std::string& split)
	{
		std::vector<std::string> _ret = std::vector<std::string>();
		splitString(src, _ret, split);
		return _ret;
	}

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
			std::vector<std::string> dataMean = splitString(buffer, "|");
			if (dataMean.size() != 3)
				continue;
			std::string name = dataMean[0];
			std::vector<std::string> dataBuffer = splitString(dataMean[2], ";");
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
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
					float3 vertexPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					meshPosBuffer->m_Data.push_back(vertexPos);
				}				
			}
			else if (name == "primList")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
					int2 prim = make_int2(std::stoi(data[0]), std::stoi(data[1]));
					meshPrimList->m_Data.push_back(prim);
				}
			}
			else if (name == "velBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
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
				std::vector<std::string> data = splitString(dataBuffer[0], ",");
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
			std::vector<std::string> dataMean = splitString(buffer, "|");
			if (dataMean.size() != 3)
				continue;
			std::string name = dataMean[0];
			std::vector<std::string> dataBuffer = splitString(dataMean[2], ";");
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
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
					float3 vertexPos = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					constrPosBuffer->m_Data.push_back(vertexPos);
				}
			}
			else if (name == "primList")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
					int2 prim = make_int2(std::stoi(data[0]), std::stoi(data[1]));
					constrPrimList->m_Data.push_back(prim);
				}
			}
			else if (name == "prdPBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
					float3 prdP = make_float3(std::stof(data[0]), std::stof(data[1]), std::stof(data[2]));
					constrPrdpBuffer->m_Data.push_back(prdP);
				}
			}
			else if (name == "restPosBuffer")
			{
				for (int i = 0; i < dataBuffer.size(); ++i)
				{
					std::vector<std::string> data = splitString(dataBuffer[i], ",");
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
			std::string dataLine = splitString(buffer, firstSplit)[1];
			dataStringList = splitString(dataLine, secondSplit);
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
			std::string dataLine = splitString(buffer, firstSplit)[1];
			dataStringList = splitString(dataLine, secondSplit);

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