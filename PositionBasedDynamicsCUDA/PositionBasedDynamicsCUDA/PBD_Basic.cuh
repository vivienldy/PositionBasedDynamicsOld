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
#include "cuda_runtime.h"
#include "helper_math.h"
#include "thrust/sort.h"
#include <device_launch_parameters.h>


using namespace std;

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

template<class T>
class Buffer
{
public:
	Buffer() : m_BufferSize(0)
	{
		m_sName = "";
		m_DevicePtr = nullptr;
	}

	typedef T ThisType;
	vector<T> m_Data;

	template<class S>
	S* GetDevicePtr() { return  (S*)m_DevicePtr; }
	void* GetDevicePtr() { return  m_DevicePtr; }

	void** GetDevicePtrRaw() { return (void**)&m_DevicePtr; }

	int GetSize() { return m_Data.size(); }
	int GetCopySize() { return m_Data.size() * sizeof(T); }

	int GetCapacity() { return m_Data.capacity(); }

	uint2 EvalBlockSize(int nThread) { return make_uint2(ceil((float)GetSize() / nThread), nThread); }

	void Save(std::ofstream& os)
	{
		if (m_sName.size() == 0)
			return;
		if (!os)
			return;
		os << m_sName << "|" << typeid(T).name() << "|";
		for (int i = 0; i < m_Data.size(); i++)
			os << m_Data[i] << ";";
		os << std::endl;
	}

	void SetName(std::string n)
	{
		m_sName = n;
	}

	inline bool LoadToHost()
	{

		if (!m_DevicePtr)
		{
			printf("[ %s ] Device Not Registered!\n", m_sName.c_str());
			return false;
		}
		cudaError_t cudaStatus = cudaMemcpy(m_Data.data(), this->GetDevicePtr(), this->GetCopySize(), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("[ %s ]  cudaMemcpy failed\n", m_sName.c_str());
			return false;
		}
		return true;
	}

	inline bool LoadToDevice(bool printSuccInfo = false)
	{
		if (!m_DevicePtr)
		{
			printf("[ %s ] Device Not Registered!\n", m_sName.c_str());
			return false;
		}
		if (m_iDeviceBufferSize < m_Data.size())
		{
			cudaFree(m_DevicePtr);
			m_DevicePtr = nullptr;
			DeviceMalloc();
		}
		cudaError_t cudaStatus = cudaMemcpy(this->GetDevicePtr(), m_Data.data(), this->GetCopySize(), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf(" [%s] cudaMemcpy failed whit size [%d]\n", m_sName.c_str(), this->GetCopySize());
			return false;
		}

		m_iDeviceBufferSize = m_Data.size();
		if (printSuccInfo)
			printf("Load [%s] To Device Succeeded\n", m_sName.c_str());
		return true;
	}

	inline bool DeviceMalloc()
	{
		if (this->GetCopySize() == 0)
		{
			printf("[%s] Buffer Size = 0 !!! \n", m_sName.c_str());
			return false;
		}
		if (m_DevicePtr != nullptr)
		{
			printf("[%s] DeviceRegistered!\n", m_sName.c_str());
			return false;
		}

		cudaError_t cudaStatus = cudaMalloc(this->GetDevicePtrRaw(), this->GetCopySize());
		if (cudaStatus != cudaSuccess)
		{
			printf("[%s] CudaMalloc Failed\n", m_sName.c_str());
			return false;
		}
		/*else
		{
			printf("Malloc device memory succ use | %12.3f MB | [%3s ]\n",
				(float)GetCopySize() / 1024.0f / 1024.0f, m_sName.c_str());
		}*/
		return true;
	}

	inline bool MallocAndLoadToDevice()
	{
		if (!this->m_DevicePtr)
			this->DeviceMalloc();
		if (!this->LoadToDevice())
			return false;
		return true;
	}

private:
	std::string m_sName;
	void* m_DevicePtr;
	long int m_iDeviceBufferSize;
	long int m_BufferSize;
};

typedef Buffer<int> BufferInt;
typedef Buffer<unsigned int> BufferUInt;
typedef Buffer<float> BufferFloat;
typedef Buffer<float3> BufferVector3f;
typedef Buffer<int2> BufferInt2;


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
};

struct ConstraintPBD
{
	Topology topol;
	BufferFloat restLengthBuffer;  // determine the rest pos structure of constrs
	BufferVector3f prdPBuffer;
	BufferVector3f restPosBuffer;
	BufferVector3f restAngleBuffer;
	//BufferFloat thickness;
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
};

class PBDObject
{
public:
	PBDObject()
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
	void SetConstrOption(uint ct, float* stiffnessSetting);

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
	void initMeshTopol();
	void initMassVel();
	void initPosition(float2 cord);
	void initMeshTopolIndices();
	void initConstr();


	void groundTruthTest();

	void initGPUBuffers();
	void freeGPUBuffers();

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

private:
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

	inline void SaveToplogy(Topology top, std::string path)
	{
		std::ofstream ofs(path);
		if (!ofs.is_open())
			return;
		for (int i = 0; i < top.indices.GetSize(); ++i)
		{
			printf("%d-", top.indices.m_Data[i]);
		}
		SaveBuffer<int>(top.indices, ofs);
		SaveBuffer<float3>(top.posBuffer, ofs);
		SaveBuffer<int2>(top.primList, ofs);

		ofs.flush();
		ofs.close();
	}
}
#endif