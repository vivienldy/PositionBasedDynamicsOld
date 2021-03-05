#ifndef PBD_BASIC_H
#define PBD_BASIC_H

#include <iostream>
#include<fstream>
#include<vector>
#include <map>
#include <string>
#include <set>
#include<assert.h>
#include<vector_types.h>
#include<vector_functions.h>
#include "cuda_runtime.h"
#include "helper_math.h"
#include "device_launch_parameters.h"

using namespace std;

#define KERNEL_FUNC __device__ __host__ 
#define FORBIDBITS 64
#define MAXFORBID 18446744073709551615ul

//KERNEL_FUNC float Length(float3 vec) pow(vec.x*vec.x + vec.y *vec.y + vec.z *vec.z, 0.5)

enum ConstraintType
{
	DISTANCE,
	BENDING,
	ANCHOR,
	NUMCONSTR
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
	typedef T ThisType;
	void LoadToGPU();
	void LoadToCPU();
	vector<T> m_Data;
	template<class S>
	S* GetDevicePtr() { return  (S*)m_DevicePtr; }

	int GetSize() { return m_Data.size(); }

	int GetCapacity() { return m_Data.capacity(); }

	uint2 EvalBlockSize(int nThread) { return make_uint2(ceil((float)m_Size / nThread), nThread); }

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
private:
	std::string m_sName;
	void* m_DevicePtr;
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
	BufferFloat restLength;
	BufferVector3f prdPBuffer;
	BufferVector3f velBuffer;
	BufferFloat mass;
	//BufferFloat thickness;
	BufferFloat stiffness;
	BufferInt constraintType;
	// Util Maps
	std::map<int, std::vector<int>> Point2PrimsMap;
	P2P Prim2PrimsMap; 

	BufferInt color; // reserved for edge coloring
	BufferInt prdColor; // reserved for edge coloring
	BufferInt sortedColor;  // reserved for edge coloring
};


class PBDObject
{
public:
	PBDObject()
	{

	}
	PBDObject(float dampingRate, float3 gravity, int resX, int resY, float sizeX, float sizeY) :
		dampingRate(dampingRate), gravity(gravity), resX(resX), resY(resY), sizeX(sizeX), sizeY(sizeY)
	{
		Init();
	}
	~PBDObject()
	{

	}
	//static PBDObject* Create(std::string topologyPath, std::string constraintPath)
	//{
	//	PBDObject* ret = new PBDObject();
	//	//
	//	//
	//	//
	//	return ret;
	//}
	void InitConstr(int numOfConstr, float unitMass, float* stiffnesses);

	void groundTruthTest();

	Topology meshTopol;  // opengl will draw this topol
	ConstraintPBD constrPBDBuffer;
	BufferVector3f restPosBuffer;
	float dampingRate;
	float3 gravity;
	int resX, resY;
	float sizeX, sizeY;
	
	
private:
	int detectConflict;
	void Init();
	void CreatePosition(BufferVector3f& positionBuffer, float2 cord, float sizeX, float sizeY, int resY, int resX);
	void CreateOpenGLIndices(BufferInt& openGLIndices, int resY, int resX);
	void CreateDistanceIndices(BufferInt& indices, int resY, int resX);
	void CreateSingleDistConstr(BufferVector3f& positionBuffer, BufferInt2& primList, BufferInt& indices, float stiffness, float unitMass);

	// init two maps for edge coloring
	void GenePoint2PrimsMap(Topology topol);
	void GenePrim2PrimsMap(Topology topol);

	// TODO: Edge Coloring 
	void AssignColorsCPU();
	void ResolveConflictsCPU();
	void EdgeColoring(int iterations);
	void SortEdgesColors();

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
	void SetTarget(PBDObject* pbdObj) { this->pbdObj = pbdObj; }
	void Advect(float dt, HardwareType ht);
	void ProjectConstraint(HardwareType ht, SolverType st, int iterations);
	void Integration(float dt, HardwareType ht);

private:
	PBDObject* pbdObj;
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
//
// 
//namespace Kernel
//{
//
//	__global__ void projectConstraintD(
//		float3* posBuf,
//		int* induices,
//		int2* sortedPrims,
//		int* typeBuf,
//		int nPrims)
//	{
//		if (primId >= nPrims)
//			return;
//		int type = typeBuf[primId];
//		if (type == ConstraintType::DISTANCE)
//		{
//			projectDistance();
//			posBuf[]
//			posBuf[]
//		}
//		if (type == ConstraintType::Attach)
//		{
//			projectAttach();
//			posBuf[]
//				
//		}
//
//		if (type == ConstraintType::BENDING)
//		{
//			projectAttach();
//			posBuf[]
//			posBuf[]
//			posBuf[]
//			posBuf[]
//
//		}
//	}
//
//	__device__ inline void projectDistance(float3& p0,float3& p1);
//	__device__ inline void projectAttach();
//}

#endif

