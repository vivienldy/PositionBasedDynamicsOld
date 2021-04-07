#pragma once
#ifndef CCD_BASIC_H
#define CCD_BASIC_H

#include"PBD_Basic.cuh"
#include "SpatialHashSystem.h"
#include "Utility.h"

#define __ISDUPLICATE 1
#define __USE_EPSILON 1

struct Contact
{
	enum Type { VF, EE, VFDCD } type;
	//bool thicknessVF;
	bool colliFreeDir;
	float t; // time
	double w[4];  // weight
	float3 n; // normal
};

typedef Buffer<Contact> BufferContact;

struct ContactData
{
	BufferContact ctxs;
	BufferInt ctxIndices; // vector<Node*> nodes;
	BufferInt2 ctxStartNum;
	//bool active;
	void Save(std::ofstream& ofs);
};

class CollisionSolver
{
public:
	CollisionSolver()
	{

	}
	~CollisionSolver()
	{

	}
	void SetTarget(PBDObject* pbdObj)
	{
		this->m_pbdObj = pbdObj;
	}
	void SetTargetTest(Topology& topol, BufferVector3f prdPBuffer)
	{
		this->m_topol = topol;
		this->m_prdPBuffer = prdPBuffer;
	}
	void SetAcceStruct(SpatialHashSystem* shs)
	{
		this->m_shs = shs;
	}
	void SetThickness(float thickness)
	{
		this->m_thickness = thickness;
	}
	void SetIterations(int iterations)
	{
		this->m_iterations = iterations;
	}
	BufferVector3f GetPrdPBuffer()
	{
		return this->m_pbdObj->constrPBDBuffer.prdPBuffer;
	}
	ContactData contactData;
	ContactData repulsionContact;
	ContactData testRepulsionContact;
	void CCD_N2();
	void CCD_SH();
	void CCD_SH_Extended();
	void CCD_SH_Repulsion();
	void CCD_SH_Narrow();
	void ResolveRepulsion(float deltaTime);
	void CollisionResolve();
	void CollisionResolveNew(BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, int debug, int iteration, int& debugFrameId);
	void CollisionResolveWithoutItreation();
	void SaveResult();
	void CCD_N2Test();
	void CollisionResolveTest();
	void SetTimer(Timer* timer) { this->m_ccdSolverTimer = timer; }
	void SaveCollision(string path); // for collision debugging
	void SaveContact(string path); // for collision debugging
	void SavePrdPBeforeCCD(string path); // for collision debugging
	void SavePrdPAfterCCD(string path); // for collision debugging

	// for collision debugging
	BufferVector3f beforeColliPrdPBuffer;
	BufferVector3f afterColliPrdPBuffer;
	BufferInt m_nContact;
	std::map<int, BufferFloat> m_resolveDepths;
	int m_debugFrameID;
	BufferInt vfIndices;
	BufferVector3f afterProjPrdpBuffer;

	std::map<int, std::set<int>> m_contactVFMap;
	std::map<int, BufferVector3f> m_contactVHitMap;

private:
	PBDObject* m_pbdObj;
	// Timer
	Timer* m_ccdSolverTimer;
	//bool m_timerStatus;  // 1 - on; 0 - off

	SpatialHashSystem* m_shs;
	int m_iterations;
	float m_thickness;
	BufferInt m_vfResolveTimes;

	// for ccd testing
	Topology m_topol;
	BufferVector3f m_prdPBuffer;

	bool ExtendedVFTest(
		float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
		Contact& contact, int i0, int i1, int i2, int i3);
	bool RepulsionVFTest(
		float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
		Contact& contact, int i0, int i1, int i2, int i3);
	void AddRepulsion(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
		float3 vtxVel, float3 p1Vel, float3 p2Vel, float3 p3Vel,
		float vtxM, float p1M, float p2M, float p3M, float deltaTime,
		Contact& contact, int  i0, int i1, int i2, int i3);
	void AddRepulsion(
		float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
		float3 vtxVel, float3 p1Vel, float3 p2Vel, float3 p3Vel,
		BufferVector3f& repulsionResolveVel, BufferInt& repulsionResolveTime,
		float vtxM, float p1M, float p2M, float p3M, float deltaTime,
		Contact& contact, int  i0, int i1, int i2, int i3);
	bool NarrowVFCCDTest(
		float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
		Contact& contact, int i0, int i1, int i2, int i3);
	void VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
		Contact ctx, int  i0, int i1, int i2, int i3);
	void VFResolveNew(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
		float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
		Contact& contact, int  i0, int i1, int i2, int i3,
		BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, int debug, int iteration);
	bool VFTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
		float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
		Contact& contact, int i0, int i1, int i2, int i3);
	bool VFResolveTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
		float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
		Contact& contact, int i0, int i1, int i2, int i3);
	bool VFDCDTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
		float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
		Contact& contact);
	bool VFCCDTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
		float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
		Contact::Type type, Contact& contact, int i0, int i1, int i2, int i3);
};

class CCDTest
{
public:
	Topology topol;
	BufferVector3f prdPBuffer;
	void PrepareTestData(string topolPath, string prdPath);
	void writePointsToFile()
	{
		fstream file;
		string path = "D://0310ContinuousCollisionDectection//ResolvedPrdP.obj";
		file.open(path, ios::out);
		file << "g" << endl;
		for (int i = 0; i < prdPBuffer.GetSize(); i++)
		{
			file << "v " << prdPBuffer.m_Data[i].x << "  " << prdPBuffer.m_Data[i].y << "  " << prdPBuffer.m_Data[i].z << endl;
		}
		file.close();
	}
private:
	void readBufferFromTxt(string filename);
	void readMeshFromTxt(string filename);
	void splitString(const std::string& src, std::vector<std::string>& v, const std::string& split);
	std::vector<std::string> splitString(const std::string& src, const std::string& split);
	void createPrimList(int count);

};

const double infinity = numeric_limits<double>::infinity();

void CCDTestMain();
int solve_cubic(double a3, double a2, double a1, double a0, double t[3]);
int solve_quadratic(double a, double b, double c, double x[2]);
double newtons_method(double a, double b, double c, double d, double x0, int init_dir);
bool ComputeDir(float3 vtx, float3 pvtx, float3 n);
float Point2Plane(float3 vtx, float3 p1, float3 p2, float3 p3);
float3 pos(const float3 pos, const float3 prdPos, double t);
float3 BarycentricCoord(float3 pos, float3 p0, float3 p1, float3 p2);
bool IsDuplicated(BufferInt2& primList, BufferInt& indices, int i0, int i1, int i2, int i3);

//inline float Distance(float3 p1, float3 p2)
//{
//	return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
//}

inline auto stp(const float3& u, const float3& v, const float3& w) -> decltype(dot(u, cross(v, w))) { return dot(u, cross(v, w)); }
template <typename T> inline  T sgn(const T& x) { return x < 0 ? -1 : 1; }
inline auto norm2(const float3& u)-> decltype(dot(u, u)) { return dot(u, u); }
template <typename T> inline T min(const T& a, const T& b, const T& c) { return std::min(a, std::min(b, c)); }
template <typename T> inline T min(const T& a, const T& b, const T& c, const T& d) { return std::min(std::min(a, b), std::min(c, d)); }

// -------------- Data Orienated -------------------
void CCD_SH(ContactData& contactData, SpatialHashSystem& shs, Topology meshTopol, BufferVector3f prdPBuffer, float thickness);
bool VFTest(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos, float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP, Contact& contact, float thickness, int i0, int i1, int i2, int i3);
bool VFDCDTest(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos, float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP, Contact& contact, float thickness, int i0, int i1, int i2, int i3);
bool VFCCDTest(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos, float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP, Contact::Type type, Contact& contact, float thickness, int i0, int i1, int i2, int i3);
void CollisionResolve(Topology meshTopol, BufferVector3f& prdPBuffer, ContactData contactData, int itereations, float thickness, int& debugFrameId, BufferInt& vfIndices, BufferInt& resolveTimes);
void VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos, float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd, Contact contact, float thickness, BufferInt& resolveTimes, int i0, int i1, int i2, int i3);
void VFResolveNew(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
	Contact contact, float thickness,
	int i0, int i1, int i2, int i3);
bool VFResolveTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3, float thickness);
void CollisionResolveNew(
	Topology meshTopol,
	BufferVector3f& prdPBuffer,
	ContactData contactData,
	int itereations, // for later resolve iteration times
	float thickness,
	int& debugFrameId,
	BufferInt& vfIndices);
bool ExtendedVFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3);
void CCD_SH_Extended(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol,
	BufferVector3f prdPBuffer,
	float thickness);
bool NarrowVFCCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3);
void CCD_SH_Narrow(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol,
	BufferVector3f prdPBuffer,
	float thickness,
	std::map<int, std::set<int>>& contactVFMap,
	std::map<int, BufferVector3f>& contactVHitMap);
void CCD_SH_Narrow(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol,
	BufferVector3f prdPBuffer,
	float thickness);
void SavePrdPBeforeCCD(string path, Topology meshTopol, BufferVector3f prdP);
void SavePrdPAfterCCD(string path, Topology meshTopol, BufferVector3f prdP);
void SaveContact(string path, std::map<int, std::set<int>>& contactVFMap, std::map<int, BufferVector3f>& contactVHitMap);
// --------- repulsion -----------
bool RepulsionVFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness,
	int i0, int i1, int i2, int i3);
void CCD_SH_Repulsion(
	ContactData& repulsionContact,
	SpatialHashSystem& shs,
	Topology meshTopol,
	BufferVector3f& prdPBuffer,
	float thickness);
void ResolveRepulsion(
	Topology meshTopol,
	BufferVector3f& prdPBuffer, BufferVector3f& velBuffer, BufferVector3f& repulsionResolve,
	BufferFloat& massBuffer, BufferInt& repulsionResolveTime,
	ContactData repulsionContact,
	float thickness, float deltaTime);
void AddRepulsion(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
	float3& vtxVel, float3& p1Vel, float3& p2Vel, float3& p3Vel,
	float vtxM, float p1M, float p2M, float p3M, float deltaTime,
	Contact& contact, int  i0, int i1, int i2, int i3);
void AddRepulsion(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
	float3& vtxVel, float3& p1Vel, float3& p2Vel, float3& p3Vel,
	BufferVector3f& repulsionResolveVel, BufferInt& repulsionResolveTime,
	float vtxM, float p1M, float p2M, float p3M, float deltaTime,
	Contact& contact, int  i0, int i1, int i2, int i3);
// for collision test 
void readMeshFromTxt(string filename, Topology& topol);
void readBufferFromTxt(string filename, BufferVector3f& prdPBuffer);
#endif