#include<vector>
#include <string>
#include "cuda_runtime.h"
# include <thrust/host_vector.h>
#include"PBD_Basic.cuh"
#include"CCD_Basic.h"
#include  "BufferDebugModule.h"
#include <sstream>

using namespace std;

static const int pbdObjTimer = 0,
pbdSolverTimer = 1,
colliSolverTimer = 2,
shsTimer = 3,
globalTimer = 4;


int gSubStep = 2;
int gPBDIteration = 100; 
int gCollisionPasses = 5;
int gFPS = 24;
float gDeltaTime = 1.0f / gFPS / (float)gSubStep;
int gCollisionResolves = 7;
int gStartFrame = 1;
int gEndFrame = 50;

static const int nModules = 5;

void ContinueSim()
{
	Timer timers[nModules];
	
	HardwareType ht = CPU;
	PBDObject pbdObj;
	pbdObj.ContinueSimInit("D://0319CCDTest//continueSimData//meshTopol//TestDLargeClothWithSphereMeshTopol.54.cache",
											"D://0319CCDTest//continueSimData//constraint//TestDLargeClothWithSphereConstraint.54.cache", ht);
	SolverType st = GAUSSSEIDEL;
	pbdObj.SetTimer(&timers[pbdObjTimer]);

	SolverPBD pbdSolver;
	pbdSolver.SetTarget(&pbdObj);
	pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	// three cloth with sphere
	//float3 cellSize = make_float3(2.5f, 2.5f, 2.5f);
	//float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	//float3 gridCenter = make_float3(0.0f, -1.0f, 0.0f);
	//int3 gridSize = make_int3(6, 4, 6);
	// 5 cloth with sphere
	//float3 cellSize = make_float3(0.5f, 0.5f, 0.5f);
	//float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	//int3 gridSize = make_int3(30, 10, 30);
	// 1 large cloth with sphere
	float3 cellSize = make_float3(0.35f, 0.35f, 0.35f);
	float3 gridCenter = make_float3(0.0f, -3.0f, 0.0f);
	uint3 gridSize = make_uint3(32, 23, 32);

	// initialize SH
	SpatialHashSystem shs(pbdObj.constrPBDBuffer.prdPBuffer, pbdObj.meshTopol.indices, CPU);
	shs.SetGridCenter(gridCenter);
	shs.SetGridSize(gridSize);
	shs.SetDivision(cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	CollisionSolver colliSolver;
	colliSolver.SetTarget(&pbdObj);
	colliSolver.SetThickness(0.03f);
	colliSolver.SetIterations(7);
	colliSolver.SetAcceStruct(&shs);
	colliSolver.SetTimer(&timers[colliSolverTimer]);

	// for collision debugging
	colliSolver.m_debugFrameID = 1;
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	string meshPath = "D://0324CCDTest//continueSimData//testMesh//TestCLargeClothWithSphereMeshTopol.";
	string constrPath = "D://0324CCDTest//continueSimData//testConstraint//TestCLargeClothWithSphereConstraint.";
	string collisionPath = "D://0324CCDTest//continueSimData//testCollision//TestCLargeClothWithSphereCollision.";

	int substep = gSubStep;
	int fps = 24;
	float dt = 1.0 / fps / (float)substep;
	int contiCookTimes = 0;
	for (size_t i = gStartFrame; i <= 1; i++)
	{
		timers[globalTimer].Tick(); 
		
		BufferVector3f fixedBuffer;
		fixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		BufferVector3f vFixedBuffer;
		vFixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		BufferVector3f fFixedBuffer;
		fFixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraintWithColli(st, gPBDIteration, &colliSolver, fixedBuffer, vFixedBuffer, fFixedBuffer, s);
			//for (int i = 0; i < colliSolver.vfIndices.GetSize(); ++i)
			//{
			//	printf("%d,", colliSolver.vfIndices.m_Data[i]);
			//}
			//printf("\n");
			pbdSolver.Integration(dt);
			IO::SaveToplogy(pbdObj.meshTopol, "D://0324CCDTest//continueSimData//12NewcontiSimData." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "topol saved" << endl;
			//string path = to_string((i - 1) * gSubStep + s + 1) + ".cache";
			//pbdObj.SaveMeshTopol(meshPath + path);
			//pbdObj.SaveConstraint(constrPath + path);
			//colliSolver.SaveCollision(collisionPath + path);
			fixedBuffer.SetName("P");
			IO::SaveBuffer(fixedBuffer, "D://0324CCDTest//continueSimData//12TestCLargeClothWithSphereFixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "fixed buffer saved" << endl;
			vFixedBuffer.SetName("P");
			IO::SaveBuffer(vFixedBuffer, "D://0324CCDTest//continueSimData//12TestCLargeClothWithSphereVFixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "v fixed buffer saved" << endl;
			fFixedBuffer.SetName("P");
			IO::SaveBuffer(fFixedBuffer, "D://0324CCDTest//continueSimData//12TestCLargeClothWithSphereFFixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "f fixed buffer saved" << endl;
			pbdObj.velBuffer.SetName("P");
			IO::SaveBuffer(pbdObj.velBuffer, "D://0324CCDTest//continueSimData//12TestCLargeClothWithSpherevelBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "vel buffer saved" << endl;
		}
		
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}

void CollisionTest(int colliPasses) // for later data oriented continue sim
{
	float thickness = 0.03f;
	int debugFrameId = 1;
	int colliResolveIterations = 2;
	Topology meshTopol; // 上一帧collision free导出的meshtopol
	BufferVector3f prdPBuffer; // 当前帧导出的posBuffer
	readMeshFromTxt("D://0319CCDTest//singleResolveResult//frame14.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//singleResolveResult//frame17.txt", prdPBuffer);
	
	float3 cellSize = make_float3(0.35f, 0.35f, 0.35f);
	float3 gridCenter = make_float3(0.0f, -3.0f, 0.0f);
	uint3 gridSize = make_uint3(32, 23, 32);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	BufferInt vfIndices;
	BufferInt resolveTimes;
	for (int cp = 0; cp < colliPasses; ++cp)
	{
		ContactData contactData;
		CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
		CollisionResolve(meshTopol, prdPBuffer, contactData, colliResolveIterations, thickness, debugFrameId, vfIndices, resolveTimes);
		printf("--------%d collision resolve----------\n", cp);
	}
	for (int i = 0; i < vfIndices.GetSize(); ++i)
	{
		printf("%d,", vfIndices.m_Data[i]);
	}
	printf("\n");
	Topology afterResolveTopol;
	prdPBuffer.SetName("P");
	afterResolveTopol.indices = meshTopol.indices;
	afterResolveTopol.primList = meshTopol.primList;
	afterResolveTopol.posBuffer = prdPBuffer;
	afterResolveTopol.indices.SetName("indices");
	afterResolveTopol.primList.SetName("primList");
	afterResolveTopol.posBuffer.SetName("P");
	IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//singleResolveResult//singleResolveResult." + to_string(17) + ".cache");
	printf("---------------------------resolved topol saved--------------------\n");
	
}

static bool gEnableSelf = true;
void RegularSim();
void RegularSimColliWithProj();
void Cloth2StaticTest();
void Cloth2StaticLargeTest();

int main()
{
	BufferDebugModule::GetInstance()->Load();
	//RegularSim();
	//RegularSimColliWithProj();
	//ContinueSim();
	//CollisionTest(gCollisionPasses);

	//Cloth2StaticTest();
	Cloth2StaticLargeTest();
	return 0;
}

void RegularSim()
{
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	float stiffnessSetting[1] = { 1.0f };
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;

	//string topolFileName = "D://0310ContinuousCollisionDectection//ccdTestData//InitTopolLowHard.txt";
	//string distConstrFileName = "D://0310ContinuousCollisionDectection//ccdTestData//DistanceConstrLowHard.txt";
	//string topolFileName = "D://0310ContinuousCollisionDectection//ccdTestData//5SheetsInitTopol.txt";
	//string distConstrFileName = "D://0310ContinuousCollisionDectection//ccdTestData//5SheetsDistanceConstr.txt";
	string topolFileName = "D://0319CCDTest//1ClothWithSphereTopol.txt";
	string distConstrFileName = "D://0319CCDTest//1ClothWithSphereConstr.txt";

	PBDObject pbdObj(dampingRate, gravity, ht);
	pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	pbdObj.SetTimer(&timers[pbdObjTimer]);
	pbdObj.Init(topolFileName, distConstrFileName);

	SolverPBD pbdSolver;
	pbdSolver.SetTarget(&pbdObj);
	pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	// three cloth with sphere
	//float3 cellSize = make_float3(2.5f, 2.5f, 2.5f);
	//float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	//float3 gridCenter = make_float3(0.0f, -1.0f, 0.0f);
	//int3 gridSize = make_int3(6, 4, 6);
	// 5 cloth with sphere
	//float3 cellSize = make_float3(0.5f, 0.5f, 0.5f);
	//float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	//int3 gridSize = make_int3(30, 10, 30);
	// 1 large cloth with sphere
	float3 cellSize = make_float3(0.35f, 0.35f, 0.35f);
	float3 gridCenter = make_float3(0.0f, -3.0f, 0.0f);
	uint3 gridSize = make_uint3(32, 23, 32);

	// initialize SH
	SpatialHashSystem shs(pbdObj.constrPBDBuffer.prdPBuffer, pbdObj.meshTopol.indices, CPU);
	shs.SetGridCenter(gridCenter);
	shs.SetGridSize(gridSize);
	shs.SetDivision(cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	CollisionSolver colliSolver;
	colliSolver.SetTarget(&pbdObj);
	colliSolver.SetThickness(0.03f);
	colliSolver.SetIterations(gCollisionResolves);
	colliSolver.SetAcceStruct(&shs);
	colliSolver.SetTimer(&timers[colliSolverTimer]);

	// for collision debugging
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	string meshPath = "D://0319CCDTest//continueSimData//meshTopol//5NewLargeClothMeshTopol.";
	string constrPath = "D://0319CCDTest//continueSimData//constraint//5NewLargeClothConstraint.";
	string collisionPath = "D://0319CCDTest//continueSimData//collision//5NewLargeClothCollision.";

	int cookTimes = 0;
	int contiCookTimes = 0;
	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
		//if (i >= 20)
		//{
		   // BufferDebugModule::GetInstance()->PushStack(contiCookTimes, "ENTRY_PRDB_检查输入", pbdObj.constrPBDBuffer.prdPBuffer);
		//}
		timers[globalTimer].Tick();  // global timer
		for (size_t s = 0; s < gSubStep; s++)
		{
			pbdSolver.Advect(gDeltaTime);
			pbdSolver.ProjectConstraint(st, gPBDIteration);
			//if (i >= 20)
			//{
			   // BufferDebugModule::GetInstance()->PushStack( contiCookTimes, "算完_ProjectConstraint", pbdObj.constrPBDBuffer.prdPBuffer);
			//}
			//POST Rendering
			if (gEnableSelf)
			{
				for (int col = 0; col < gCollisionPasses; col++)
				{
					colliSolver.CCD_SH();
					colliSolver.CollisionResolve();
				}
			}

			pbdSolver.Integration(gDeltaTime);
			string path = to_string((i - 1) * gSubStep + s + 1) + ".cache";
			pbdObj.SaveMeshTopol(meshPath + path);
			pbdObj.SaveConstraint(constrPath + path);
			colliSolver.SaveCollision(collisionPath + path);
			//IO::SaveToplogy(pbdObj.meshTopol, "D://0310ContinuousCollisionDectection//ContinuousSimData//5ccdTestLow." + to_string((i-1)*10 + s +1)  + ".cache");
			//cout << "topol saved" << endl;

			cookTimes++;
		}
		IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//1largeClothOutput//5NewLargeClothWithSphere." + to_string(i) + ".cache");
		printf("---------------------------frame %d topol saved--------------------\n", i);
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}


void RegularSimColliWithProj()
{
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	float stiffnessSetting[1] = { 1.0f };
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;

	string topolFileName = "D://0319CCDTest//1ClothWithSphereTopol.txt";
	string distConstrFileName = "D://0319CCDTest//1ClothWithSphereConstr.txt";

	PBDObject pbdObj(dampingRate, gravity, ht);
	pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	pbdObj.SetTimer(&timers[pbdObjTimer]);
	pbdObj.Init(topolFileName, distConstrFileName);

	SolverPBD pbdSolver;
	pbdSolver.SetTarget(&pbdObj);
	pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	// 1 large cloth with sphere
	float3 cellSize = make_float3(0.35f, 0.35f, 0.35f);
	float3 gridCenter = make_float3(0.0f, -3.0f, 0.0f);
	uint3 gridSize = make_uint3(32, 23, 32);

	// initialize SH
	SpatialHashSystem shs(pbdObj.constrPBDBuffer.prdPBuffer, pbdObj.meshTopol.indices, CPU);
	shs.SetGridCenter(gridCenter);
	shs.SetGridSize(gridSize);
	shs.SetDivision(cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	CollisionSolver colliSolver;
	colliSolver.SetTarget(&pbdObj);
	colliSolver.SetThickness(0.03f);
	colliSolver.SetIterations(gCollisionResolves);
	colliSolver.SetAcceStruct(&shs);
	colliSolver.SetTimer(&timers[colliSolverTimer]);

	// for collision debugging
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	string meshPath = "D://0319CCDTest//continueSimData//meshTopol//TestDLargeClothWithSphereMeshTopol.";
	string constrPath = "D://0319CCDTest//continueSimData//constraint//TestDLargeClothWithSphereConstraint.";
	string collisionPath = "D://0319CCDTest//continueSimData//collision//TestDLargeClothWithSphereCollision.";

	int contiCookTimes = 0;
	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
		timers[globalTimer].Tick();  // global timer
		BufferVector3f fixedBuffer;
		fixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		BufferVector3f vFixedBuffer;
		vFixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		BufferVector3f fFixedBuffer;
		fFixedBuffer.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));
		for (size_t s = 0; s < gSubStep; s++)
		{
			pbdSolver.Advect(gDeltaTime);
			pbdSolver.ProjectConstraintWithColli(st, gPBDIteration, &colliSolver, fixedBuffer, vFixedBuffer, fFixedBuffer, i);
			pbdSolver.Integration(gDeltaTime);
			string path = to_string((i - 1) * gSubStep + s + 1) + ".cache";
			pbdObj.SaveMeshTopol(meshPath + path);
			pbdObj.SaveConstraint(constrPath + path);
			//colliSolver.SaveCollision(collisionPath + path);
		}
		IO::SaveToplogy(pbdObj.meshTopol, "D://TestResult//0326NewCCDMethod2//TestDLargeClothWithSphere." + to_string(i) + ".cache");
		printf("---------------------------frame %d topol saved--------------------\n", i);
		fixedBuffer.SetName("P");
		IO::SaveBuffer(fixedBuffer, "D://TestResult//0326NewCCDMethod2//TestDLargeClothWithSphereFixedBuffer." + to_string(i) + ".cache");
		printf("---------------------------frame %d fixedBuffer saved--------------------\n", i);
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}
// --------------- version test data -------------------
void Cloth2StaticTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 10;
	Topology meshTopol; // 上一帧collision free导出的meshtopol
	BufferVector3f prdPBuffer; // 当前帧导出的posBuffer
	readMeshFromTxt("D://VersionTestData//2clothStatic//2clothMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://VersionTestData//2clothStatic//2clothPrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(6.0f, 6.0f, 6.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(2, 1, 2);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	BufferInt vfIndices;
	BufferInt resolveTimes;
	for (int cp = 1; cp <= 10; ++cp)
	{
		ContactData contactData;
		CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result
	
		printf("contact size:%d\n", contactData.ctxs.GetSize());
		std::set<int> idList;
		for (int i = 0; i < contactData.ctxStartNum.GetSize(); ++i)
		{
			int id = contactData.ctxStartNum.m_Data[i].x;
			idList.insert(contactData.ctxIndices.m_Data[id]);
		}
		std::set<int>::iterator it;
		printf("id: ");
		for (it = idList.begin(); it != idList.end(); it++)
		{
			printf("%d -- ", *it);
		}
		printf("\n");
		
		CollisionResolveNew(meshTopol, prdPBuffer, contactData, colliResolveIterations, thickness, debugFrameId, vfIndices);
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://VersionTestData//0327ResolveTest//2cloth//result." + to_string(cp) + ".cache");
		printf("---------------------------resolved topol saved--------------------\n");
	}
}

void Cloth2StaticLargeTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 10;
	Topology meshTopol; // 上一帧collision free导出的meshtopol
	BufferVector3f prdPBuffer; // 当前帧导出的posBuffer
	readMeshFromTxt("D://VersionTestData//2clolthLarger//2clothLargeMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://VersionTestData//2clolthLarger//2clothLargePrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(6.0f, 6.0f, 6.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(5, 1, 5);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	BufferInt vfIndices;
	BufferInt resolveTimes;
	for (int cp = 1; cp <= 20; ++cp)
	{
		ContactData contactData;
		CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result
		
		printf("contact size:%d\n", contactData.ctxs.GetSize());
		std::set<int> idList;
		for (int i = 0; i < contactData.ctxStartNum.GetSize(); ++i)
		{
			int id = contactData.ctxStartNum.m_Data[i].x;
			idList.insert(contactData.ctxIndices.m_Data[id]);
		/*	if (contactData.ctxIndices.m_Data[id] == 729)
			{
				printf("%d %d %d %d\n", contactData.ctxIndices.m_Data[id], contactData.ctxIndices.m_Data[id + 1],
					contactData.ctxIndices.m_Data[id + 2], contactData.ctxIndices.m_Data[id + 3]);
			}*/
		}
		std::set<int>::iterator it;
		printf("id: ");
		for (it = idList.begin(); it != idList.end(); it++)
		{
			printf("%d -- ", *it);
		}
		printf("\n");
		
		CollisionResolveNew(meshTopol, prdPBuffer, contactData, colliResolveIterations, thickness, debugFrameId, vfIndices);
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://VersionTestData//0327ResolveTest//2clothLarge//result." + to_string(cp) + ".cache");
		printf("---------------------------frame %d resolved topol saved--------------------\n", cp);
	}
}

