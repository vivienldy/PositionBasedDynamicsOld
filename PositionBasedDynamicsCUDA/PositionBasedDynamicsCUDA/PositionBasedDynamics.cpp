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


int gSubStep = 1;
int gPBDIteration = 705; //??
int gCollisionPasses = 5;
int gFPS = 24;
float gDeltaTime = 1.0f / gFPS / (float)gSubStep;
int gCollisionResolves = 7;
int gStartFrame = 1;
int gEndFrame = 20;

//

static const int nModules = 5;

void ContinueSim()
{
	Timer timers[nModules];
	
	HardwareType ht = CPU;
	PBDObject pbdObj;
	pbdObj.ContinueSimInit("D://0319CCDTest//continueSimData//meshTopol//NewLargeClothMeshTopol.19.cache",
											"D://0319CCDTest//continueSimData//constraint//NewLargeClothConstraint.19.cache", ht);
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
	colliSolver.m_debugFrameID = 0;
	// for collision debugging
	for (int i = 0; i < pbdObj.meshTopol.posBuffer.GetSize(); ++i)
	{
		colliSolver.m_resolveTimes.m_Data.push_back(make_int2(i, 0));
	}
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	int startFrame = 1;
	int endFrame = 50;
	int substep = 1;
	int iteration = 705;
	int fps = 24;
	float dt = 1.0 / fps / (float)substep;
	int colPasses = 5;
	int contiCookTimes = 0;
	for (size_t i = startFrame; i <= endFrame; i++)
	{
		timers[globalTimer].Tick(); 

		//BufferDebugModule::GetInstance()->CompareStack(contiCookTimes, "ENTRY_PRDB_检查输入", pbdObj.constrPBDBuffer.prdPBuffer,0.0f);

		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraint(st, iteration);
			//BufferDebugModule::GetInstance()->CompareStack(contiCookTimes, "算完_ProjectConstraint", pbdObj.constrPBDBuffer.prdPBuffer, 0.0f);
			for (int col = 0; col < colPasses; col++)
			{
				colliSolver.CCD_SH();
				colliSolver.CollisionResolve();
			}
			printf("------------------------------------frame: %d-------------------------------", i);	
			printf("after prdp: (%f, %f, %f)\n",
					  colliSolver.afterColliPrdPBuffer.m_Data[3806].x, colliSolver.afterColliPrdPBuffer.m_Data[3806].y, colliSolver.afterColliPrdPBuffer.m_Data[3806].z);
			pbdSolver.Integration(dt);
			IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//continueSimData//NewcontiSimData." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "topol saved" << endl;
			contiCookTimes++;
		}
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
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
	 for (int i = 0; i < pbdObj.meshTopol.posBuffer.GetSize(); ++i)
	 {
		 colliSolver.m_resolveTimes.m_Data.push_back(make_int2(i, 0));
	 }
	 colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	 string meshPath = "D://0319CCDTest//continueSimData//meshTopol//NewLargeClothMeshTopol.";
	 string constrPath = "D://0319CCDTest//continueSimData//constraint//NewLargeClothConstraint.";
	 string collisionPath = "D://0319CCDTest//continueSimData//collision//NewLargeClothCollision.";

	 int substep = gSubStep;
	 int cookTimes = 0;
	 int contiCookTimes = 0;
	 for (size_t i = gStartFrame; i <= gEndFrame; i++)
	 {
		 //if (i >= 20)
		 //{
			// BufferDebugModule::GetInstance()->PushStack(contiCookTimes, "ENTRY_PRDB_检查输入", pbdObj.constrPBDBuffer.prdPBuffer);
		 //}
		 timers[globalTimer].Tick();  // global timer
		 for (size_t s = 0; s < substep; s++)
		 {
			 pbdSolver.Advect(gDeltaTime);
			 pbdSolver.ProjectConstraint(st, gPBDIteration);
			 //if (i >= 20)
			 //{
				// BufferDebugModule::GetInstance()->PushStack( contiCookTimes, "算完_ProjectConstraint", pbdObj.constrPBDBuffer.prdPBuffer);
			 //}
			 for (int col = 0; col < gCollisionPasses; col++)
			 {
				 colliSolver.CCD_SH();
				 if (i >= 20)
				 {
					 //BufferDebugModule::GetInstance()->PushStack(contiCookTimes*col, "detectContact", colliSolver.contactData.ctxIndices);
				 }
				 colliSolver.CollisionResolve();
				 //recode
				 //std::stringstream sstr;
				 //sstr << col;
				 //BufferDebugModule::GetInstance()->PushStack(contiCookTimes, std::string("resolve_")+ sstr.str(), pbdObj.constrPBDBuffer.prdPBuffer);

			 }
			 if (i >= 20)
			 {
				 //BufferDebugModule::GetInstance()->PushStack(contiCookTimes, "afterResolvePos", colliSolver.afterColliPrdPBuffer);
				 //contiCookTimes++;
			 }
			 pbdSolver.Integration(gDeltaTime);
			 string path = to_string((i-1) * substep + s + 1) + ".cache";
			 pbdObj.SaveMeshTopol(meshPath + path);
			 pbdObj.SaveConstraint(constrPath + path);
			 colliSolver.SaveCollision(collisionPath + path);
			 //IO::SaveToplogy(pbdObj.meshTopol, "D://0310ContinuousCollisionDectection//ContinuousSimData//5ccdTestLow." + to_string((i-1)*10 + s +1)  + ".cache");
			 //cout << "topol saved" << endl;

			 cookTimes++;
		 }
		 IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//1largeClothOutput//NewLargeClothWithSphere." + to_string(i) + ".cache");
		 printf("---------------------------frame %d topol saved--------------------\n", i);
		 timers[globalTimer].Tock();  // global timer
		 PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());


	 }
}

int main()
{
	BufferDebugModule::GetInstance()->Load();
	RegularSim();
	//ContinueSim();
	
	return 0;
}