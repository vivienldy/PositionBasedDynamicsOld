#include<vector>
#include <string>
#include "cuda_runtime.h"
# include <thrust/host_vector.h>
#include"PBD_Basic.cuh"
#include"CCD_Basic.h"


using namespace std;

static const int pbdObjTimer = 0,
pbdSolverTimer = 1,
colliSolverTimer = 2,
shsTimer = 3,
globalTimer = 4;

static const int nModules = 5;

void WritePointsToFile(BufferVector3f positionBuffer, int frame)
{
	fstream file;
	string path = "D://0310ContinuousCollisionDectection//0315Test//mG." + to_string(frame) + ".obj";
	file.open(path, ios::out);
	file << "g" << endl;
	for (int i = 0; i < positionBuffer.GetSize(); i++)
	{
		file << "v " << positionBuffer.m_Data[i].x << "  " << positionBuffer.m_Data[i].y << "  " << positionBuffer.m_Data[i].z << endl;
	}
	file.close();
}

void ContinueSim()
{
	HardwareType ht = CPU;
	PBDObject pbdObj;
	pbdObj.ContinueSimInit("D://0319CCDTest//continueSimData//meshTopol//NewLargeClothMeshTopol.19.cache",
										   "D://0319CCDTest//continueSimData//constraint//NewLargeClothConstraint.19.cache", ht);

	Timer timers[nModules];
	int startFrame = 1;
	int endFrame = 50;
	int substep = 1;
	int iteration = 75;
	
	SolverType st = GAUSSSEIDEL;
	pbdObj.SetTimer(&timers[pbdObjTimer]);

	SolverPBD pbdSolver;
	pbdSolver.SetTarget(&pbdObj);
	pbdSolver.SetTimer(&timers[pbdSolverTimer]);


	//float3 cellSize = make_float3(2.5f, 2.5f, 2.5f);
	//float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	//uint3 gridSize = make_uint3(6, 4, 6);

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

	int colPasses = 5;
	// for collision debugging
	for (int i = 0; i < pbdObj.meshTopol.posBuffer.GetSize(); ++i)
	{
		colliSolver.m_resolveTimes.m_Data.push_back(make_int2(i, 0));		
	}
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);
	colliSolver.m_debugFrameID = 0;

	int fps = 24;
	float dt = 1.0 / fps / (float)substep;

	for (size_t i = startFrame; i <= endFrame; i++)
	{
		timers[globalTimer].Tick();  // global timer
		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraint(st, iteration);
			for (int col = 0; col < colPasses; col++)
			{
				colliSolver.CCD_SH();
				colliSolver.CollisionResolve();
			}
			printf("------------------------------------frame: %d-------------------------------", i);
			printf("after prdp: (%f, %f, %f)\n",
					  colliSolver.m_afterColliPrdPBuffer.m_Data[3806].x, colliSolver.m_afterColliPrdPBuffer.m_Data[3806].y, colliSolver.m_afterColliPrdPBuffer.m_Data[3806].z);
			pbdSolver.Integration(dt);
			IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//continueSimData//NewcontiSimData." + to_string((i - 1) * 10 + s + 1) + ".cache");
			cout << "topol saved" << endl;
			if (0)
			{
				
				for (int i = 0; i < colliSolver.m_resolveTimes.GetSize(); ++i)
				{
					if (i == 3086 /*|| i == 356 || i == 387 || i == 985 || i == 1286*/)
					{
						printf("before prdp: (%f, %f, %f)\n",
							colliSolver.m_beforeColliPrdPBuffer.m_Data[i].x, colliSolver.m_beforeColliPrdPBuffer.m_Data[i].y, colliSolver.m_beforeColliPrdPBuffer.m_Data[i].z);
						printf("resolve depth: ");
						for (int j = 0; j < colliSolver.m_resolveDepths[i].GetSize(); ++j)
						{
							printf("%f ", colliSolver.m_resolveDepths[i].m_Data[j]);
						}
						cout << endl;
						printf("id: %d contact fs: %d\n", colliSolver.m_resolveTimes.m_Data[i].x, colliSolver.m_nContact.m_Data[i]);
						printf("id: %d resolveTimes: %d\n", colliSolver.m_resolveTimes.m_Data[i].x, colliSolver.m_resolveTimes.m_Data[i].y);
						printf("after prdp: (%f, %f, %f)\n",
							colliSolver.m_afterColliPrdPBuffer.m_Data[i].x, colliSolver.m_afterColliPrdPBuffer.m_Data[i].y, colliSolver.m_afterColliPrdPBuffer.m_Data[i].z);
						printf("-------------------------------------------------\n");
					}
				}
			}
		}
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}

}

// PBD CCD&SH Test

/// Save everystep
///  |
/// continue sim
///  
int main()
{
	ContinueSim();


/*
	//CCDTestMain();
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	int startFrame = 1;
	int endFrame = 50;
	int substep = 1;
	int iteration = 705;
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;
	float stiffnessSetting[1] = { 1.0f };

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
	colliSolver.SetIterations(7);
	colliSolver.SetAcceStruct(&shs);
	colliSolver.SetTimer(&timers[colliSolverTimer]);

	// for collision debugging
	for (int i = 0; i < pbdObj.meshTopol.posBuffer.GetSize(); ++i)
	{
		colliSolver.m_resolveTimes.m_Data.push_back(make_int2(i, 0));
	}
	colliSolver.m_nContact.m_Data.resize(pbdObj.meshTopol.posBuffer.GetSize(), 0);

	
	int fps = 24;
	float dt = 1.0 / fps / (float)substep;

	string meshPath = "D://0319CCDTest//continueSimData//meshTopol//NewLargeClothMeshTopol.";
	string constrPath = "D://0319CCDTest//continueSimData//constraint//NewLargeClothConstraint.";
	string collisionPath = "D://0319CCDTest//continueSimData//collision//NewLargeClothCollision.";
	int colPasses = 5;
	for (size_t i = startFrame; i <= endFrame; i++)
	{
		timers[globalTimer].Tick();  // global timer
		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraint(st, iteration);
			for (int col = 0; col < colPasses; col++)
			{
				colliSolver.CCD_SH();
				colliSolver.CollisionResolve();
			}
			//colliSolver.CCD_SH(); 
			//colliSolver.CCD_N2();
			//colliSolver.CollisionResolve();
			pbdSolver.Integration(dt);
			string path = to_string((i-1) * substep + s + 1) + ".cache";
			pbdObj.SaveMeshTopol(meshPath + path);
			pbdObj.SaveConstraint(constrPath + path);
			colliSolver.SaveCollision(collisionPath + path);
			//IO::SaveToplogy(pbdObj.meshTopol, "D://0310ContinuousCollisionDectection//ContinuousSimData//5ccdTestLow." + to_string((i-1)*10 + s +1)  + ".cache");
			//cout << "topol saved" << endl;
		}
		IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//1largeClothOutput//NewLargeClothWithSphere." + to_string(i) + ".cache");
		cout << "topol saved" << endl;
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
		//WritePointsToFile(pbdObj.constrPBDBuffer.prdPBuffer, i);	
	}
	*/
	return 0;
}