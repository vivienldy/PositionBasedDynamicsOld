#include<vector>
#include <string>
#include "cuda_runtime.h"
# include <thrust/host_vector.h>
#include"PBD_Basic.cuh"
#include"CCD_Basic.h"


using namespace std;


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

// PBD CCD&SH Test
int main()
{

	//CCDTestMain();

	//auto start = chrono::steady_clock::now();
	//clock_t tStart = clock();


	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	int startFrame = 1;
	int endFrame = 35;
	int substep = 10;
	int iteration = 10;
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;
	float stiffnessSetting[1] = { 1.0f };

	string topolFileName = "D://0310ContinuousCollisionDectection//ccdTestData//InitTopol.txt";
	string distConstrFileName = "D://0310ContinuousCollisionDectection//ccdTestData//DistanceConstr.txt";

	PBDObject pbdObj(dampingRate, gravity, ht);
	pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	pbdObj.Init(topolFileName, distConstrFileName);

	printf("primList:(%d) \n", pbdObj.meshTopol.primList.GetSize());
	//printf("primList:(%d, %d) \n", pbdObj.meshTopol.primList.m_Data[0].x, pbdObj.meshTopol.primList.m_Data[0].y);

	SolverPBD pbdSolver;
	pbdSolver.SetTarget(&pbdObj);

	float3 cellSize = make_float3(3.0f, 3.0f, 3.0f);
	float3 gridCenter = make_float3(0.0f, 4.0f, 0.0f);
	uint3 gridSize = make_uint3(5, 3, 5);

	// initialize SH
	SpatialHashSystem shs(pbdObj.meshTopol.posBuffer, pbdObj.meshTopol.indices, CPU);
	shs.SetGridCenter(gridCenter);
	shs.SetGridSize(gridSize);
	shs.SetDivision(cellSize);
	shs.InitSH();
	shs.UpdateSH(0.0f);

	CollisionSolver colliSolver;
	colliSolver.SetTarget(&pbdObj);
	colliSolver.SetThickness(0.03f);
	//colliSolver.SetIterations(2);
	colliSolver.SetAcceStruct(&shs);

	//pbdObj.meshTopol.indices = pbdObj.constrPBDBuffer.topol.indices;
	//pbdObj.meshTopol.primList = pbdObj.constrPBDBuffer.topol.primList;

	//IO::SaveToplogy(pbdObj.meshTopol, "D:/3SheetsofCloth.cache");
	//cout << "topol saved" << endl;
	//IO::SaveBuffer(pbdObj.constrPBDBuffer.color, "D:/GPUcolor.cache");
	//cout << "color saved" << endl;

	int fps = 24;
	float dt = 1.0 / fps / (float)substep;

	for (size_t i = startFrame; i <= endFrame; i++)
	{
		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraint(st, iteration);
			colliSolver.CCD_SH(); 
			for (int iteration = 0; iteration < 2; ++iteration)
			{
				shs.UpdateSH(dt);
				colliSolver.CollisionResolve();
				colliSolver.ColliWithShpGrd();
			}
			pbdSolver.Integration(dt);
		}
		//WritePointsToFile(pbdObj.constrPBDBuffer.prdPBuffer, i);
		IO::SaveToplogy(pbdObj.meshTopol, "D://0310ContinuousCollisionDectection//0315Test//mG." + to_string(i) + ".cache");
		cout << "topol saved" << endl;
	}

	//auto end = chrono::steady_clock::now();
	//auto diff = end - start;
	//cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;
	return 0;
}