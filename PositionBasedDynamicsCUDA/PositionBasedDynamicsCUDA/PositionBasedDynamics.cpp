#include<vector>
#include <string>
#include "cuda_runtime.h"
# include <thrust/host_vector.h>
#include <chrono>
#include"PBD_Basic.cuh"

using namespace std;


void WritePointsToFile(BufferVector3f positionBuffer, int frame)
{
	fstream file;
	string path = "D://0301PBDCUDA//Test0308//GPUPrdp." + to_string(frame) + ".obj";
	file.open(path, ios::out);
	file << "g" << endl;
	for (int i = 0; i < positionBuffer.GetSize(); i++)
	{
		file << "v " << positionBuffer.m_Data[i].x << "  " << positionBuffer.m_Data[i].y << "  " << positionBuffer.m_Data[i].z << endl;
	}
	file.close();
}

int main()
{
	auto start = chrono::steady_clock::now();
	clock_t tStart = clock();

	int resY = 256;
	int resX = 256;
	float dampingRate = 0.9f;
	float sizeX = 10.0f;
	float sizeY = 10.0f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	int startFrame = 1;
	int endFrame = 20;
	int substep = 4;
	int iteration = 10;
	HardwareType ht = GPU;
	SolverType st = GAUSSSEIDEL;
	float stiffnessSetting[1] = { 1.0f };

	PBDObject pbdObj(dampingRate, gravity, resX, resY, sizeX, sizeY, ht);
	pbdObj.setConstrOption(DISTANCE | ANCHOR, stiffnessSetting);
	pbdObj.Init();

	SolverPBD solver;
	solver.SetTarget(&pbdObj);

	pbdObj.meshTopol.indices = pbdObj.constrPBDBuffer.topol.indices;
	pbdObj.meshTopol.primList = pbdObj.constrPBDBuffer.topol.primList;

	//IO::SaveToplogy(pbdObj.meshTopol, "D:/GPUtopol.cache");
	//cout << "topol saved" << endl;
	//IO::SaveBuffer(pbdObj.constrPBDBuffer.color, "D:/GPUcolor.cache");
	//cout << "color saved" << endl;

	int fps = 24;
	float dt = 1.0 / fps / (float)substep;
	////int frame = 0
	//solver.Advect(dt);
	//solver.ProjectConstraint(GAUSSSEIDEL, iteration);
	//WritePointsToFile(pbdObj.constrPBDBuffer.prdPBuffer, 0);

	for (size_t i = startFrame; i < endFrame; i++)
	{
		for (size_t s = 0; s < substep; s++)
		{
			solver.Advect(dt);
			solver.ProjectConstraint(st, iteration);
			solver.Integration(dt);
		}
		// WritePointsToFile(pbdObj.meshTopol.posBuffer, i);
	}

	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;

}