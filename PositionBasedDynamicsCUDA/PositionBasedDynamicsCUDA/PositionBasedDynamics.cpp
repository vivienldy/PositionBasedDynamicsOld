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

static const int nModules = 5;

<<<<<<< Updated upstream
int gSubStep = 1;
int gPBDIteration = 705; //??
=======
int gSubStep = 2;
int gPBDIteration = 100;
>>>>>>> Stashed changes
int gCollisionPasses = 5;
int gFPS = 24;
float gDeltaTime = 1.0f / gFPS / (float)gSubStep;
int gCollisionResolves = 7;
int gStartFrame = 1;
int gEndFrame = 20;

static bool gEnableSelf = true;
void RegularSim();
void RegularSimColliWithProj();
void ContinueSim();
void Cloth2StaticTest();
void Cloth2ChildStaticTest();
void Cloth2StaticLargeTest();
void EETest1();  // two simple triangles
void EETest2(); // two mirror/parallel triangles
void EETest3(); // many triangles
void EETest4();  // two small clothes thickness test
void EETest5();  // two small clothes penetration test
void CollisionTest(int colliPasses); // just for collision test
bool EESingleTest(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, float* constrVal, Contact& contact, float thickness);
void EECheck();

int main()
{
	//BufferDebugModule::GetInstance()->Load();
	//RegularSim();
	//CollisionTest(gCollisionPasses);
    //ContinueSim();
	RegularSimColliWithProj();
	//Cloth2StaticTest();
	//Cloth2ChildStaticTest();
	//Cloth3StaticTest();
	//EECheck();

	// EE Test 1
	//EETest1();
	// EE Test 2
	//EETest2();
	// EE Test3
	//Cloth2StaticLargeTest();
	// EE Test4
	//EETest4();
	// EE Test5
	//EETest5();
	return 0;
}

void EECheck()
{
	// 74 36 58 40
	float3 e0Pos = { 0.566259,0.0,1.68779 };
	float3 e1Pos = { -0.0,0.0,1.08795 };
	float3 e2Pos = { 0.476035,0.1,0.70308 };
	float3 e3Pos = { -0.059103,0.1,1.36439 };

	float3 e0PrdP = { 0.556494, -0.143569, 1.673545 };
	float3 e1PrdP = { -0.006841, 0.182226, 1.077975 };
	float3 e2PrdP = { 0.476035, 0.579123, 0.703080 };
	float3 e3PrdP = { -0.049618, 0.166053, 1.378224 };
	float constrVal = 0.0f;
	Contact contact;
	float restConstr = calcRestEEConstr(e0Pos, e1Pos, e2Pos, e3Pos, contact);
	contact.colliFreeDir = restConstr > 0.0 ? true : false;
	bool shouldResolveEE = contact.colliFreeDir ? 1 : -1;
	printf("1: constrVal is %f\n", restConstr);
	shouldResolveEE = EESingleTest(e0PrdP, e1PrdP, e2PrdP, e3PrdP, &constrVal, contact, 0.05f);
	printf("2: shouldResolveEE is %d; constrVal is %f\n", shouldResolveEE, constrVal);
}

bool EESingleTest(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, float* constrVal, Contact& contact, float thickness)
{
	float3 prdPN = cross(normalize(e1PrdP - e0PrdP), normalize(e3PrdP - e2PrdP));
	if (norm2(prdPN) <= 1e-6)
		return false;
	prdPN = normalize(prdPN);
	contact.n = prdPN;  // update normal of the contact

	double a0 = stp(e3PrdP - e1PrdP, e2PrdP - e1PrdP, prdPN), a1 = stp(e2PrdP - e0PrdP, e3PrdP - e0PrdP, prdPN),
		b0 = stp(e0PrdP - e3PrdP, e1PrdP - e3PrdP, prdPN), b1 = stp(e1PrdP - e2PrdP, e0PrdP - e2PrdP, prdPN);
	contact.w[0] = a0 / (a0 + a1);
	contact.w[1] = a1 / (a0 + a1);
	contact.w[2] = -b0 / (b0 + b1);
	contact.w[3] = -b1 / (b0 + b1);

	// TODO: clamp the weights
	// Recalculate current constr value
	//*constrVal = calcEEConstr(e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness);
	//printf("function constr %f\n", *constrVal);
	printf("colliFreeDir is %d\n", contact.colliFreeDir);
	*constrVal = dot(prdPN * (contact.colliFreeDir ? 1 : -1), (contact.w[0] * e0PrdP + contact.w[1] * e1PrdP) - ((-contact.w[2]) * e2PrdP + (-contact.w[3]) * e3PrdP)) - 2.0 * thickness;
	printf("constrVal is %f\n", *constrVal);
	if ((*constrVal) >= 1e-6)
	{
		return false;
	}
	else
	{
		printf("\tEEResolveTest\n");
		printf("\tprdPN: (%f, %f, %f)\n", prdPN.x, prdPN.y, prdPN.z);
		printf("\te0PrdP: (%f, %f, %f); e1PrdP: (%f, %f, %f); e2PrdP: (%f, %f, %f); e3PrdP: (%f, %f, %f)\n",
			e0PrdP.x, e0PrdP.y, e0PrdP.z,
			e1PrdP.x, e1PrdP.y, e1PrdP.z,
			e2PrdP.x, e2PrdP.y, e2PrdP.z,
			e3PrdP.x, e3PrdP.y, e3PrdP.z);
		printf("\ta0: %f; a1: %f; b0: %f; b1: %f\n", a0, a1, b0, b1);
		printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
		
		return true;
	}
}

// continue simulation
void ContinueSim()
{
	Timer timers[nModules];

	HardwareType ht = CPU;
	PBDObject pbdObj;
<<<<<<< Updated upstream
	pbdObj.ContinueSimInit("D://0319CCDTest//continueSimData//meshTopol//NewLargeClothMeshTopol.19.cache",
		"D://0319CCDTest//continueSimData//constraint//NewLargeClothConstraint.19.cache", ht);
=======
	pbdObj.ContinueSimInit("D://0319CCDTest//VFEE//meshPath//VFEEClothWithSphereMeshTopol.14.cache",
		"D://0319CCDTest//VFEE//constrPath//VFEEClothWithSphereConstr.14.cache", ht);
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
	int substep = 1;
	int iteration = 705;
=======
	string meshPath = "D://0319CCDTest//VFEE//meshPath//VFEEClothWithSphereMeshTopol.";
	string constrPath = "D://0319CCDTest//VFEE//constrPath//VFEEClothWithSphereConstr.";
	string collisionPath = "D://0319CCDTest//VFEE//collisionPath//VFEEClothWithSphereCollision.";

	int substep = gSubStep;
>>>>>>> Stashed changes
	int fps = 24;
	float dt = 1.0 / fps / (float)substep;
	int colPasses = 5;
	int contiCookTimes = 0;
<<<<<<< Updated upstream
	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
		timers[globalTimer].Tick();

		//BufferDebugModule::GetInstance()->CompareStack(contiCookTimes, "ENTRY_PRDB_¼ì²éÊäÈë", pbdObj.constrPBDBuffer.prdPBuffer,0.0f);

		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(dt);
			pbdSolver.ProjectConstraint(st, iteration);
			//BufferDebugModule::GetInstance()->CompareStack(contiCookTimes, "ËãÍê_ProjectConstraint", pbdObj.constrPBDBuffer.prdPBuffer, 0.0f);
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
=======
	for (size_t i = gStartFrame; i <= 5; i++)
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
			pbdSolver.Integration(dt);
			IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//VFEE//continueSimData//results//VFEEClothWithSphereResults." + to_string((i - 1) * substep + s + 1) + ".cache");
			cout << "topol saved" << endl;
			//fixedBuffer.SetName("P");
			//IO::SaveBuffer(fixedBuffer, "D://0319CCDTest//VFEE//continueSimData//12TestCLargeClothWithSphereFixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			//cout << "fixed buffer saved" << endl;
			//vFixedBuffer.SetName("P");
			//IO::SaveBuffer(vFixedBuffer, "D://0319CCDTest//VFEE//continueSimData//12TestCLargeClothWithSphereVFixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			//cout << "v fixed buffer saved" << endl;
			//fFixedBuffer.SetName("P");
			//IO::SaveBuffer(fFixedBuffer, "D://0319CCDTest//VFEE//continueSimData//fixedBuffer//VFEEfixedBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			//cout << "f fixed buffer saved" << endl;
			//pbdObj.velBuffer.SetName("P");
			//IO::SaveBuffer(pbdObj.velBuffer, "D://0319CCDTest//VFEE//continueSimData//velBuffer//VFEEvelBuffer." + to_string((i - 1) * substep + s + 1) + ".cache");
			//cout << "vel buffer saved" << endl;
		}

>>>>>>> Stashed changes
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}

<<<<<<< Updated upstream
=======
// just for collision test
void CollisionTest(int colliPasses) // for later data oriented continue sim
{
	float thickness = 0.03f;
	int debugFrameId = 1;
	int colliResolveIterations = 2;
	Topology meshTopol;
	BufferVector3f prdPBuffer;
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
>>>>>>> Stashed changes

// colli as post processing
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

<<<<<<< Updated upstream
	int substep = gSubStep;
	int cookTimes = 0;
	int contiCookTimes = 0;
	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
		//if (i >= 20)
		//{
		   // BufferDebugModule::GetInstance()->PushStack(contiCookTimes, "ENTRY_PRDB_¼ì²éÊäÈë", pbdObj.constrPBDBuffer.prdPBuffer);
		//}
=======
	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
>>>>>>> Stashed changes
		timers[globalTimer].Tick();  // global timer
		for (size_t s = 0; s < substep; s++)
		{
			pbdSolver.Advect(gDeltaTime);
			pbdSolver.ProjectConstraint(st, gPBDIteration);
<<<<<<< Updated upstream
			//if (i >= 20)
			//{
			   // BufferDebugModule::GetInstance()->PushStack( contiCookTimes, "ËãÍê_ProjectConstraint", pbdObj.constrPBDBuffer.prdPBuffer);
			//}
			for (int col = 0; col < gCollisionPasses; col++)
=======

			//POST processing
			if (gEnableSelf)
>>>>>>> Stashed changes
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
			string path = to_string((i - 1) * substep + s + 1) + ".cache";
			pbdObj.SaveMeshTopol(meshPath + path);
			pbdObj.SaveConstraint(constrPath + path);
			colliSolver.SaveCollision(collisionPath + path);
			//IO::SaveToplogy(pbdObj.meshTopol, "D://0310ContinuousCollisionDectection//ContinuousSimData//5ccdTestLow." + to_string((i-1)*10 + s +1)  + ".cache");
			//cout << "topol saved" << endl;
		}
<<<<<<< Updated upstream
		IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//1largeClothOutput//NewLargeClothWithSphere." + to_string(i) + ".cache");
		printf("---------------------------frame %d topol saved--------------------\n", i);
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());


	}
}

int main()
{
	BufferDebugModule::GetInstance()->Load();
	//RegularSim();
	ContinueSim();

	return 0;
=======
		IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//1largeClothOutput//5NewLargeClothWithSphere." + to_string(i) + ".cache");
		printf("---------------------------frame %d topol saved--------------------\n", i);
		timers[globalTimer].Tock();  // global timer
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}

// colli as constraint
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

	string meshPath = "D://0319CCDTest//VFEE//meshPath//VFEEClothWithSphereMeshTopol.";
	string constrPath = "D://0319CCDTest//VFEE//constrPath//VFEEClothWithSphereConstr.";
	string collisionPath = "D://0319CCDTest//VFEE//collisionPath//VFEEClothWithSphereCollision.";

	for (size_t i = gStartFrame; i <= gEndFrame; i++)
	{
		timers[globalTimer].Tick();
		// for collision debugging
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
		IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//VFEE//TestResults//VFEEClothWithSphereResult." + to_string(i) + ".cache");
		printf("---------------------------frame %d topol saved--------------------\n", i);
		fixedBuffer.SetName("P");
		IO::SaveBuffer(fixedBuffer, "D://0319CCDTest//VFEE//FixedBuffer//VFEEClothWithSphereFixed." + to_string(i) + ".cache");
		printf("---------------------------frame %d fixedBuffer saved--------------------\n", i);
		timers[globalTimer].Tock();
		PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());
	}
}

// --------------- test data -------------------
void Cloth2StaticTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 10;
	Topology meshTopol;
	BufferVector3f prdPBuffer;
	readMeshFromTxt("D://0319CCDTest//VFEESimple//2clothMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//VFEESimple//2clothPrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(6.0f, 6.0f, 6.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(2, 1, 2);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	meshTopol.indices.SetName("indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");
	IO::SaveToplogy(meshTopol, "D://0319CCDTest//VFEESimple//results//results." + to_string(-1) + ".cache");
	printf("---------------------------rest pos saved--------------------\n");
	
	Topology afterResolveTopol;
	prdPBuffer.SetName("P");
	afterResolveTopol.indices = meshTopol.indices;
	afterResolveTopol.primList = meshTopol.primList;
	afterResolveTopol.posBuffer = prdPBuffer;
	afterResolveTopol.indices.SetName("indices");
	afterResolveTopol.primList.SetName("primList");
	afterResolveTopol.posBuffer.SetName("P");
	IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//VFEESimple//results//results." + to_string(0) + ".cache");
	printf("---------------------------prdp pos saved--------------------\n");

	BufferInt vfIndices;
	BufferInt resolveTimes;
	for (int cp = 1; cp <= 10; ++cp)
	{
		printf("CCD Pass %d\n", cp);
		ContactData contactData;
		CCD_SH_Extended(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result for debugging
		printf("contact size:%d\n", contactData.ctxs.GetSize());
		std::set<int> idList;
		//for (int i = 0; i < contactData.ctxStartNum.GetSize(); ++i)
		//{
		//	int id = contactData.ctxStartNum.m_Data[i].x;
		//	if (contactData.ctxs.m_Data[id].type == Contact::VF)
		//		printf("v: %d, f: %d %d %d\n",
		//			contactData.ctxIndices.m_Data[id],
		//			contactData.ctxIndices.m_Data[id + 1],
		//			contactData.ctxIndices.m_Data[id + 2],
		//			contactData.ctxIndices.m_Data[id + 3]);
		//	// idList.insert(contactData.ctxIndices.m_Data[id]);
		//}
		//std::set<int>::iterator it;
		//printf("id: ");
		//for (it = idList.begin(); it != idList.end(); it++)
		//{
		//	printf("%d -- ", *it);
		//}
		//printf("\n");

		CollisionResolveNew(meshTopol, prdPBuffer, contactData, colliResolveIterations, thickness, debugFrameId, vfIndices);

		ContactData contactData1;
		CCD_SH_Narrow(contactData1, shs, meshTopol, prdPBuffer, thickness);
		printf("contact size:%d\n", contactData1.ctxs.GetSize());
		//for (int i = 0; i < contactData1.ctxStartNum.GetSize(); ++i)
		//{
		//	int id = contactData1.ctxStartNum.m_Data[i].x;
		//	if(contactData1.ctxs.m_Data[id].type == Contact::VF)
		//	printf("e0: %d-%d; e1: %d-%d\n",
		//		contactData1.ctxIndices.m_Data[id],
		//		contactData1.ctxIndices.m_Data[id + 1],
		//		contactData1.ctxIndices.m_Data[id + 2],
		//		contactData1.ctxIndices.m_Data[id + 3]);
		//	//idList.insert(contactData.ctxIndices.m_Data[id]);
		//}

		// save result
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//VFEESimple//results//results." + to_string(cp) + ".cache");
		printf("---------------------------resolved topol saved--------------------\n");
	}
}

void Cloth2ChildStaticTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 5;
	Topology meshTopol;
	BufferVector3f prdPBuffer;
	readMeshFromTxt("D://0319CCDTest//VFEEChild//2childClothMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//VFEEChild//2childClothPrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(3.0f, 3.0f, 3.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(3, 2, 3);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	meshTopol.indices.SetName("indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");
	IO::SaveToplogy(meshTopol, "D://0319CCDTest//VFEEChild//results//results." + to_string(-1) + ".cache");
	printf("---------------------------rest pos saved--------------------\n");

	Topology afterResolveTopol;
	prdPBuffer.SetName("P");
	afterResolveTopol.indices = meshTopol.indices;
	afterResolveTopol.primList = meshTopol.primList;
	afterResolveTopol.posBuffer = prdPBuffer;
	afterResolveTopol.indices.SetName("indices");
	afterResolveTopol.primList.SetName("primList");
	afterResolveTopol.posBuffer.SetName("P");
	IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//VFEEChild//results//results." + to_string(0) + ".cache");
	printf("---------------------------prdp pos saved--------------------\n");

	BufferInt vfIndices;
	BufferInt resolveTimes;
	for (int cp = 1; cp <= 10; ++cp)
	{
		printf("CCD Pass %d\n", cp);
		ContactData contactData;
		CCD_SH_Extended(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result for debugging
		printf("contact size:%d\n", contactData.ctxs.GetSize());
		std::set<int> idList;
		//for (int i = 0; i < contactData.ctxStartNum.GetSize(); ++i)
		//{
		//	int id = contactData.ctxStartNum.m_Data[i].x;
		//	if (contactData.ctxs.m_Data[id].type == Contact::VF)
		//		printf("v: %d, f: %d %d %d\n",
		//			contactData.ctxIndices.m_Data[id],
		//			contactData.ctxIndices.m_Data[id + 1],
		//			contactData.ctxIndices.m_Data[id + 2],
		//			contactData.ctxIndices.m_Data[id + 3]);
		//	// idList.insert(contactData.ctxIndices.m_Data[id]);
		//}
		//std::set<int>::iterator it;
		//printf("id: ");
		//for (it = idList.begin(); it != idList.end(); it++)
		//{
		//	printf("%d -- ", *it);
		//}
		//printf("\n");

		CollisionResolveNew(meshTopol, prdPBuffer, contactData, colliResolveIterations, thickness, debugFrameId, vfIndices);

		ContactData contactData1;
		CCD_SH_Narrow(contactData1, shs, meshTopol, prdPBuffer, thickness);
		printf("contact size:%d\n", contactData1.ctxs.GetSize());
		//for (int i = 0; i < contactData1.ctxStartNum.GetSize(); ++i)
		//{
		//	int id = contactData1.ctxStartNum.m_Data[i].x;
		//	if(contactData1.ctxs.m_Data[id].type == Contact::VF)
		//	printf("e0: %d-%d; e1: %d-%d\n",
		//		contactData1.ctxIndices.m_Data[id],
		//		contactData1.ctxIndices.m_Data[id + 1],
		//		contactData1.ctxIndices.m_Data[id + 2],
		//		contactData1.ctxIndices.m_Data[id + 3]);
		//	//idList.insert(contactData.ctxIndices.m_Data[id]);
		//}

		// save result
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//VFEEChild//results//results." + to_string(cp) + ".cache");
		printf("---------------------------resolved topol saved--------------------\n");
	}
}


void Cloth2StaticLargeTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 10;
	Topology meshTopol;
	BufferVector3f prdPBuffer;
	readMeshFromTxt("D://VersionTestData//2clolthLarger//2clothLargeMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://VersionTestData//2clolthLarger//2clothLargePrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(6.0f, 6.0f, 6.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(5, 1, 5);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	BufferInt vfIndices; // for collision debugging
	for (int cp = 1; cp <= 20; ++cp)
	{
		ContactData contactData;
		CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result for debugging
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

		// save result
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

void Cloth3StaticLargeTest()
{
	float thickness = 0.05f;
	int debugFrameId = 1;
	int colliResolveIterations = 10;
	Topology meshTopol;
	BufferVector3f prdPBuffer;
	readMeshFromTxt("D://VersionTestData//2clolthLarger//2clothLargeMeshTopol.txt", meshTopol);
	readBufferFromTxt("D://VersionTestData//2clolthLarger//2clothLargePrdp.txt", prdPBuffer);

	float3 cellSize = make_float3(6.0f, 6.0f, 6.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(5, 1, 5);
	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.InitSH();

	BufferInt vfIndices; // for collision debugging
	for (int cp = 1; cp <= 20; ++cp)
	{
		ContactData contactData;
		CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
		// ccd test result for debugging
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

		// save result
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


// EE Tests

// two simple triangles
void EETest1()
{
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	float stiffnessSetting[1] = { 1.0f };
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;

	Topology meshTopol; // ��һ֡collision free������meshtopol
	BufferVector3f prdPBuffer; // ��ǰ֡������posBuffer
	readMeshFromTxt("D://0319CCDTest//EETests//EETest1//restPos.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//EETests//EETest1//prdP.txt", prdPBuffer);

	//PBDObject pbdObj(dampingRate, gravity, ht);
	//pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	//pbdObj.SetTimer(&timers[pbdObjTimer]);
	//pbdObj.Init(topolFileName, distConstrFileName);

	//SolverPBD pbdSolver;
	//pbdSolver.SetTarget(&pbdObj);
	//pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	float3 cellSize = make_float3(15.0f, 15.0f, 15.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(1, 1, 1);

	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	//CollisionSolver colliSolver;
	//colliSolver.SetThickness(0.03f);
	//colliSolver.SetIterations(gCollisionResolves);
	//colliSolver.SetAcceStruct(&shs);
	//colliSolver.SetTimer(&timers[colliSolverTimer]);

	//colliSolver.CCD_SH(meshTopol.indices, meshTopol.posBuffer, meshTopol.primList, prdPBuffer);
	//colliSolver.CollisionResolve(meshTopol.indices, meshTopol.primList, meshTopol.posBuffer, prdPBuffer, colliSolver.contactData);
	float thickness = 0.03f;
	int iterations = 1;
	ContactData contactData;
	BufferInt ccdIndices;
	int debugFrameId = 0;

	timers[globalTimer].Tick();  // global timer
	CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
	CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);

	//IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//EETests//EETest1//EETest1." + to_string(1) + ".cache");
	//printf("---------------------------frame %d topol saved--------------------\n", 1);
	timers[globalTimer].Tock();  // global timer
	PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());

}

// two mirror triangles
void EETest2()
{
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	float stiffnessSetting[1] = { 1.0f };
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;

	Topology meshTopol; // ��һ֡collision free������meshtopol
	BufferVector3f prdPBuffer; // ��ǰ֡������posBuffer
	readMeshFromTxt("D://0319CCDTest//EETests//EETest2//restPos.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//EETests//EETest2//prdP.txt", prdPBuffer);

	//PBDObject pbdObj(dampingRate, gravity, ht);
	//pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	//pbdObj.SetTimer(&timers[pbdObjTimer]);
	//pbdObj.Init(topolFileName, distConstrFileName);

	//SolverPBD pbdSolver;
	//pbdSolver.SetTarget(&pbdObj);
	//pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	float3 cellSize = make_float3(15.0f, 15.0f, 15.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(1, 1, 1);

	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	//CollisionSolver colliSolver;
	//colliSolver.SetThickness(0.03f);
	//colliSolver.SetIterations(gCollisionResolves);
	//colliSolver.SetAcceStruct(&shs);
	//colliSolver.SetTimer(&timers[colliSolverTimer]);

	//colliSolver.CCD_SH(meshTopol.indices, meshTopol.posBuffer, meshTopol.primList, prdPBuffer);
	//colliSolver.CollisionResolve(meshTopol.indices, meshTopol.primList, meshTopol.posBuffer, prdPBuffer, colliSolver.contactData);

	float thickness = 0.03f;
	int iterations = 1;
	ContactData contactData;
	BufferInt ccdIndices;
	int debugFrameId = 0;

	timers[globalTimer].Tick();  // global timer
	CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
	CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);

	//IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//EETests//EETest1//EETest1." + to_string(1) + ".cache");
	//printf("---------------------------frame %d topol saved--------------------\n", 1);
	timers[globalTimer].Tock();  // global timer
	PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());

}

// many triangles
void EETest3()
{
	Timer timers[nModules];
	float dampingRate = 0.9f;
	float3 gravity = make_float3(0.0, -10.0, 0.0);
	float stiffnessSetting[1] = { 1.0f };
	HardwareType ht = CPU;
	SolverType st = GAUSSSEIDEL;

	Topology meshTopol; // ��һ֡collision free������meshtopol
	BufferVector3f prdPBuffer; // ��ǰ֡������posBuffer
	readMeshFromTxt("D://0319CCDTest//EETests//EETest3//restPos.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//EETests//EETest3//prdP.txt", prdPBuffer);

	//PBDObject pbdObj(dampingRate, gravity, ht);
	//pbdObj.SetConstrOption(DISTANCE, stiffnessSetting);
	//pbdObj.SetTimer(&timers[pbdObjTimer]);
	//pbdObj.Init(topolFileName, distConstrFileName);

	//SolverPBD pbdSolver;
	//pbdSolver.SetTarget(&pbdObj);
	//pbdSolver.SetTimer(&timers[pbdSolverTimer]);

	float3 cellSize = make_float3(15.0f, 15.0f, 15.0f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(1, 1, 1);

	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();


	//CollisionSolver colliSolver;
	//colliSolver.SetThickness(0.03f);
	//colliSolver.SetIterations(gCollisionResolves);
	//colliSolver.SetAcceStruct(&shs);
	//colliSolver.SetTimer(&timers[colliSolverTimer]);

	//colliSolver.CCD_SH(meshTopol.indices, meshTopol.posBuffer, meshTopol.primList, prdPBuffer);
	//colliSolver.CollisionResolve(meshTopol.indices, meshTopol.primList, meshTopol.posBuffer, prdPBuffer, colliSolver.contactData);

	float thickness = 0.03f;
	int iterations = 1;
	ContactData contactData;
	BufferInt ccdIndices;
	int debugFrameId = 0;

	timers[globalTimer].Tick();  // global timer
	// CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
	// CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
	int i;
	for (i = 0; i < 500; ++i)
	{
		printf("Collision Pass %d\n", i);
		if (i % 10 == 0)
		{
			CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
			if (contactData.ctxs.GetSize() == 0)
				break;
		}

		CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
	}
	printf("num of iterations: %d\n", i);
	//colliSolver.CCD_SH(meshTopol.indices, meshTopol.posBuffer, meshTopol.primList, prdPBuffer);
	//colliSolver.CollisionResolve(meshTopol.indices, meshTopol.primList, meshTopol.posBuffer, prdPBuffer, colliSolver.contactData);

	//IO::SaveToplogy(pbdObj.meshTopol, "D://0319CCDTest//EETests//EETest1//EETest1." + to_string(1) + ".cache");
	//printf("---------------------------frame %d topol saved--------------------\n", 1);
	timers[globalTimer].Tock();  // global timer
	PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());

}

void EETest4()  // two small clothes
{

	Timer timers[nModules];

	Topology meshTopol; // ��һ֡collision free������meshtopol
	BufferVector3f prdPBuffer; // ��ǰ֡������posBuffer
	readMeshFromTxt("D://0319CCDTest//EETests//EETest4//restPos.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//EETests//EETest4//prdP.txt", prdPBuffer);

	float3 cellSize = make_float3(2.0f, 2.0f, 2.0f);
	float3 gridCenter = make_float3(0.0f, 0.5f, 0.0f);
	uint3 gridSize = make_uint3(3, 2, 3);

	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	meshTopol.indices.SetName("indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");
	IO::SaveToplogy(meshTopol, "D://0319CCDTest//EETests//EETest4//results//results." + to_string(-1) + ".cache");
	printf("---------------------------rest pos saved--------------------\n");

	Topology beforeResolveTopol;
	prdPBuffer.SetName("P");
	beforeResolveTopol.indices = meshTopol.indices;
	beforeResolveTopol.primList = meshTopol.primList;
	beforeResolveTopol.posBuffer = prdPBuffer;
	beforeResolveTopol.indices.SetName("indices");
	beforeResolveTopol.primList.SetName("primList");
	beforeResolveTopol.posBuffer.SetName("P");
	IO::SaveToplogy(beforeResolveTopol, "D://0319CCDTest//EETests//EETest4//results//results." + to_string(0) + ".cache");
	printf("---------------------------prdp pos saved--------------------\n");

	int debugFrameId = 0;
	float thickness = 0.05f;
	int iterations = 1;
	ContactData contactData;
	BufferInt ccdIndices;

	timers[globalTimer].Tick();  // global timer
	// CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
	// CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
	int i;
	for (i = 1; i <= 1; ++i)
	{
		// printf("Collision Pass %d\n", i);

		CCD_SH_Extended(contactData, shs, meshTopol, prdPBuffer, thickness);
		if (contactData.ctxs.GetSize() == 0)
			break;

		CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//EETests//EETest4//results//results." + to_string(i) + ".cache");
		printf("---------------------------resolved topol saved--------------------\n");
	}
	//printf("num of iterations: %d\n", i);
	
	timers[globalTimer].Tock();  // global timer
	PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());

}


void EETest5()  // two small clothes
{

	Timer timers[nModules];

	Topology meshTopol; // ��һ֡collision free������meshtopol
	BufferVector3f prdPBuffer; // ��ǰ֡������posBuffer
	readMeshFromTxt("D://0319CCDTest//EETests//EETest5//restPos.txt", meshTopol);
	readBufferFromTxt("D://0319CCDTest//EETests//EETest5//prdP.txt", prdPBuffer);

	float3 cellSize = make_float3(2.0f, 2.0f, 2.0f);
	float3 gridCenter = make_float3(0.0f, 0.5f, 0.0f);
	uint3 gridSize = make_uint3(3, 2, 3);

	SpatialHashSystem shs(prdPBuffer, meshTopol.indices, CPU, gridCenter, gridSize, cellSize);
	shs.SetTimer(&timers[shsTimer]);
	shs.InitSH();

	meshTopol.indices.SetName("indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");
	IO::SaveToplogy(meshTopol, "D://0319CCDTest//EETests//EETest5//results//results." + to_string(-1) + ".cache");
	printf("---------------------------rest pos saved--------------------\n");

	Topology beforeResolveTopol;
	prdPBuffer.SetName("P");
	beforeResolveTopol.indices = meshTopol.indices;
	beforeResolveTopol.primList = meshTopol.primList;
	beforeResolveTopol.posBuffer = prdPBuffer;
	beforeResolveTopol.indices.SetName("indices");
	beforeResolveTopol.primList.SetName("primList");
	beforeResolveTopol.posBuffer.SetName("P");
	IO::SaveToplogy(beforeResolveTopol, "D://0319CCDTest//EETests//EETest4//results//results." + to_string(0) + ".cache");
	printf("---------------------------prdp pos saved--------------------\n");

	int debugFrameId = 0;
	float thickness = 0.05f;
	int iterations = 5;
	ContactData contactData;
	BufferInt ccdIndices;

	timers[globalTimer].Tick();  // global timer
	// CCD_SH(contactData, shs, meshTopol, prdPBuffer, thickness);
	// CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
	int i;
	for (i = 1; i <= 10; ++i)
	{
		// printf("Collision Pass %d\n", i);

		CCD_SH_Extended(contactData, shs, meshTopol, prdPBuffer, thickness);
		if (contactData.ctxs.GetSize() == 0)
			break;

		CollisionResolveNew(meshTopol, prdPBuffer, contactData, iterations, thickness, debugFrameId, ccdIndices);
		Topology afterResolveTopol;
		prdPBuffer.SetName("P");
		afterResolveTopol.indices = meshTopol.indices;
		afterResolveTopol.primList = meshTopol.primList;
		afterResolveTopol.posBuffer = prdPBuffer;
		afterResolveTopol.indices.SetName("indices");
		afterResolveTopol.primList.SetName("primList");
		afterResolveTopol.posBuffer.SetName("P");
		IO::SaveToplogy(afterResolveTopol, "D://0319CCDTest//EETests//EETest5//results//results." + to_string(i) + ".cache");
		printf("---------------------------resolved topol saved--------------------\n");
	}
	//printf("num of iterations: %d\n", i);

	timers[globalTimer].Tock();  // global timer
	PBD_DEBUGTIME(timers[globalTimer].GetFuncTime());

>>>>>>> Stashed changes
}