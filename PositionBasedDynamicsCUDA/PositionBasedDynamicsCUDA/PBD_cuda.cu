#include<vector>
#include<assert.h>
#include<vector_types.h>
#include<vector_functions.h>
#include "cuda_runtime.h"
#include "helper_math.h"
#include "device_launch_parameters.h"
#include "PBD_Basic.cuh"


float Distance(float3 p1, float3 p2)
{
	return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
}
// ----------------------------------------------------------------------------------------

void PBDObject::Init()
{
	// OpenGL Topology
	CreatePosition(meshTopol.posBuffer, make_float2(0.0, 0.0), sizeX, sizeY, resY, resX);
	CreateOpenGLIndices(meshTopol.indices, resY, resX);
	float stiffnessBuff[1] = { 1.0f };
	InitConstr(1, 1.0f, stiffnessBuff);
}

void PBDObject::InitConstr(int numOfConstr, float unitMass, float* stiffnesses)
{
	constrPBDBuffer.topol.posBuffer = meshTopol.posBuffer;
	constrPBDBuffer.prdPBuffer = meshTopol.posBuffer;
	memset(&constrPBDBuffer.velBuffer.m_Data, 0, sizeof(constrPBDBuffer.topol.posBuffer.m_Data));
	fill(constrPBDBuffer.velBuffer.m_Data.begin(), constrPBDBuffer.velBuffer.m_Data.end(), 0);
	
	// assert(numOfConstr == stiffnesses.size());
	for (int i = 0; i < numOfConstr; ++i)
	{
		switch (i)
		{
		case DISTANCE:
			CreateDistanceIndices(constrPBDBuffer.topol.indices, resY, resX);
			CreateSingleDistConstr(constrPBDBuffer.topol.posBuffer, constrPBDBuffer.topol.primList, constrPBDBuffer.topol.indices, stiffnesses[DISTANCE], unitMass);
			break;
		case BENDING:
			break;
		case ANCHOR:
			break;
		default:
			break;
		}
	}
}

void PBDObject::CreateDistanceIndices(BufferInt& indices, int resY, int resX)
{
	int num = resY * resX;
	for (int i = 0; i < num - resX; i++)
	{
		if (i % resX == 0)
		{
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + 1);
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + resX);
		}
		else if (i % resX == resX - 1)
		{
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + resX);
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + resX - 1);
		}
		else
		{
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + 1);
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + resX);
			indices.m_Data.push_back(i);
			indices.m_Data.push_back(i + resX - 1);
		}
	}
}

void PBDObject::CreateSingleDistConstr(BufferVector3f& positionBuffer, BufferInt2& primList, BufferInt& indices, float stiffness, float unitMass)
{
	int count = 2;
	for (int i = 0; i < indices.GetSize() / count; i++)
	{
		//init primList
		int2 p;
		p.x = i * count;
		p.y = count;
		primList.m_Data.push_back(p);
		//init stiffness
		constrPBDBuffer.stiffness.m_Data.push_back(stiffness);
		//init rest length
		int i0 = indices.m_Data[p.x];
		int i1 = indices.m_Data[p.x + 1];
		float d = Distance(positionBuffer.m_Data[i0], positionBuffer.m_Data[i1]);
		constrPBDBuffer.restLength.m_Data.push_back(d);
		//init mass
		constrPBDBuffer.mass.m_Data.push_back(unitMass);
		constrPBDBuffer.mass.m_Data.push_back(unitMass);
		constrPBDBuffer.constraintType.m_Data.push_back(DISTANCE);
	}
}


void PBDObject::CreatePosition(BufferVector3f& positionBuffer, float2 cord, float sizeX, float sizeY, int resY, int resX)
{
	float lengthInterval = sizeX / (resX - 1);
	float heightInterval = sizeY / (resY - 1);
	int num = resY * resX;
	int index = 0;
	for (int i = 0; i < resY; i++)
	{
		for (int j = 0; j < resX; j++)
		{
			float3 p;
			p.x = cord.x + j * lengthInterval;
			p.y = 0;
			p.z = cord.y + i * heightInterval;
			positionBuffer.m_Data.push_back(p);
			index++;
		}
	}
	restPosBuffer.m_Data.push_back(positionBuffer.m_Data[0]) ; 
	restPosBuffer.m_Data.push_back(positionBuffer.m_Data[resX-1]) ; 
}

void PBDObject::CreateOpenGLIndices(BufferInt& openGLIndices, int resY, int resX)
{
	int num = resY * resX;
	for (int i = 0; i < num - resX; i++)
	{
		if (i % resX == resX - 1)
			continue;
		openGLIndices.m_Data.push_back(i);
		openGLIndices.m_Data.push_back(i + resX + 1);
		openGLIndices.m_Data.push_back(i + resX);
		openGLIndices.m_Data.push_back(i);
		openGLIndices.m_Data.push_back(i + 1);
		openGLIndices.m_Data.push_back(i + resX + 1);

	}
}

void SolverPBD::Advect(float dt, HardwareType ht)
{
	switch (ht)
	{
	case CPU:
		break;
	case GPU:
		break;
	default:
		break;
	}
}

void SolverPBD::ProjectConstraint(HardwareType ht, SolverType st, int iterations)
{
	auto primList = pbdObj->constrPBDBuffer.topol.primList;
	auto prdPBuffer = pbdObj->constrPBDBuffer.prdPBuffer;
	auto massBuffer = pbdObj->constrPBDBuffer.mass;
	auto restLengthBuffer = pbdObj->constrPBDBuffer.restLength;
	auto stiffnessBuffer = pbdObj->constrPBDBuffer.stiffness;
	auto restPosBuffer = pbdObj->restPosBuffer;
	auto indices = pbdObj->constrPBDBuffer.topol.indices;
	switch (ht)
	{
	case CPU:
		for (size_t ii = 0; ii < iterations; ii++)
		{
			for (size_t i = 0; i < primList.GetSize(); i++)
			{
				if (primList.m_Data[i].y != 2)
					continue;
				int i0 = indices.m_Data[primList.m_Data[i].x];
				int i1 = indices.m_Data[primList.m_Data[i].x + 1];
				float3 dp1;
				float3 dp2;
				float d = Distance(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1]);
				float3 v;
				v = prdPBuffer.m_Data[i0] - prdPBuffer.m_Data[i1];
				dp1.x = -massBuffer.m_Data[i0] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.x / d;
				dp1.y = -massBuffer.m_Data[i0] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.y / d;
				dp1.z = -massBuffer.m_Data[i0] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.z / d;
				dp2.x = massBuffer.m_Data[i1] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.x / d;
				dp2.y = massBuffer.m_Data[i1] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.y / d;
				dp2.z = massBuffer.m_Data[i1] / (massBuffer.m_Data[i0] + massBuffer.m_Data[i1]) * (d - restLengthBuffer.m_Data[i]) * v.z / d;
				float k = 1 - powf(1 - stiffnessBuffer.m_Data[i], 1.0 / (ii + 1));
				dp1 *= k;
				dp2 *= k;
				prdPBuffer.m_Data[i0] += dp1;
				prdPBuffer.m_Data[i1] += dp2;
			}

			for (size_t j = 0; j < prdPBuffer.GetSize(); j++)
			{
				//attach points
				if (j == 0)
				{
					prdPBuffer.m_Data[j] = restPosBuffer.m_Data[0];
				}
				if ( j == pbdObj->resX - 1)
				{
					prdPBuffer.m_Data[j] = restPosBuffer.m_Data[1];
				}

				////point collide with sphere
				//bool isCollideSphere = ColliderSphere(prdPBuffer.m_Data[j], sphereOrigin, sphereRadius, j);
				//if (isCollideSphere) //move the point to the point which intersect with sphere
				//{
				//	float3 moveVector = GenerateMoveVectorSphere(sphereOrigin, sphereRadius, prdPBuffer.m_Data[j], j);
				//	prdPBuffer.m_Data[j] += moveVector;
				//}
				////point collide with ground
				//bool isCollideGoround = CollideGround(prdPBuffer.m_Data[j], groundCenter);
				//if (isCollideGoround)
				//{
				//	prdPBuffer.m_Data[j].y = groundCenter.y;
				//}
			}
		}
		break;
	case GPU:
		for (int i = 0; i < iterations; i++)
		{
			// TODO: edge coloring
			/*
			int* color
			int* sortedColor
			int* sortedPrim

			int main
			{
				sort(color, sortedColor, sortedPrim)
				//color: 11220001110222
				//sortedColor:000111222
				//sortedPrim:
				int2 workSets = Eval(SortedColor)
				//workSets.x start
				//workSets.y num
				for(auto workset: workSets)
				{
					int start  = workSets.x
					int num = workSets.y
					kernel project<<<numBlock, numThread>>>(num, sortedPrim+start, prdPbuffer, restlength)
					{
						if(index > = num)
							return;
						int2 prim = sortedPrimId[index]
						int i0 = primList[prim.x]
						int i1 = primList[prim.x + 1]

						float3 p0 = prdPbuffer[i0]
						float3 p1 = prdPbuffer[i1]

					}
				}
			}

			*/
		}
		break;
	default:
		break;
	}
}

void SolverPBD::Integration(float dt, HardwareType ht)
{
	auto positionBuffer = pbdObj->meshTopol.posBuffer;
	auto velBuffer = pbdObj->constrPBDBuffer.velBuffer;
	auto prdPBuffer= pbdObj->constrPBDBuffer.prdPBuffer;
	printf("vel: %d, pos: %d, prd: %d\n", velBuffer.GetSize(), positionBuffer.GetSize(), prdPBuffer.GetSize());
	switch (ht)
	{
	case CPU:
		for (size_t i = 0; i < positionBuffer.GetSize(); i++)
		{
			//velBuffer.m_Data[i] = (prdPBuffer.m_Data[i] - positionBuffer.m_Data[i]) / dt;
			//positionBuffer.m_Data[i] = prdPBuffer.m_Data[i];
		}
		break;
	case GPU:
		break;
	default:
		break;
	}
}


