#include "PBD_Basic.cuh"


float Distance(float3 p1, float3 p2)
{
	return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
}
// ----------------------------------------------------------------------------------------

void PBDObject::groundTruthTest()
{
	vector<int> arr0 = { 0, 1, 3, 3, 2, 0 };
	vector<int> arr1 = {0,3,0,1,0,2,1,3,2,3};
	constrPBDBuffer.topol.indices.m_Data = arr1;
	vector<int2> arr2 = { make_int2(0,2),  make_int2(2,2),  make_int2(4,2),  make_int2(6,2),  make_int2(8,2) };
	constrPBDBuffer.topol.primList.m_Data = arr2;
	vector<int> arr3 = { -1,-1,-1,-1,-1 };
	constrPBDBuffer.color.m_Data = arr3;
	constrPBDBuffer.prdColor.m_Data = arr3;
}

// ----------------------------------------------------------------------------------------



void PBDObject::Init()
{
	// OpenGL Topology
	meshTopol.indices.SetName("Indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");

	CreatePosition(meshTopol.posBuffer, make_float2(0.0, 0.0), sizeX, sizeY, resY, resX);
	//CreateOpenGLIndices(meshTopol.indices, resY, resX);
	float stiffnessBuff[1] = { 1.0f };
	groundTruthTest();
	InitConstr(1, 1.0f, stiffnessBuff);
	
}

void PBDObject::InitConstr(int constrNumSetting, float unitMass, float* stiffnesses)
{
	constrPBDBuffer.topol.posBuffer = meshTopol.posBuffer;
	constrPBDBuffer.prdPBuffer = meshTopol.posBuffer;
	constrPBDBuffer.velBuffer.m_Data.resize(constrPBDBuffer.prdPBuffer.GetSize(), make_float3(0.0f, 0.0f, 0.0f));

	constrPBDBuffer.topol.primList.SetName("primList");
	constrPBDBuffer.topol.indices.SetName("Indices");
	constrPBDBuffer.color.SetName("color");
	constrPBDBuffer.prdColor.SetName("color");
	/*for (int i = 0; i < 3; i++)
	{
		printf("postion Buffer: %f, %f, %f \n", constrPBDBuffer.topol.posBuffer.m_Data[i].x, constrPBDBuffer.topol.posBuffer.m_Data[i].y, constrPBDBuffer.topol.posBuffer.m_Data[i].z);
	}*/

	std::set<unsigned long>::iterator it;
	// assert(numOfConstr == stiffnesses.size());
	for (int i = 0; i < constrNumSetting; ++i)
	{
		switch (i)
		{
		case DISTANCE:
			//CreateDistanceIndices(constrPBDBuffer.topol.indices, resY, resX);
			//CreateSingleDistConstr(constrPBDBuffer.topol.posBuffer, constrPBDBuffer.topol.primList, constrPBDBuffer.topol.indices, stiffnesses[DISTANCE], unitMass);
			GenePoint2PrimsMap(constrPBDBuffer.topol);
			GenePrim2PrimsMap(constrPBDBuffer.topol);
			EdgeColoring(1000);
			/*printf("Indices: \n");
			printf("\t");
			for (int i = 0; i < constrPBDBuffer.topol.indices.GetSize(); ++i)
			{
				printf("%d-", constrPBDBuffer.topol.indices.m_Data[i]);
			}
			printf("\nPoint2PrimsMap: \n");
			
			for (int i = 0; i < constrPBDBuffer.Point2PrimsMap.size(); ++i)
			{
				printf("\t");
				printf("point %d (size: %d): ", i, constrPBDBuffer.Point2PrimsMap[i].size());
				for (int j = 0; j < constrPBDBuffer.Point2PrimsMap[i].size(); ++j)
				{
					printf("%d-", constrPBDBuffer.Point2PrimsMap[i][j]);
				}
				printf("\n");
			}
*/
			printf("\nPrim2PrimsMap: \n");
			printf("size of prim2prims indices%d-", constrPBDBuffer.Prim2PrimsMap.indices.GetSize());
			printf("\n");
			for (int i = 0; i < constrPBDBuffer.Prim2PrimsMap.indices.GetSize(); ++i)
			{
				printf("%d-", constrPBDBuffer.Prim2PrimsMap.indices.m_Data[i]);
			}
			printf("\n");
			for (int i = 0; i < constrPBDBuffer.Prim2PrimsMap.startNumList.GetSize(); ++i)
			{
				printf("startIdx: %d, neighbour Num: %d\n", constrPBDBuffer.Prim2PrimsMap.startNumList.m_Data[i].x, constrPBDBuffer.Prim2PrimsMap.startNumList.m_Data[i].y);
			}

			/*printf("primMap (size: %d): ",constrPBDBuffer.Prim2PrimsMap.size());
			for (int i = 0; i < constrPBDBuffer.Prim2PrimsMap.size(); ++i)
			{
				printf("\t");
				printf("prim %d (size: %d): ", i, constrPBDBuffer.Prim2PrimsMap[i].size());
				for  (auto it = constrPBDBuffer.Prim2PrimsMap[i].begin(); it != constrPBDBuffer.Prim2PrimsMap[i].end(); ++it)
				{
					printf("%d-", *it);
				}
				printf("\n");
			}*/
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
	for (int i = num - resX; i < num - 1; ++i)
	{
		indices.m_Data.push_back(i);
		indices.m_Data.push_back(i + 1);
	}
}

void PBDObject::CreateSingleDistConstr(BufferVector3f& positionBuffer, BufferInt2& primList, BufferInt& indices, float stiffness, float unitMass)
{
	int count = 2;
	for (int i = 0; i < indices.GetSize() / count; i++)
	{
		// init primList
		int2 p;
		p.x = i * count;
		p.y = count;
		primList.m_Data.push_back(p);
		// init stiffness
		constrPBDBuffer.stiffness.m_Data.push_back(stiffness);
		// init rest length
		int i0 = indices.m_Data[p.x];
		int i1 = indices.m_Data[p.x + 1];
		float d = Distance(positionBuffer.m_Data[i0], positionBuffer.m_Data[i1]);
		constrPBDBuffer.restLength.m_Data.push_back(d);
		// init mass
		constrPBDBuffer.mass.m_Data.push_back(unitMass);
		constrPBDBuffer.mass.m_Data.push_back(unitMass);
		// init contraint type
		constrPBDBuffer.constraintType.m_Data.push_back(DISTANCE);
		// init color
		constrPBDBuffer.color.m_Data.push_back(-1);
		constrPBDBuffer.prdColor.m_Data.push_back(-1);
		constrPBDBuffer.sortedColor.m_Data.push_back(-1);
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

void PBDObject::GenePoint2PrimsMap(Topology topol)
{
	printf("entered gen point2prims\n");
	auto primList =&(topol.primList);
	auto indices = &(topol.indices);
	for (int primId = 0; primId < primList->GetSize(); ++primId)
	{
		int2 currPrim = primList->m_Data[primId];
		for (int i = 0; i < currPrim.y; ++i)
		{
			constrPBDBuffer.Point2PrimsMap[indices->m_Data[currPrim.x+i]].push_back(primId);
		}
	}
	printf("end gen point2prims\n");
}


void PBDObject::GenePrim2PrimsMap(Topology topol)
{
	printf("entered gen Prim2PrimsMap\n");
	auto primList = &(topol.primList);
	auto indices = &(topol.indices);
	auto Point2PrimsMap = constrPBDBuffer.Point2PrimsMap;
	auto Prim2PrimsMap = &(constrPBDBuffer.Prim2PrimsMap);
	for (int primId = 0; primId < primList->GetSize(); ++primId)
	{
		//std::set<int> linkedPrimsSet;
		std::vector<int> linkedPrimsSet;
		int nptPrims = primList->m_Data[primId].y;
		for (int ptId = 0; ptId < nptPrims; ++ptId)
		{
			int currPtId = indices->m_Data[primList->m_Data[primId].x + ptId];
			// printf("primId: %d; ", primId);
			auto linkedPrims = Point2PrimsMap[currPtId];
			// printf("linked primtive id: ");
			for (int nlp = 0; nlp < linkedPrims.size(); nlp++)
			{
				int linkPrimId = linkedPrims[nlp];
				
				if (linkPrimId == primId)
					continue;
				linkedPrimsSet.push_back(linkPrimId);
			}			
		}
		int startIdx = Prim2PrimsMap->indices.GetSize();
		Prim2PrimsMap->indices.m_Data.insert(
			std::end(Prim2PrimsMap->indices.m_Data),
			std::begin(linkedPrimsSet),
			std::end(linkedPrimsSet));
		Prim2PrimsMap->startNumList.m_Data.push_back(make_int2(startIdx, linkedPrimsSet.size()));
	}
	printf("end gen Prim2PrimsMap\n");
}

void PBDObject::AssignColorsCPU()
{
	auto p2pIndices = &(constrPBDBuffer.Prim2PrimsMap.indices);
	auto p2pStartNumList = &(constrPBDBuffer.Prim2PrimsMap.startNumList);
	auto color = &(constrPBDBuffer.color);
	auto prdColor = &(constrPBDBuffer.prdColor);
	for (int idx = 0; idx < constrPBDBuffer.topol.primList.GetSize(); idx++)
	{
		if (idx == 0)
			detectConflict = 0;

		if (color->m_Data[idx] >= 0)
			continue;

		int nidx = p2pStartNumList->m_Data[idx].x;
		//int nlast = p2pStartNumList->m_Data[idx].x + p2pStartNumList->m_Data[idx].y;
		int nlast = p2pStartNumList->m_Data[idx].x+ p2pStartNumList->m_Data[idx].y;
		int c = -1, offset = 0;
		while (c < 0)
		{
			// Bit flag to store seen neighbor colors.
			unsigned long forbidden = 0;
			for (int i = nidx; i < nlast; ++i)
			{
				int n = p2pIndices->m_Data[i];
				int nc = color->m_Data[n] - offset;
				if (nc >= 0 && nc < FORBIDBITS)
					forbidden |= (1ul << nc);
			}
			// Check if there's an open color in the current bitrange.
			if (forbidden ^ MAXFORBID)
			{
				unsigned long x = forbidden;
				c = offset;
				// Find position of first zero bit.
					x = (~x) & (x + 1);
				// Color is log2(x)
				while (x > 1)
				{
					x /= 2;
					c++;
				}
			}
			else
			{
				// Otherwise we need to try again with the next range of colors.
				offset += FORBIDBITS;
			}
		}
		// Record speculative coloring.
		prdColor->m_Data[idx] = c;
	}
}

void PBDObject::ResolveConflictsCPU()
{
	auto primList = &(constrPBDBuffer.topol.primList);
	auto p2pIndices = &(constrPBDBuffer.Prim2PrimsMap.indices);
	auto p2pStartNumList = &(constrPBDBuffer.Prim2PrimsMap.startNumList);
	auto color = &(constrPBDBuffer.color);
	auto prdColor = &(constrPBDBuffer.prdColor);

	for (int idx = 0; idx < primList->GetSize(); ++idx)
	{
		// Nothing to do if already colored.
		if (color->m_Data[idx] >= 0)
			continue;

		int nidx = p2pStartNumList->m_Data[idx].x;
		int nlast = p2pStartNumList->m_Data[idx].x + p2pStartNumList->m_Data[idx].y;
		int c = prdColor->m_Data[idx];
		//int c = newcolor[idx];
		int conflict = 0;
		int npt = primList->m_Data[idx].y;

		for (int i = nidx; i < nlast; ++i)
		{
			int n = p2pIndices->m_Data[i];
			int pc = prdColor->m_Data[n];

			// Check for conflict.
			if (pc == c)
			{
				int nnpt = primList->m_Data[n].y;
				// Resolution gives preference to primitives with more points, 
				// otherwise the one that comes first.
				// (NOTE: We color in fewer iterations if we prioritize by number
				// of graph neighbors, but then assuming we process prims by color,
				// we're usually farther away from original point order leading to many
				// cache misses and slower downstream processing.)
				if (nnpt > npt ||
					(nnpt == npt && n < idx))
				{
					conflict = 1;
					break;
				}
			}
		}

		// If there's a conflict then reset sizes for more work,
		// otherewise accept speculative color.
		if (conflict)
		{
			detectConflict = primList->GetSize();
			break;
		}
		else
			color->m_Data[idx] = c;

	}
}

void PBDObject::EdgeColoring(int iterations)
{
	for (int i = 0; i < iterations; ++i)
	{
		printf("iteration %d color: ", i);
		for (int i = 0; i < constrPBDBuffer.color.GetSize(); ++i)
		{
			printf("%d // ", constrPBDBuffer.color.m_Data[i]);
		}
		printf("\n");
		printf("iteration %d prd color: ", i);
		for (int i = 0; i < constrPBDBuffer.color.GetSize(); ++i)
		{
			printf("%d // ", constrPBDBuffer.color.m_Data[i]);
		}
		printf("\n");
		AssignColorsCPU();
		ResolveConflictsCPU();
		printf("iteration %d color: ", i);
		for (int i = 0; i < constrPBDBuffer.color.GetSize(); ++i)
		{
			printf("%d // ", constrPBDBuffer.color.m_Data[i]);
		}
		printf("\n");
		printf("iteration %d prd color: ", i);
		for (int i = 0; i < constrPBDBuffer.color.GetSize(); ++i)
		{
			printf("%d // ", constrPBDBuffer.color.m_Data[i]);
		}
		printf("\n");
		if (detectConflict == 0)
			break;
		string fileName = "D:/colorEdge/color." + to_string(i) + ".cache";
		IO::SaveBuffer(constrPBDBuffer.prdColor, fileName);
		cout << "color saved" << endl;	
	}
}

//
//  iteartion
//	
//


void SolverPBD::Advect(float dt, HardwareType ht)
{
	auto velBuffer = &(pbdObj->constrPBDBuffer.velBuffer);
	auto prdPBuffer = &(pbdObj->constrPBDBuffer.prdPBuffer);
	auto positionBuffer = &(pbdObj->meshTopol.posBuffer);
	switch (ht)
	{
	case CPU:
	
		for (int i = 0; i < velBuffer->GetSize(); i++)
		{
			/*if (i == 30)
				printf("old velocity Buffer: %f, %f, %f \n", velBuffer->m_Data[i].x, velBuffer->m_Data[i].y, velBuffer->m_Data[i].z);*/
			velBuffer->m_Data[i] += pbdObj->gravity * dt;
			/*if(i == 30)
				printf("new velocity Buffer: %f, %f, %f \n", velBuffer->m_Data[i].x, velBuffer->m_Data[i].y, velBuffer->m_Data[i].z);*/
			velBuffer->m_Data[i] *= powf(pbdObj->dampingRate, dt);
			
		}

		for (int j = 0; j < prdPBuffer->GetSize(); j++)
		{
			prdPBuffer->m_Data[j] = positionBuffer->m_Data[j] + velBuffer->m_Data[j] * dt;

			//printf("postion Buffer: %f, %f, %f \n", prdPBuffer.m_Data[j].x, prdPBuffer.m_Data[j].y, prdPBuffer.m_Data[j].z);

		}
		break;
	case GPU:
		break;
	default:
		break;
	}
}

void SolverPBD::ProjectConstraint(HardwareType ht, SolverType st, int iterations)
{
	auto primList = &(pbdObj->constrPBDBuffer.topol.primList);
	auto prdPBuffer = &(pbdObj->constrPBDBuffer.prdPBuffer);
	auto massBuffer = &(pbdObj->constrPBDBuffer.mass);
	auto restLengthBuffer = &(pbdObj->constrPBDBuffer.restLength);
	auto stiffnessBuffer = &(pbdObj->constrPBDBuffer.stiffness);
	auto restPosBuffer = &(pbdObj->restPosBuffer);
	auto indices = &(pbdObj->constrPBDBuffer.topol.indices);
	switch (ht)
	{
	case CPU:
		for (size_t ii = 0; ii < iterations; ii++)
		{
			for (size_t i = 0; i < primList->GetSize(); i++)
			{
				if (primList->m_Data[i].y != 2)
					continue;
				int i0 = indices->m_Data[primList->m_Data[i].x];
				int i1 = indices->m_Data[primList->m_Data[i].x + 1];
				float3 dp1;
				float3 dp2;
				float d = Distance(prdPBuffer->m_Data[i0], prdPBuffer->m_Data[i1]);
				float3 v;
				v = prdPBuffer->m_Data[i0] - prdPBuffer->m_Data[i1];
				dp1.x = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.x / d;
				dp1.y = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.y / d;
				dp1.z = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.z / d;
				dp2.x = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.x / d;
				dp2.y = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.y / d;
				dp2.z = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * v.z / d;
				float k = 1 - powf(1 - stiffnessBuffer->m_Data[i], 1.0 / (ii + 1));
				dp1 *= k;
				dp2 *= k;
				prdPBuffer->m_Data[i0] += dp1;
				prdPBuffer->m_Data[i1] += dp2;
			}

			for (size_t j = 0; j < prdPBuffer->GetSize(); j++)
			{
				//attach points
				if (j == 0)
				{
					prdPBuffer->m_Data[j] = restPosBuffer->m_Data[0];
				}
				if ( j == pbdObj->resX - 1)
				{
					prdPBuffer->m_Data[j] = restPosBuffer->m_Data[1];
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
	auto positionBuffer = &(pbdObj->meshTopol.posBuffer);
	auto velBuffer = &(pbdObj->constrPBDBuffer.velBuffer);
	auto prdPBuffer= &(pbdObj->constrPBDBuffer.prdPBuffer);
	switch (ht)
	{
	case CPU:
		for (size_t i = 0; i < positionBuffer->GetSize(); i++)
		{
			velBuffer->m_Data[i] = (prdPBuffer->m_Data[i] - positionBuffer->m_Data[i]) / dt;
			positionBuffer->m_Data[i] = prdPBuffer->m_Data[i];
		}
		break;
	case GPU:
		break;
	default:
		break;
	}
}
