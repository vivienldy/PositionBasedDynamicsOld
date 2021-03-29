#include "PBD_Basic.cuh"
#include"CCD_Basic.h"
KERNEL_FUNC float Distance(float3 p1, float3 p2)
{
	return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
}

// ----------------------------------------------------------------------------------------
// ConstraintPBD 
void ConstraintPBD::InitDistanceConstr(BufferVector3f& meshPosBuffer, float stiffness, int resY, int resX)
{
	// cout << __FUNCTION__ << "   resolution:" << resY << "----" << resX << std::endl;

	InitDistanceIndices(resY, resX);

	InitDistanceInfo(meshPosBuffer, stiffness);
}

void ConstraintPBD::InitDistanceIndices(int resY, int resX)
{
	int num = resY * resX;

	for (int i = 0; i < num - resX; i++)
	{
		if (i % resX == 0)
		{
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + 1);
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + resX);
		}
		else if (i % resX == resX - 1)
		{
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + resX);
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + resX - 1);
		}
		else
		{
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + 1);
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + resX);
			topol.indices.m_Data.push_back(i);
			topol.indices.m_Data.push_back(i + resX - 1);
		}
	}
	for (int i = num - resX; i < num - 1; ++i)
	{
		topol.indices.m_Data.push_back(i);
		topol.indices.m_Data.push_back(i + 1);
	}
}

void ConstraintPBD::InitDistanceInfo(BufferVector3f& meshPosBuffer, float stiffness)
{
	//
	// constraint allocation 
	// 
	int count = 2;
	int prims = topol.indices.GetSize() / count;
	//cout << __FUNCTION__ <<"   prims:"<< prims << std::endl;
	for (int i = 0; i < prims; i++)
	{
		// init primList & sortedPrimID 
		int2 p;
		p.x = i * count;
		p.y = count;
		topol.primList.m_Data.push_back(p);
		sortedPrimId.m_Data.push_back(i); // WARNING: sortedPrimID JUST for now (without collision)
		// init stiffness
		stiffnessBuffer.m_Data.push_back(stiffness);
		// init rest length
		int i0 = topol.indices.m_Data[p.x];
		int i1 = topol.indices.m_Data[p.x + 1];
		float d = Distance(meshPosBuffer.m_Data[i0], meshPosBuffer.m_Data[i1]);
		restLengthBuffer.m_Data.push_back(d);
		// init contraint type
		constraintType.m_Data.push_back(DISTANCE);
		// init color
		color.m_Data.push_back(-1);
		prdColor.m_Data.push_back(-1);
	}

	/*
	for (int a = 0; a < 2; a++)
	{
		int2 p;
		p.x = indices.GetSize();
		p.y = 1;

		constrPBDBuffer.color.m_Data.push_back(-1);
		constrPBDBuffer.prdColor.m_Data.push_back(-1);
		constrPBDBuffer.sortedColor.m_Data.push_back(-1);
		constrPBDBuffer.constraintType.m_Data.push_back(ANCHOR);
		constrPBDBuffer.stiffness.m_Data.push_back(1.0);

		indices.m_Data.push_back(a == 0 ? 0 : resY - 1);
		primList.m_Data.push_back(p);
		sortedPrimID.m_Data.push_back(prims);

		std::cout << "primList:" << primList.GetSize() << std::endl;
		std::cout << "indices:" << indices.m_Data[indices.GetSize()] << std::endl;
	}
	*/

}

void ConstraintPBD::InitBendingConstr()
{
	//todo
	//init restAngleBuffer;
		/*
	for (auto edge : edge2primMap)
	{
		//edge2primMap numEdge*2
		// ��������
		// edge 2 prim
		// int2 List;
	}
	*/
}

void ConstraintPBD::InitAnchorConstr(BufferVector3f& meshPosBuffer, float stiffness, int resY)
{

	int currPrimSize = sortedPrimId.GetSize();
	//	std::cout << __FUNCTION__ << std::endl;
	for (int a = 0; a < 2; a++)
	{
		int idx = (a == 0 ? 0 : resY - 1);
		// init primList
		int2 p;
		p.x = topol.indices.GetSize();
		p.y = 1;
		topol.primList.m_Data.push_back(p);
		// init indices
		topol.indices.m_Data.push_back(idx);
		// init restPosBuffer
		restPosBuffer.m_Data.push_back(meshPosBuffer.m_Data[idx]);
		// init stiffnessBuffer
		stiffnessBuffer.m_Data.push_back(stiffness);
		// init constraintType
		constraintType.m_Data.push_back(ANCHOR);
		//init color prdColor sortedPrimID
		color.m_Data.push_back(stiffness);
		prdColor.m_Data.push_back(stiffness);
		sortedPrimId.m_Data.push_back(currPrimSize + a);
		/*std::cout << "primList:" << primList.GetSize() << std::endl;
		std::cout << "indices:" << indices.m_Data[indices.GetSize()] << std::endl;*/
	}
}

void ConstraintPBD::GenePoint2PrimsMap(Topology topol)
{
	auto primList = &(topol.primList);
	auto indices = &(topol.indices);
	for (int primId = 0; primId < primList->GetSize(); ++primId)
	{
		int2 currPrim = primList->m_Data[primId];
		for (int i = 0; i < currPrim.y; ++i)
		{
			Point2PrimsMap[indices->m_Data[currPrim.x + i]].push_back(primId);
		}
	}
	// printf("generated point2prims map\n");
}

void ConstraintPBD::GenePrim2PrimsMap(Topology topol)
{
	auto primList = &(topol.primList);
	auto indices = &(topol.indices);
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
		int startIdx = Prim2PrimsMap.indices.GetSize();
		Prim2PrimsMap.indices.m_Data.insert(
			std::end(Prim2PrimsMap.indices.m_Data),
			std::begin(linkedPrimsSet),
			std::end(linkedPrimsSet));
		Prim2PrimsMap.startNumList.m_Data.push_back(make_int2(startIdx, linkedPrimsSet.size()));
	}
	//printf("generated prim2prims map\n");
}

void ConstraintPBD::AssignColorsCPU()
{
	auto p2pIndices = &(Prim2PrimsMap.indices);
	auto p2pStartNumList = &(Prim2PrimsMap.startNumList);
	for (int idx = 0; idx < topol.primList.GetSize(); idx++)
	{
		if (idx == 0)
			detectConflict = 0;

		if (color.m_Data[idx] >= 0)
			continue;

		int nidx = p2pStartNumList->m_Data[idx].x;
		//int nlast = p2pStartNumList->m_Data[idx].x + p2pStartNumList->m_Data[idx].y;
		int nlast = p2pStartNumList->m_Data[idx].x + p2pStartNumList->m_Data[idx].y;
		int c = -1, offset = 0;
		while (c < 0)
		{
			// Bit flag to store seen neighbor colors.
			unsigned long forbidden = 0;
			for (int i = nidx; i < nlast; ++i)
			{
				int n = p2pIndices->m_Data[i];
				int nc = color.m_Data[n] - offset;
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
		prdColor.m_Data[idx] = c;
	}
}

void ConstraintPBD::ResolveConflictsCPU()
{
	auto primList = &(topol.primList);
	auto p2pIndices = &(Prim2PrimsMap.indices);
	auto p2pStartNumList = &(Prim2PrimsMap.startNumList);

	for (int idx = 0; idx < primList->GetSize(); ++idx)
	{
		// Nothing to do if already colored.
		if (color.m_Data[idx] >= 0)
			continue;

		int nidx = p2pStartNumList->m_Data[idx].x;
		int nlast = p2pStartNumList->m_Data[idx].x + p2pStartNumList->m_Data[idx].y;
		int c = prdColor.m_Data[idx];
		//int c = newcolor[idx];
		int conflict = 0;
		int npt = primList->m_Data[idx].y;

		for (int i = nidx; i < nlast; ++i)
		{
			int n = p2pIndices->m_Data[i];
			int pc = prdColor.m_Data[n];

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
			color.m_Data[idx] = c;

	}
}

// kernel functions for PBDObj class
void __global__ AssignColorsGPU(
	int2* topolPrimList,
	int* p2pIndices,
	int2* p2pStartNumList,
	int* color,
	int* prdColor,
	int* detectConflict,
	int worksetLength)
{

	int primId = blockIdx.x * blockDim.x + threadIdx.x;
	// if(primId < 5) printf("%d: entering kernel1\n", primId);

	if (primId >= worksetLength)
		return;

	if (primId == 0)
		*detectConflict = 0;

	if (color[primId] >= 0)
		return;

	int nidx = p2pStartNumList[primId].x;
	int nlast = p2pStartNumList[primId].x + p2pStartNumList[primId].y;
	int c = -1, offset = 0;
	while (c < 0)
	{
		// Bit flag to store seen neighbor colors.
		unsigned long forbidden = 0;
		for (int i = nidx; i < nlast; ++i)
		{
			int n = p2pIndices[i];
			int nc = color[n] - offset;
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
	prdColor[primId] = c;
	// if (primId < 5) printf("\tAssignColorsGPU-%d: prdColor %d\n", primId, prdColor[primId]);
}

void __global__ ResolveConflictsGPU(
	int2* topolPrimList,
	int* p2pIndices,
	int2* p2pStartNumList,
	int* color,
	int* prdColor,
	int* detectConflict,
	int worksetLength)
{

	int primId = blockIdx.x * blockDim.x + threadIdx.x;
	// if(primId < 5) printf("\tResolveConflictsGPU-%d: prdColor %d\n", primId, prdColor[primId]);

	if (primId >= worksetLength)
		return;

	if (color[primId] >= 0)
		return;

	int nidx = p2pStartNumList[primId].x;
	int nlast = p2pStartNumList[primId].x + p2pStartNumList[primId].y;
	int c = prdColor[primId];

	int conflict = 0;
	int npt = topolPrimList[primId].y;

	for (int i = nidx; i < nlast; ++i)
	{
		int n = p2pIndices[i];
		int pc = prdColor[n];

		// Check for conflict.
		if (pc == c)
		{
			int nnpt = topolPrimList[n].y;
			// Resolution gives preference to primitives with more points, 
			// otherwise the one that comes first.
			// (NOTE: We color in fewer iterations if we prioritize by number
			// of graph neighbors, but then assuming we process prims by color,
			// we're usually farther away from original point order leading to many
			// cache misses and slower downstream processing.)
			if (nnpt > npt ||
				(nnpt == npt && n < primId))
			{
				conflict = 1;
				break;
			}
		}
	}

	// If there's a conflict then reset sizes for more work,
	// otherewise accept speculative color.
	if (conflict)
		*detectConflict = worksetLength;
	else
		color[primId] = c;

}

void ConstraintPBD::EdgeColoring(int iterations)
{
	switch (ht)
	{
	case CPU:
		EdgeColoringCPU(iterations);
		break;
	case GPU:
		EdgeColoringGPU(iterations);
		break;
	default:
		break;
	}
}

void ConstraintPBD::EdgeColoringCPU(int iterations)
{
	for (int i = 0; i < iterations; ++i)
	{
		/*printf("Before iteration %d color: ", i);
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
		printf("\n");*/
		AssignColorsCPU();
		ResolveConflictsCPU();
		/*printf("After iteration %d color: ", i);
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
		printf("\n");*/
	}
}

void ConstraintPBD::EdgeColoringGPU(int iterations)
{
	auto primList = &(topol.primList);
	int primNum = primList->GetSize();
	auto p2pIndices = &(Prim2PrimsMap.indices);
	auto p2pStartNumList = &(Prim2PrimsMap.startNumList);

	uint2 blockSize = primList->EvalBlockSize(512);
	//printf("Edge Color GPU: block dim = %d, thread dim = %d\n", blockSize.x, blockSize.y);

	//// Host To Device: load buffers for Edge Coloring
	//primList->MallocAndLoadToDevice();  // topolPrimList
	//p2pIndices->MallocAndLoadToDevice();  // p2pIndices
	//p2pStartNumList->MallocAndLoadToDevice();  // p2pStartNumList
	//color->MallocAndLoadToDevice();  // color
	//prdColor->MallocAndLoadToDevice();  // prdColor
	//
	//cudaMalloc((void **)&dFlag, sizeof(int));
	//cudaError_t cudaStatus = cudaMemcpy(dFlag, &detectConflict, sizeof(int), cudaMemcpyHostToDevice);
	//if(cudaStatus == cudaSuccess)  printf("\nMalloc and loads succeeded!\n");
	// edge coloring SPMD
	for (int i = 0; i < iterations; ++i)
	{
		AssignColorsGPU << < blockSize.x, blockSize.y >> > ((int2*)primList->GetDevicePtr(),
			(int*)p2pIndices->GetDevicePtr(),
			(int2*)p2pStartNumList->GetDevicePtr(),
			(int*)color.GetDevicePtr(),
			(int*)prdColor.GetDevicePtr(),
			dFlag,
			primNum);
		ResolveConflictsGPU << < blockSize.x, blockSize.y >> > ((int2*)primList->GetDevicePtr(),
			(int*)p2pIndices->GetDevicePtr(),
			(int2*)p2pStartNumList->GetDevicePtr(),
			(int*)color.GetDevicePtr(),
			(int*)prdColor.GetDevicePtr(),
			dFlag,
			primNum);
	}
	// Device To Host: load buffers back
	color.LoadToHost();
	/*if (color.LoadToHost())
		printf("Edge Color GPU: Load color Back Succeeded!\n");*/

		/*printf("\nAfter iteration(color): ");
		for (int i = 0; i < constrPBDBuffer.color.GetSize(); ++i)
		{
			printf("%d // ", constrPBDBuffer.color.m_Data[i]);
		}
		printf("\n");*/

		// Free GPU memory
		/*cudaFree(primList->GetDevicePtr());
		cudaFree(p2pIndices->GetDevicePtr());
		cudaFree(p2pStartNumList->GetDevicePtr());
		cudaFree(color->GetDevicePtr());
		cudaFree(prdColor->GetDevicePtr());
		cudaFree(dFlag);*/
}

void ConstraintPBD::SortEdgesColors()
{
	// cout << "--------" << __FUNCTION__ << "--------" <<  endl;
	/*for (int i = 0; i < color->GetSize(); ++i)
	{
		printf("%d - ", color->m_Data[i]);
	}
	printf("\n");
	for (int i = 0; i < sortedPrimID->GetSize(); ++i)
	{
		printf("%d - ", sortedPrimID->m_Data[i]);
	}
	printf("\n");*/
	// cout << __FUNCDNAME__ << endl;
	thrust::sort_by_key(color.m_Data.begin(), color.m_Data.end(), sortedPrimId.m_Data.begin());
	auto dColor = &(color);
	auto dSortedPrimId = &(sortedPrimId);
	dColor->LoadToDevice();
	dSortedPrimId->LoadToDevice();
	/*for (int i = 0; i < color.GetSize(); ++i)
	{
		printf("color : %d - ", color.m_Data[i]);
	}
	printf("\n");
	for (int i = 0; i < sortedPrimId.GetSize(); ++i)
	{
		printf("sortedPrimId: %d - ", sortedPrimId.m_Data[i]);
	}
	printf("\n");*/
	/*for (int i = 0; i < sortedColor->GetSize(); ++i)
	{
		printf("%d - ", sortedColor->m_Data[i]);
	}
	printf("\n");

	printf("\n");
	for (int i = 0; i < color->GetSize(); ++i)
	{
		printf("%d - ", color->m_Data[i]);
	}
	printf("\n");*/
}

void ConstraintPBD::EvalWorksets()
{
	//cout << "--------" << __FUNCTION__ << "--------" << endl;
	colorWorksets.m_Data.clear();
	int count = 1;
	for (int i = 1; i < color.GetSize(); ++i)
	{
		if (i == color.GetSize() - 1 && color.m_Data[i] == color.m_Data[i - 1])
		{
			count++;
			colorWorksets.m_Data.push_back(make_int2(i - count + 1, count));
		}
		else if (i == color.GetSize() - 1 && color.m_Data[i] != color.m_Data[i - 1])
		{
			colorWorksets.m_Data.push_back(make_int2(i - count, count));
			colorWorksets.m_Data.push_back(make_int2(i, 1));
		}
		else if (i != color.GetSize() - 1 && color.m_Data[i] != color.m_Data[i - 1])
		{
			colorWorksets.m_Data.push_back(make_int2(i - count, count));
			count = 1;
		}
		else
		{
			count++;
		}
	}
	/*for (int i = 0; i < sortedColor->GetSize(); ++i)
	{
		printf("%d - ", sortedColor->m_Data[i]);
	}*/
	/*for (int i = 0; i < colorWorksets.GetSize(); ++i)
	{
		printf("start: %d, num: %d  ", colorWorksets.m_Data[i].x, colorWorksets.m_Data[i].y);
	}
	printf("\n");*/
}

void ConstraintPBD::Save(std::ofstream& ofs)
{
	topol.Save(ofs);
	prdPBuffer.SetName("prdPBuffer");
	IO::SaveBuffer(prdPBuffer, ofs);
	restPosBuffer.SetName("restPosBuffer"); // empty now, initial in initPosition which is not used in read data from houdini
	IO::SaveBuffer(restPosBuffer, ofs);
	restLengthBuffer.SetName("restLengthBuffer");
	IO::SaveBuffer(restLengthBuffer, ofs);
	stiffnessBuffer.SetName("stiffnessBuffer");
	IO::SaveBuffer(stiffnessBuffer, ofs);
	constraintType.SetName("constraintType");
	IO::SaveBuffer(constraintType, ofs);
}

// ----------------------------------------------------------------------------------------
// PBDObject class
void PBDObject::Init()
{
	initMeshTopol();
	initConstr();
}

void PBDObject::Init(string topolFileName, string distConstrFileName)
{
	/*if ((IO::readTopolFromTxt(topolFileName, this)) && (IO::readDistConstrFromTxt(distConstrFileName, this)))
		printf("PBD Object was initialized successfully\n");*/
	bool readTopol = IO::ReadTopolFromTxt(topolFileName, this);
	bool readConstr = IO::ReadDistConstrFromTxt(distConstrFileName, this);
	if (readTopol && readConstr)
		printf("PBD Object was initialized successfully\n");
}

void PBDObject::ContinueSimInit(string meshTopolPath, string constrPath, HardwareType hardwareType)
{
	ht = hardwareType;
	bool readTopol = IO::ReadTopolFromCache(meshTopolPath, this);
	bool readConstr = IO::ReadConstraintFromCache(constrPath, this);
	if (readTopol && readConstr)
		printf("PBD Object was initialized successfully\n");
}

void PBDObject::SetConstrOption(uint ct, float* stiffnessSetting)
{
	this->ct = ct;
	this->stiffnessSetting = stiffnessSetting;
}

void PBDObject::initMeshTopol()
{
	// OpenGL Topology
	meshTopol.indices.SetName("Indices");
	meshTopol.posBuffer.SetName("P");
	meshTopol.primList.SetName("primList");

	initPosition(make_float2(0.0, 0.0));
	initMassVel();
	initMeshTopolIndices();
}

void PBDObject::initMassVel()
{
	// init mass
	float3 initVel = make_float3(0.0f, 0.0f, 0.0f);
	int num = resY * resX;
	for (int i = 0; i < num; i++)
	{
		massBuffer.m_Data.push_back(1.0);
		/*if (i == 0 || i == resY - 1)
		{
			massBuffer.m_Data.push_back(0.0);
		}
		else
		{
			massBuffer.m_Data.push_back(1.0);
		}*/
		velBuffer.m_Data.push_back(initVel);
	}
}

void PBDObject::InitGPUBuffers()
{
	// printf("init GPU buffers\n");
	auto dPrimList = &(constrPBDBuffer.topol.primList);
	auto dP2pIndices = &(constrPBDBuffer.Prim2PrimsMap.indices);
	auto dP2pStartNumList = &(constrPBDBuffer.Prim2PrimsMap.startNumList);
	auto dColor = &(constrPBDBuffer.color);
	auto dPrdColor = &(constrPBDBuffer.prdColor);
	auto dsortedPrimId = &(constrPBDBuffer.sortedPrimId);
	auto dConstrType = &(constrPBDBuffer.constraintType);

	auto dVelBuffer = &(velBuffer);
	auto dPrdPBuffer = &(constrPBDBuffer.prdPBuffer);
	auto dPositionBuffer = &(meshTopol.posBuffer);

	auto dMassBuffer = &(massBuffer);
	auto dRestLengthBuffer = &(constrPBDBuffer.restLengthBuffer);
	auto dStiffnessBuffer = &(constrPBDBuffer.stiffnessBuffer);
	auto dRestPosBuffer = &(constrPBDBuffer.restPosBuffer);
	auto dIndices = &(constrPBDBuffer.topol.indices);

	// Host To Device: load buffers for Edge Coloring
	cudaMalloc((void**)&constrPBDBuffer.dFlag, sizeof(int));
	cudaError_t cudaStatus = cudaMemcpy(constrPBDBuffer.dFlag, &(constrPBDBuffer.detectConflict), sizeof(int), cudaMemcpyHostToDevice);
	//if (cudaStatus == cudaSuccess)  printf("\nMalloc and loads succeeded!\n");

	dPrimList->MallocAndLoadToDevice();  // topolPrimList
	dP2pIndices->MallocAndLoadToDevice();  // p2pIndices
	dP2pStartNumList->MallocAndLoadToDevice();  // p2pStartNumList
	dColor->MallocAndLoadToDevice();  // color
	dPrdColor->MallocAndLoadToDevice();  // prdColor
	dConstrType->MallocAndLoadToDevice();  // constrType
	dsortedPrimId->DeviceMalloc();

	dVelBuffer->MallocAndLoadToDevice(); // velocity Buffer
	dPrdPBuffer->MallocAndLoadToDevice();  // predicted position buffer
	dPositionBuffer->MallocAndLoadToDevice();  // point real position buffer

	dMassBuffer->MallocAndLoadToDevice();
	dRestLengthBuffer->MallocAndLoadToDevice();
	dStiffnessBuffer->MallocAndLoadToDevice();
	dRestPosBuffer->MallocAndLoadToDevice();
	dIndices->MallocAndLoadToDevice();
}

void PBDObject::freeGPUBuffers()
{
	auto dPrimList = &(constrPBDBuffer.topol.primList);
	auto dP2pIndices = &(constrPBDBuffer.Prim2PrimsMap.indices);
	auto dP2pStartNumList = &(constrPBDBuffer.Prim2PrimsMap.startNumList);
	auto dColor = &(constrPBDBuffer.color);
	auto dPrdColor = &(constrPBDBuffer.prdColor);
	auto dSortedPrimId = &(constrPBDBuffer.sortedPrimId);
	auto dConstrType = &(constrPBDBuffer.constraintType);

	auto dVelBuffer = &(velBuffer);
	auto dPrdPBuffer = &(constrPBDBuffer.prdPBuffer);
	auto dPositionBuffer = &(meshTopol.posBuffer);

	auto dMassBuffer = &(massBuffer);
	auto dRestLengthBuffer = &(constrPBDBuffer.restLengthBuffer);
	auto dStiffnessBuffer = &(constrPBDBuffer.stiffnessBuffer);
	auto dRestPosBuffer = &(constrPBDBuffer.restPosBuffer);
	auto dIndices = &(constrPBDBuffer.topol.indices);

	dPositionBuffer->LoadToHost();
	dVelBuffer->LoadToHost();

	cudaFree(constrPBDBuffer.dFlag);

	cudaFree(dPrimList->GetDevicePtr());
	cudaFree(dP2pIndices->GetDevicePtr());
	cudaFree(dP2pStartNumList->GetDevicePtr());
	cudaFree(dColor->GetDevicePtr());
	cudaFree(dPrdColor->GetDevicePtr());
	cudaFree(dSortedPrimId->GetDevicePtr());
	cudaFree(dConstrType->GetDevicePtr());

	cudaFree(dVelBuffer->GetDevicePtr());
	cudaFree(dPrdPBuffer->GetDevicePtr());
	cudaFree(dPositionBuffer->GetDevicePtr());

	cudaFree(dMassBuffer->GetDevicePtr());
	cudaFree(dRestLengthBuffer->GetDevicePtr());
	cudaFree(dStiffnessBuffer->GetDevicePtr());
	cudaFree(dRestPosBuffer->GetDevicePtr());
	cudaFree(dIndices->GetDevicePtr());
}

// init : allocation
// setvalue
void PBDObject::initConstr()
{
	constrPBDBuffer.ht = ht;
	constrPBDBuffer.topol.posBuffer = meshTopol.posBuffer;
	constrPBDBuffer.prdPBuffer = meshTopol.posBuffer;

	constrPBDBuffer.topol.primList.SetName("primList");
	constrPBDBuffer.topol.indices.SetName("Indices");
	constrPBDBuffer.color.SetName("color");
	constrPBDBuffer.prdColor.SetName("prdcolor");
	constrPBDBuffer.sortedPrimId.SetName("sortedPrimId");
	if ((DISTANCE & ct) == DISTANCE)
	{
		constrPBDBuffer.InitDistanceConstr(meshTopol.posBuffer, stiffnessSetting[0], resY, resX);
	}
	if ((BENDING & ct) == BENDING)
	{
		constrPBDBuffer.InitBendingConstr();
	}
	if ((ANCHOR & ct) == ANCHOR)
	{
		constrPBDBuffer.InitAnchorConstr(meshTopol.posBuffer, -1.0f, resY);
	}
	constrPBDBuffer.GenePoint2PrimsMap(constrPBDBuffer.topol);
	constrPBDBuffer.GenePrim2PrimsMap(constrPBDBuffer.topol);
	if (ht == GPU)
	{
		InitGPUBuffers();
		constrPBDBuffer.EdgeColoring(20000);
	}

}

void PBDObject::initPosition(float2 cord)
{
	auto positionBuffer = &(meshTopol.posBuffer);
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
			positionBuffer->m_Data.push_back(p);
			index++;
		}
	}
	constrPBDBuffer.restPosBuffer.m_Data.push_back(positionBuffer->m_Data[0]);
	constrPBDBuffer.restPosBuffer.m_Data.push_back(positionBuffer->m_Data[resX - 1]);
}

void PBDObject::initMeshTopolIndices()
{
	auto meshTopolIndicies = &(meshTopol.indices);
	int num = resY * resX;
	for (int i = 0; i < num - resX; i++)
	{
		if (i % resX == resX - 1)
			continue;
		meshTopolIndicies->m_Data.push_back(i);
		meshTopolIndicies->m_Data.push_back(i + resX);
		meshTopolIndicies->m_Data.push_back(i + resX + 1);
		meshTopolIndicies->m_Data.push_back(i);
		meshTopolIndicies->m_Data.push_back(i + resX + 1);
		meshTopolIndicies->m_Data.push_back(i + 1);
	}
}

void PBDObject::groundTruthTest()
{
	vector<int> arr0 = { 0, 1, 3, 3, 2, 0 };
	vector<int> arr1 = { 0,3,0,1,0,2,1,3,2,3 };
	constrPBDBuffer.topol.indices.m_Data = arr1;
	vector<int2> arr2 = { make_int2(0,2),  make_int2(2,2),  make_int2(4,2),  make_int2(6,2),  make_int2(8,2) };
	constrPBDBuffer.topol.primList.m_Data = arr2;
	vector<int> arr3 = { -1,-1,-1,-1,-1 };
	constrPBDBuffer.color.m_Data = arr3;
	constrPBDBuffer.prdColor.m_Data = arr3;
}

void PBDObject::Save(string path)
{
	PBD_DEBUG;
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;

	ofs << "Header|float3Buffer,6;floatBuffer,2;int2Buffer,2;intBuffer,3;float3,1;float,1" << endl; //HEADER
	meshTopol.indices.SetName("meshTopol indices");
	meshTopol.posBuffer.SetName("meshTopol posBuffer");
	meshTopol.primList.SetName("meshTopol primList");
	meshTopol.Save(ofs);
	constrPBDBuffer.Save(ofs);
	velBuffer.SetName("velBuffer");
	IO::SaveBuffer(velBuffer, ofs);
	massBuffer.SetName("massBuffer");
	IO::SaveBuffer(massBuffer, ofs);
	IO::SaveData(dampingRate, "dampingRate", ofs);
	IO::SaveData(gravity, "gravity", ofs);

	ofs.flush();
	ofs.close();
}

void PBDObject::SaveMeshTopol(string path)
{
	PBD_DEBUG;
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;

	//ofs << "Header|float3Buffer,3;int2Buffer,1;intBuffer,1;float3,1;float,1" << endl; //HEADER
	meshTopol.Save(ofs);
	velBuffer.SetName("velBuffer");
	IO::SaveBuffer(velBuffer, ofs);
	massBuffer.SetName("massBuffer");
	IO::SaveBuffer(massBuffer, ofs);
	IO::SaveData(dampingRate, "dampingRate", ofs);
	IO::SaveData(gravity, "gravity", ofs);

	ofs.flush();
	ofs.close();
}

void PBDObject::SaveConstraint(string path)
{
	PBD_DEBUG;
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;

	//ofs << "Header|float3Buffer,3;floatBuffer,2;int2Buffer,1;intBuffer,2" << endl; //HEADER
	constrPBDBuffer.Save(ofs);

	ofs.flush();
	ofs.close();
}

void PBDObject::Read(string path)
{

}
// ----------------------------------------------------------------------------------------
// SolverPBD class

// kernel functions for SolverPBD class
void __global__ AdvectGPUKernel(
	int pointNum,
	float dt,
	float dampingRate,
	float3 gravity,
	float* mass,
	float3* velBuffer,
	float3* prdPBuffer,
	float3* positionBuffer)
{
	int pointId = blockIdx.x * blockDim.x + threadIdx.x;
	// if(primId < 5) printf("\tResolveConflictsGPU-%d: prdColor %d\n", primId, prdColor[primId]);

	if (pointId >= pointNum)
		return;


	/*if (i == 30)
		printf("old velocity Buffer: %f, %f, %f \n", velBuffer->m_Data[i].x, velBuffer->m_Data[i].y, velBuffer->m_Data[i].z);*/
	velBuffer[pointId] += gravity * dt * mass[pointId];
	/*if(i == 30)
		printf("new velocity Buffer: %f, %f, %f \n", velBuffer->m_Data[i].x, velBuffer->m_Data[i].y, velBuffer->m_Data[i].z);*/
	velBuffer[pointId] *= powf(dampingRate, dt);
	prdPBuffer[pointId] = positionBuffer[pointId] + velBuffer[pointId] * dt;
	//printf("postion Buffer: %f, %f, %f \n", prdPBuffer.m_Data[j].x, prdPBuffer.m_Data[j].y, prdPBuffer.m_Data[j].z);

	/*printf("Advect: prdPBuffer: \n");
	printf("(%f,%f,%f)", prdPBuffer[pointId].x, prdPBuffer[pointId].y, prdPBuffer[pointId].z);*/


}

void SolverPBD::Advect(float dt)
{
	PBD_DEBUG;
	switch (m_ht)    // TODO: change back to ht
	{
	case CPU:
		advectCPU(dt);
		break;
	case GPU:
		advectGPU(dt);
		break;
	default:
		break;
	}
}

void SolverPBD::advectCPU(float dt)
{
	auto velBuffer = &(m_pbdObj->velBuffer);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto positionBuffer = &(m_pbdObj->meshTopol.posBuffer);

	for (int i = 0; i < velBuffer->GetSize(); i++)
	{
		velBuffer->m_Data[i] += m_pbdObj->gravity * dt;
		prdPBuffer->m_Data[i] = positionBuffer->m_Data[i] + velBuffer->m_Data[i] * dt;

	}
}

void SolverPBD::advectGPU(float dt)
{
	auto velBuffer = &(m_pbdObj->velBuffer);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto positionBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto massBuffer = &(m_pbdObj->massBuffer);

	//printf("Before Advect GPU:");
	//printf("point 0: %f, %f, %f; point col-1: %f, %f, %f\n",
	//	prdPBuffer->m_Data[0].x, prdPBuffer->m_Data[0].y, prdPBuffer->m_Data[0].z,
	//	prdPBuffer->m_Data[pbdObj->resY-1].x, prdPBuffer->m_Data[pbdObj->resProjectConstraintGPUY - 1].y,
	//	prdPBuffer->m_Data[pbdObj->resY - 1].z );

	int pointNum = prdPBuffer->GetSize();
	float dampingRate = m_pbdObj->dampingRate;
	float3 gravity = m_pbdObj->gravity;

	uint2 blockSize = positionBuffer->EvalBlockSize(512);
	//printf("Advect GPU: block dim = %d, thread dim = %d\n", blockSize.x, blockSize.y);
	// Host To Device: load buffers for Advect

	// edge coloring SPMD
	AdvectGPUKernel << < blockSize.x, blockSize.y >> > (pointNum,
		dt,
		dampingRate,
		gravity,
		(float*)massBuffer->GetDevicePtr(),
		(float3*)velBuffer->GetDevicePtr(),
		(float3*)prdPBuffer->GetDevicePtr(),
		(float3*)positionBuffer->GetDevicePtr());
	// Device To Host: load buffers back
	//if (prdPBuffer->LoadToHost())
	//	printf("Advect GPU: Load prdPBuffer Back Succeeded!\n");
	//if (velBuffer->LoadToHost())
	//	printf("Advect GPU: Load velBuffer Back Succeeded!\n");

	//printf("After Advect GPU:");
	//printf("point 0: %f, %f, %f; point col-1: %f, %f, %f\n",
	//	prdPBuffer->m_Data[0].x, prdPBuffer->m_Data[0].y, prdPBuffer->m_Data[0].z,
	//	prdPBuffer->m_Data[pbdObj->resY - 1].x, prdPBuffer->m_Data[pbdObj->resY - 1].y,
	//	prdPBuffer->m_Data[pbdObj->resY - 1].z);
}

void SolverPBD::ProjectConstraint(SolverType st, int iterations)
{
	m_pbdSolverTimer->Tick();
	switch (m_ht)   
	{
	case CPU:
		projectConstraintCPU(st, iterations);
		break;
	case GPU:
		projectConstraintGPU(st, iterations);
		break;
	default:
		break;
	}
	m_pbdSolverTimer->Tock();
	PBD_DEBUGTIME(m_pbdSolverTimer->GetFuncTime());
}

void SolverPBD::ProjectConstraintWithColli(SolverType st, int iterations, CollisionSolver* colliSolver, 
	BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, int debug)
{
	m_pbdSolverTimer->Tick();
	switch (m_ht)
	{
	case CPU:
		projectConstraintWithColliCPU(st, iterations, colliSolver, fixedBuffer, vFixedBuffer, fFixedBuffer, debug);
		break;
	case GPU:
		projectConstraintGPU(st, iterations);
		break;
	default:
		break;
	}
	m_pbdSolverTimer->Tock();
}

void SolverPBD::projectConstraintCPU(SolverType st, int iterations)
{
	auto primList = &(m_pbdObj->constrPBDBuffer.topol.primList);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto massBuffer = &(m_pbdObj->massBuffer);
	auto restLengthBuffer = &(m_pbdObj->constrPBDBuffer.restLengthBuffer);
	auto stiffnessBuffer = &(m_pbdObj->constrPBDBuffer.stiffnessBuffer);
	// auto restPosBuffer = &(m_pbdObj->constrPBDBuffer.restPosBuffer);
	auto indices = &(m_pbdObj->constrPBDBuffer.topol.indices);
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
			float3 r = prdPBuffer->m_Data[i0] - prdPBuffer->m_Data[i1];
			r = normalize(r);

			dp1.x = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) *r.x ;
			dp1.y = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) *r.y ;
			dp1.z = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) *r.z ;
			dp2.x = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i])  *r.x ;
			dp2.y = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i])  *r.y ;
			dp2.z = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i])  *r.z ;
			float k = 1;// -powf(1 - stiffnessBuffer->m_Data[i], 1.0 / (ii + 1));
			dp1 *= k;
			dp2 *= k;
			prdPBuffer->m_Data[i0] += dp1;
			prdPBuffer->m_Data[i1] += dp2;
		}
		ColliWithShpGrd();
	}
}


void SolverPBD::projectConstraintWithColliCPU(SolverType st, int iterations, CollisionSolver* colliSolver, 
	BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, int debug)
{
	PBD_DEBUG;
	int debugFrameId = 1;
	auto primList = &(m_pbdObj->constrPBDBuffer.topol.primList);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto massBuffer = &(m_pbdObj->massBuffer);
	auto restLengthBuffer = &(m_pbdObj->constrPBDBuffer.restLengthBuffer);
	auto stiffnessBuffer = &(m_pbdObj->constrPBDBuffer.stiffnessBuffer);
	// auto restPosBuffer = &(m_pbdObj->constrPBDBuffer.restPosBuffer);
	auto indices = &(m_pbdObj->constrPBDBuffer.topol.indices);
	colliSolver->afterProjPrdpBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
	for (size_t ii = 0; ii < iterations; ii++)
	{
		if ((ii % 10 ==0) || (ii == iterations - 1)) // 0 10 20 30 || == -1
		{
			colliSolver->CCD_SH();
			printf("contact size: %d\n", colliSolver->contactData.ctxs.GetSize());			
		}
		//printInfo("--- in project", prdPBuffer->m_Data[1]);
		for (size_t i = 0; i < primList->GetSize(); i++)
		{
			if (primList->m_Data[i].y != 2)
				continue;
			int i0 = indices->m_Data[primList->m_Data[i].x];
			int i1 = indices->m_Data[primList->m_Data[i].x + 1];
			float3 dp1;
			float3 dp2;
			float d = Distance(prdPBuffer->m_Data[i0], prdPBuffer->m_Data[i1]);
			float3 r = prdPBuffer->m_Data[i0] - prdPBuffer->m_Data[i1];
			r = normalize(r);

			dp1.x = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.x;
			dp1.y = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.y;
			dp1.z = -massBuffer->m_Data[i0] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.z;
			dp2.x = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.x;
			dp2.y = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.y;
			dp2.z = massBuffer->m_Data[i1] / (massBuffer->m_Data[i0] + massBuffer->m_Data[i1]) * (d - restLengthBuffer->m_Data[i]) * r.z;
			float k = 1;// -powf(1 - stiffnessBuffer->m_Data[i], 1.0 / (ii + 1));
			dp1 *= k;
			dp2 *= k;
			prdPBuffer->m_Data[i0] += dp1;
			prdPBuffer->m_Data[i1] += dp2;
		}
		//colliSolver->CollisionResolve();

		string beforeResolvePath = "D://0326Test//testData//testBeforeResolve." + to_string(debug * iterations + ii) + ".cache";
		m_pbdObj->constrPBDBuffer.prdPBuffer.SetName("P");
		Topology tempBeforeResolve;
		tempBeforeResolve.indices = m_pbdObj->meshTopol.indices;
		tempBeforeResolve.primList = m_pbdObj->meshTopol.primList;
		tempBeforeResolve.posBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
		tempBeforeResolve.indices.SetName("Indices");
		tempBeforeResolve.primList.SetName("primList");
		tempBeforeResolve.posBuffer.SetName("P");
		IO::SaveToplogy(tempBeforeResolve, beforeResolvePath);

		colliSolver->CollisionResolveNew(fixedBuffer, vFixedBuffer, fFixedBuffer, (debug * iterations + ii), ii, debugFrameId);
		
		string path = "D://0326Test//testData//test." + to_string(debug*iterations + ii) + ".cache";
		m_pbdObj->constrPBDBuffer.prdPBuffer.SetName("P");
		Topology temp;
		temp.indices = m_pbdObj->meshTopol.indices;
		temp.primList = m_pbdObj->meshTopol.primList;
		temp.posBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
		temp.indices.SetName("Indices");
		temp.primList.SetName("primList");
		temp.posBuffer.SetName("P");
		IO::SaveToplogy(temp, path);

		//printInfo("--- after resolve", prdPBuffer->m_Data[1]);
		//printf("--------------------itreation %d-------------------\n", ii);
		ColliWithShpGrd();
	}
}

// Attach Points
//for (size_t j = 0; j < prdPBuffer->GetSize(); j++)
//{
//	//attach points
//	if (j == 0)
//	{
//		prdPBuffer->m_Data[j] = restPosBuffer->m_Data[0];
//	}
//	if (j == m_pbdObj->resY - 1)
//	{
//		prdPBuffer->m_Data[j] = restPosBuffer->m_Data[1];
//	}

//	////point collide with sphere
//	//bool isCollideSphere = ColliderSphere(prdPBuffer.m_Data[j], sphereOrigin, sphereRadius, j);
//	//if (isCollideSphere) //move the point to the point which intersect with sphere
//	//{
//	//	float3 moveVector = GenerateMoveVectorSphere(sphereOrigin, sphereRadius, prdPBuffer.m_Data[j], j);
//	//	prdPBuffer.m_Data[j] += moveVector;
//	//}
//	////point collide with ground
//	//bool isCollideGoround = CollideGround(prdPBuffer.m_Data[j], groundCenter);
//	//if (isCollideGoround)
//	//{
//	//	prdPBuffer.m_Data[j].y = groundCenter.y;
//	//}
//}

void __global__ ProjectContraintsGPUKernel(
	int resY,
	int start,
	int num,
	int iteration,
	int* sortedPrimId,
	int* indices,
	int* constraintType,
	float* massBuffer,
	float* restLengthBuffer,
	float* stiffnessBuffer,
	int2* primList,
	float3* prdPBuffer,
	float3* restBuffer)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= num)
		return;
	//printf("idx: %d ---- ", idx + start);
	auto primId = sortedPrimId[idx + start];
	//printf("idx + start: %d, primId: %d \n", idx + start, primId);


	if (constraintType[primId] == DISTANCE)
	{
		int i0 = indices[primList[primId].x];
		int i1 = indices[primList[primId].x + 1];

		float3 p0 = prdPBuffer[i0];
		float3 p1 = prdPBuffer[i1];
		/*printf("Project Constraint: prdPBuffer: ");
		printf("i0: (%f,%f,%f)", prdPBuffer[i0].x, prdPBuffer[i0].y, prdPBuffer[i0].z);
		printf("i1: (%f,%f,%f)\n", prdPBuffer[i1].x, prdPBuffer[i1].y, prdPBuffer[i1].z);*/
		float3 dp0;
		float3 dp1;
		float d = Distance(p0, p1);
		float3 v;
		v = p0 - p1;
		//printf("mass: %f", massBuffer[i0]);
		dp0 = -massBuffer[i0] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[primId]) * v / d;
		dp1 = massBuffer[i1] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[primId]) * v / d;
		float k = 1 - powf(1 - stiffnessBuffer[primId], 1.0 / (iteration + 1));
		/*printf("dp0: (%f,%f,%f) ; ", dp0.x, dp0.y, dp0.z);
		printf("dp1: (%f,%f,%f)", dp1.x, dp1.y, dp1.z);
		printf("d: %f, k: %f \n",d, k);*/
		dp0 *= k;
		dp1 *= k;

		prdPBuffer[i0] += dp0;
		prdPBuffer[i1] += dp1;
		/*printf("Project Constraint: prdPBuffer: ");
		printf("i0: (%f,%f,%f)", prdPBuffer[i0].x, prdPBuffer[i0].y, prdPBuffer[i0].z);
		printf("i1: (%f,%f,%f)\n", prdPBuffer[i1].x, prdPBuffer[i1].y, prdPBuffer[i1].z);*/
	}



	if (constraintType[primId] == ANCHOR)
	{
		int i = indices[primList[primId].x];
		if (i == 0)
			prdPBuffer[i] = restBuffer[0];
		if (i == (resY - 1))
			prdPBuffer[i] = restBuffer[1];
	}
}

void SolverPBD::projectConstraintGPU(SolverType st, int iterations)
{

	m_pbdObj->constrPBDBuffer.SortEdgesColors();
	m_pbdObj->constrPBDBuffer.EvalWorksets();
	// cout << "--------" << __FUNCTION__ << "--------" << endl;

	auto worksets = &(m_pbdObj->constrPBDBuffer.colorWorksets);
	auto primList = &(m_pbdObj->constrPBDBuffer.topol.primList);
	auto sortedPrimId = &(m_pbdObj->constrPBDBuffer.sortedPrimId);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto massBuffer = &(m_pbdObj->massBuffer);
	auto restLengthBuffer = &(m_pbdObj->constrPBDBuffer.restLengthBuffer);
	auto stiffnessBuffer = &(m_pbdObj->constrPBDBuffer.stiffnessBuffer);
	auto restPosBuffer = &(m_pbdObj->constrPBDBuffer.restPosBuffer);
	auto indices = &(m_pbdObj->constrPBDBuffer.topol.indices);
	auto constraintType = &(m_pbdObj->constrPBDBuffer.constraintType);

	for (int i = 0; i < iterations; ++i)
	{
		for (auto workset : worksets->m_Data)
		{
			int start = workset.x;
			int num = workset.y;
			int numBlock = 1;
			int numThread = 512;
			if (num > numThread && num < 1024)
				numThread = num;
			else if (num > 1024)
				numBlock = ceil(num / 512);
			//printf("----------------------------\n");
			//printf("		numBlock: %d numThread: %d\n", numBlock, numThread);
			/*printf("indices: ");
			for (int i = 0; i < indices->GetSize(); ++i)
			{
				cout << indices->m_Data[i] << "-";
			}
			printf("\n");
			printf("constraintType: ");
			for (int i = 0; i < constraintType->GetSize(); ++i)
			{
				cout << constraintType->m_Data[i] << "-";
			}
			printf("\n");
			printf("massBuffer: ");
			for (int i = 0; i < massBuffer->GetSize(); ++i)
			{
				cout << massBuffer->m_Data[i] << "-";
			}
			printf("\n");
			printf("restLengthBuffer: ");
			for (int i = 0; i < restLengthBuffer->GetSize(); ++i)
			{
				cout << restLengthBuffer->m_Data[i] << "-";
			}
			printf("\n");
			printf("stiffnessBuffer: ");
			for (int i = 0; i < stiffnessBuffer->GetSize(); ++i)
			{
				cout << stiffnessBuffer->m_Data[i] << "-";
			}*/
			/*printf("prdPBuffer: \n");
			for (int i = 0; i < prdPBuffer->GetSize(); ++i)
			{
				cout << "(" << prdPBuffer->m_Data[i].x << "," << prdPBuffer->m_Data[i].y << "," << prdPBuffer->m_Data[i].z << ")" << endl;
			}
			printf("\n");
			printf("primList: \n");
			for (int i = 0; i < primList->GetSize(); ++i)
			{
				cout << "(" << primList->m_Data[i].x << "," << primList->m_Data[i].y << ")" << endl;
			}*/

			ProjectContraintsGPUKernel << <numBlock, numThread >> > (
				m_pbdObj->resY,
				start,
				num,
				i,
				((int*)sortedPrimId->GetDevicePtr()),
				(int*)indices->GetDevicePtr(),
				(int*)constraintType->GetDevicePtr(),
				(float*)massBuffer->GetDevicePtr(),
				(float*)restLengthBuffer->GetDevicePtr(),
				(float*)stiffnessBuffer->GetDevicePtr(),
				(int2*)primList->GetDevicePtr(),
				(float3*)prdPBuffer->GetDevicePtr(),
				(float3*)restPosBuffer->GetDevicePtr());
			/*cudaError_t cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			}*/
		}
	}

	//{
	//	// TODO: edge coloring
	//	/*
	//	int* color
	//	int* sortedColor
	//	int* sortedPrim
	//	int main
	//	{
	//		sort(color, sortedColor, sortedPrim)
	//		//color: 11220001110222
	//		//sortedColor:000111222
	//		//sortedPrim:
	//		int2 workSets = Eval(SortedColor)
	//		//workSets.x start
	//		//workSets.y num
	//		for(auto workset: workSets)
	//		{
	//			int start  = workSets.x
	//			int num = workSets.y
	//			kernel project<<<numBlock, numThread>>>(num, sortedPrim+start, prdPbuffer, restlength)
	//			{
	//				if(index > = num)
	//					return;
	//				int2 prim = sortedPrimId[index]
	//				int i0 = primList[prim.x]
	//				int i1 = primList[prim.x + 1]
	//				float3 p0 = prdPbuffer[i0]
	//				float3 p1 = prdPbuffer[i1]
	//			}
	//		}
	//	}
	//}
}

void SolverPBD::Integration(float dt)
{
	PBD_DEBUG;
	switch (m_ht)    // TODO: change back to ht
	{
	case CPU:
		integrationCPU(dt);
		break;
	case GPU:
		integrationGPU(dt);
		break;
	default:
		break;
	}
}

void SolverPBD::integrationCPU(float dt)
{
	auto positionBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto velBuffer = &(m_pbdObj->velBuffer);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);

	for (size_t i = 0; i < positionBuffer->GetSize(); i++)
	{
		velBuffer->m_Data[i] = (prdPBuffer->m_Data[i] - positionBuffer->m_Data[i]) / dt;
		positionBuffer->m_Data[i] = prdPBuffer->m_Data[i];

	}
}

void __global__ IntegrationGPUKernel(
	int pointNum,
	float dt,
	float3* velBuffer,
	float3* prdPBuffer,
	float3* positionBuffer)
{

	int pointId = blockIdx.x * blockDim.x + threadIdx.x;
	// if(primId < 5) printf("\tResolveConflictsGPU-%d: prdColor %d\n", primId, prdColor[primId]);

	if (pointId >= pointNum)
		return;
	velBuffer[pointId] = (prdPBuffer[pointId] - positionBuffer[pointId]) / dt;
	positionBuffer[pointId] = prdPBuffer[pointId];
}

// TODO
void SolverPBD::integrationGPU(float dt)
{
	auto positionBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto velBuffer = &(m_pbdObj->velBuffer);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	/*printf("Before Integration GPU:");
	printf("point 0: %f, %f, %f; point col-1: %f, %f, %f\n",
		prdPBuffer->m_Data[0].x, prdPBuffer->m_Data[0].y, prdPBuffer->m_Data[0].z,
		prdPBuffer->m_Data[pbdObj->resY - 1].x, prdPBuffer->m_Data[pbdObj->resY - 1].y,
		prdPBuffer->m_Data[pbdObj->resY - 1].z);*/
	int pointNum = prdPBuffer->GetSize();

	uint2 blockSize = positionBuffer->EvalBlockSize(512);
	//printf("Integration GPU: block dim = %d, thread dim = %d\n", blockSize.x, blockSize.y);

	IntegrationGPUKernel << < blockSize.x, blockSize.y >> > (pointNum,
		dt,
		(float3*)velBuffer->GetDevicePtr(),
		(float3*)prdPBuffer->GetDevicePtr(),
		(float3*)positionBuffer->GetDevicePtr());

	// Device To Host: load buffers back

	/*if (positionBuffer->LoadToHost())
		printf("Integration GPU: Load positionBuffer Back Succeeded!\n");
	if (velBuffer->LoadToHost())
		printf("Integration GPU: Load velBuffer Back Succeeded!\n");*/

		//printf("After Integration GPU:");
		//printf("point 0: %f, %f, %f; point col-1: %f, %f, %f\n",
		//	positionBuffer->m_Data[0].x, positionBuffer->m_Data[0].y, positionBuffer->m_Data[0].z,
		//	positionBuffer->m_Data[pbdObj->resY - 1].x, positionBuffer->m_Data[pbdObj->resY - 1].y,
		//	positionBuffer->m_Data[pbdObj->resY - 1].z);
}
void SolverPBD::ColliWithShpGrd()
{
	auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	for (int vtxId = 0; vtxId < posBuffer->GetSize(); ++vtxId)
	{
		//point collide with sphere
		bool isCollideSphere = ColliderSphere(prdPBuffer->m_Data[vtxId], m_sphereCenter, m_sphereRadius);
		if (isCollideSphere) //move the point to the point which intersect with sphere
		{
			float3 moveVector = GenerateMoveVectorSphere(m_sphereCenter, m_sphereRadius, prdPBuffer->m_Data[vtxId]);
			prdPBuffer->m_Data[vtxId] += moveVector;
		}
		//point collide with ground
		//bool isCollideGoround = CollideGround(prdPBuffer->m_Data[vtxId], m_groundHeight);
		//if (isCollideGoround)
		//{
		//	prdPBuffer->m_Data[vtxId].y = m_groundHeight;
		//}
	}

}

bool SolverPBD::ColliderSphere(float3 pointPos, float3 sphereOrigin, float r)
{
	float d = Distance(pointPos, sphereOrigin);
	if (d - r > 0.001)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool SolverPBD::CollideGround(float3 pointPos, float groundHeight)
{
	if (pointPos.y - groundHeight < 0.001)
	{
		return true;
	}
	else
	{
		return false;
	}
}

float3 SolverPBD::GenerateMoveVectorSphere(float3 sphereOrigin, float sphereRadius, float3  p)
{
	float moveDistance = sphereRadius - Distance(sphereOrigin, p);
	float3 moveDirection = (p - sphereOrigin) / Distance(sphereOrigin, p);
	float3 moveLength = moveDirection * moveDistance;
	return moveLength;
}
// ------------------Topology---------------------
void Topology::Save(std::ofstream& ofs)
{
	IO::SaveBuffer(indices, ofs);
	IO::SaveBuffer(posBuffer, ofs);
	IO::SaveBuffer(primList, ofs);
}
