#include "CCD_Basic.h"

float DistanceCCD(float3 p1, float3 p2)
{
	return powf(powf((p1.x - p2.x), 2) + powf((p1.y - p2.y), 2) + powf((p1.z - p2.z), 2), 0.5);
}

float3 pos(const float3 pos, const float3 prdPos, double t)
{
	return pos + t * (prdPos - pos);
}

bool RelativePos(float3 vtx, float3 pvtx, float3 n)
{
	float3 v = vtx - pvtx;
	float d = dot(v, n);
	if (d >= 0)
	{
		return true;
	}
	else if (d < 0)
	{
		return false;
	}
}

float Point2Plane(float3 vtx, float3 p1, float3 p2, float3 p3)
{
	//float3 r = vtx - p1;
	//float3 n = cross(p2 - p1, p3 - p1);
	//d = dot(r, n) / sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	float d;
	float3 AB = p2 - p1;
	float3 AC = p3 - p1;
	float3 N = cross(AB, AC);
	float Ax = N.x;
	float By = N.y;
	float Cz = N.z;
	float D = -(Ax * p1.x + By * p1.y + Cz * p1.z);
	float mod = Ax * vtx.x + By * vtx.y + Cz * vtx.z + D;
	float area = sqrt(Ax * Ax + By * By + Cz * Cz);
	d = abs(mod) / area;
	return d;
}

// Solving cubic equations
// solves a x^3 + b x^2 + c x + d == 0
int solve_cubic(double a, double b, double c, double d, double x[3]) {
	double xc[2];
	int ncrit = solve_quadratic(3 * a, 2 * b, c, xc);
	if (ncrit == 0) {
		x[0] = newtons_method(a, b, c, d, xc[0], 0);
		return 1;
	}
	else if (ncrit == 1) {// cubic is actually quadratic
		return solve_quadratic(b, c, d, x);
	}
	else {
		double yc[2] = { d + xc[0] * (c + xc[0] * (b + xc[0] * a)),
						d + xc[1] * (c + xc[1] * (b + xc[1] * a)) };
		int i = 0;
		if (yc[0] * a >= 0)
			x[i++] = newtons_method(a, b, c, d, xc[0], -1);
		if (yc[0] * yc[1] <= 0) {
			int closer = abs(yc[0]) < abs(yc[1]) ? 0 : 1;
			x[i++] = newtons_method(a, b, c, d, xc[closer], closer == 0 ? 1 : -1);
		}
		if (yc[1] * a <= 0)
			x[i++] = newtons_method(a, b, c, d, xc[1], 1);
		return i;
	}
}

double newtons_method(double a, double b, double c, double d, double x0, int init_dir)
{
	if (init_dir != 0) {
		// quadratic approximation around x0, assuming y' = 0
		double y0 = d + x0 * (c + x0 * (b + x0 * a)),
			ddy0 = 2 * b + x0 * (6 * a);
		x0 += init_dir * sqrt(abs(2 * y0 / ddy0));
	}
	for (int iter = 0; iter < 100; iter++) {
		double y = d + x0 * (c + x0 * (b + x0 * a));
		double dy = c + x0 * (2 * b + x0 * 3 * a);
		if (dy == 0)
			return x0;
		double x1 = x0 - y / dy;
		if (abs(x0 - x1) < 1e-6)
			return x0;
		x0 = x1;
	}
	return x0;
}

int solve_quadratic(double a, double b, double c, double x[2])
{
	// http://en.wikipedia.org/wiki/Quadratic_formula#Floating_point_implementation
	double d = b * b - 4 * a * c;
	if (d < 0) {
		x[0] = -b / (2 * a);
		return 0;
	}
	double q = -(b + sgn(b) * sqrt(d)) / 2;
	int i = 0;
	if (abs(a) > 1e-12 * abs(q))
		x[i++] = q / a;
	if (abs(q) > 1e-12 * abs(c))
		x[i++] = c / q;
	if (i == 2 && x[0] > x[1])
		swap(x[0], x[1]);
	return i;
}

float3 BarycentricCoord(float3 pos, float3 p0, float3 p1, float3 p2)
{
	float3 AB = p0 - p2;
	float3 AC = p1 - p2;
	float3 AP = pos - p2;
	float dot00 = dot(AC, AC);
	float dot01 = dot(AC, AB);
	float dot02 = dot(AC, AP);
	float dot11 = dot(AB, AB);
	float dot12 = dot(AB, AP);
	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);
	float v = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	float u = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	float w = 1 - u - v;
	return make_float3(u, v, w);
}
<<<<<<< Updated upstream
// -------------------------- Collision Solver--------------------------------------
=======

bool IsDuplicated(BufferInt2& primList, BufferInt& indices, int i0, int i1, int i2, int i3)
{
	for (int i = 0; i < primList.GetSize(); ++i)
	{
		int idx = primList.m_Data[i].x;
		int vIdx = indices.m_Data[idx];
		int fIdx1 = indices.m_Data[idx + 1];
		int fIdx2 = indices.m_Data[idx + 2];
		int fIdx3 = indices.m_Data[idx + 3];
		if (vIdx == i0 && fIdx1 == i1 && fIdx2 == i2 && fIdx3 == i3)
		{
			return true;
		}
	}
	return false;
}

// -------------------------- Collision Solver--------------------------------------
bool CollisionSolver::ExtendedVFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3)
{
	if (i0 == i1 || i0 == i2 || i0 == i3)
		return false;
	const float3& x0 = vtxPos, v0 = vtxPrdP - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1Pos - x0, x2 = p2Pos - x0, x3 = p3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (p1PrdP - p1Pos) - v0, v2 = (p2PrdP - p2Pos) - v0, v3 = (p3PrdP - p3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[5];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	t[nsol + 1] = 0;
	std::sort(t, (t + nsol + 2));
	for (int i = 0; i < (nsol + 2); i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
		float3 colliPos0 = pos(vtxPos, vtxPrdP, t[i]), colliPos1 = pos(p1Pos, p1PrdP, t[i]), colliPos2 = pos(p2Pos, p2PrdP, t[i]), colliPos3 = pos(p3Pos, p3PrdP, t[i]);
		contact.n = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
		float3 colliN = normalize(cross(colliPos2 - colliPos1, colliPos3 - colliPos1));
		double dist;
		dist = dot(colliPos0 - colliPos1, colliN);
		float3 weight;
		weight = BarycentricCoord(colliPos0, colliPos1, colliPos2, colliPos3);
		contact.w[0] = 1;
		contact.w[1] = weight.x;
		contact.w[2] = weight.y;
		contact.w[3] = weight.z;
		float edge1 = IO::Distance(colliPos1, colliPos2), edge2 = IO::Distance(colliPos1, colliPos3), edge3 = IO::Distance(colliPos2, colliPos3); // Distance funtion declear many times
		float epsilon = m_thickness / ((edge1 + edge2 + edge3) / 3.0);
		bool inside;
		if (__USE_EPSILON)
		{
			bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
			bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
			bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
			if (inside1 && inside2 && inside3)
			{
				inside = true;
			}
			else
			{
				inside = false;
			}
		}
		else
		{
			inside = (min(weight.x, weight.y, weight.z) >= 1e-6);
		}
		if ((abs(dist) < 2.0f * m_thickness) && inside)
		{
			return true;
		}
	}
	return false;
}

void CollisionSolver::CCD_SH_Extended()
{
	m_ccdSolverTimer->Tick();
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	m_shs->UpdateSH(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto indices = &(m_pbdObj->meshTopol.indices);
	auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto triList = &(m_pbdObj->meshTopol.primList);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	
	// For VF
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
	// For EE
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	int tp0, tp1, tp2, np0, np1, np2;  // verticies of two triangles
	int countEE = 0;

	// for collision debugging
	m_nContact.m_Data.clear();
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection
		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		m_shs->FindNeighbors(neighborList, triId);

		if(ENABLE_VF)
		{
			for (int vtxId = start; vtxId < start + num; ++vtxId)
			{
				i0 = indices->m_Data[vtxId];
				vtxPos = posBuffer->m_Data[i0];
				vtxPrdP = prdPBuffer->m_Data[i0];
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];
					int start = triList->m_Data[nbtriId].x;
					i1 = indices->m_Data[start];
					i2 = indices->m_Data[start + 1];
					i3 = indices->m_Data[start + 2];
					triPos1 = posBuffer->m_Data[i1];
					triPos2 = posBuffer->m_Data[i2];
					triPos3 = posBuffer->m_Data[i3];
					triPrdP1 = prdPBuffer->m_Data[i1];
					triPrdP2 = prdPBuffer->m_Data[i2];
					triPrdP3 = prdPBuffer->m_Data[i3];
					Contact contact;
					if (ExtendedVFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact, i0, i1, i2, i3))
					{
						bool isAlreadyContact = IsDuplicated(contactData.ctxStartNum, contactData.ctxIndices, i0, i1, i2, i3);
						if (!isAlreadyContact)
						{
							contact.type = Contact::VF;
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(i0);
							contactData.ctxIndices.m_Data.push_back(i1);
							contactData.ctxIndices.m_Data.push_back(i2);
							contactData.ctxIndices.m_Data.push_back(i3);
						}
						else
							continue;
					}
				}
			}

		}
		

		// EE: 9 edge-edge tests for each pair of triangles
		// for each edge of the current triangle
		// x0, x1; x1, x2; x2, x0
		if (ENABLE_EE)
		{
			tp0 = indices->m_Data[start];
			tp1 = indices->m_Data[start + 1];
			tp2 = indices->m_Data[start + 2];
			int2 currTriEdges[3];
			getEdgesOfTri(tp0, tp1, tp2, currTriEdges);

			for (int cEdgeIdx = 0; cEdgeIdx < 3; ++cEdgeIdx)
			{
				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				e0 = currTriEdges[cEdgeIdx].x;
				e0Pos = posBuffer->m_Data[e0];
				e0PrdP = prdPBuffer->m_Data[e0];

				e1 = currTriEdges[cEdgeIdx].y;
				e1Pos = posBuffer->m_Data[e1];
				e1PrdP = prdPBuffer->m_Data[e1];
				//printf("neighborlist size %d\n", neighborList.GetSize());
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];

					// only compare curr triangle with triangles whose ids are larger than the current one
					if (nbtriId <= triId)  // actually shouldn't be equal
						continue;

					int nStart = triList->m_Data[nbtriId].x;
					np0 = indices->m_Data[nStart];
					np1 = indices->m_Data[nStart + 1];
					np2 = indices->m_Data[nStart + 2];
					int2 neighborTriEdges[3];
					getEdgesOfTri(np0, np1, np2, neighborTriEdges);

					for (int nEdgeIdx = 0; nEdgeIdx < 3; ++nEdgeIdx)
					{
						e2 = neighborTriEdges[nEdgeIdx].x;
						e2Pos = posBuffer->m_Data[e2];
						e2PrdP = prdPBuffer->m_Data[e2];

						e3 = neighborTriEdges[nEdgeIdx].y;
						e3Pos = posBuffer->m_Data[e3];
						e3PrdP = prdPBuffer->m_Data[e3];

						// call EETest
						Contact contact;
						if (EETest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, e0, e1, e2, e3))
						{
							countEE++;
							//printf("Adding EE contact...\n");
							//printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(e0);
							contactData.ctxIndices.m_Data.push_back(e1);
							contactData.ctxIndices.m_Data.push_back(e2);
							contactData.ctxIndices.m_Data.push_back(e3);
						}
					}
				}
			}
		}
	}
	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}

// ------------- new method -------------
bool CollisionSolver::VFResolveTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3)
{
	float3 colliFreeN = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool colliFreeDir = RelativePos(vtxPos, p1Pos, colliFreeN);
	//bool colliFreeDir = contact.colliFreeDir;
	float3 prdPN = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
	float3 prdPn = cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP);
	if (dot(prdPn, prdPn) <= 1e-6) // ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ë»ï¿½ï¿½ï¿½Ò»ï¿½ï¿½ï¿½ß£ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Î£ï¿½ï¿½ï¿½ï¿½ï¿½
	{
		return false;
	}
	float dist = dot(vtxPrdP - p1PrdP, prdPN * (colliFreeDir ? 1 : -1)) - 2.0 * m_thickness;
	if (dist >=0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void CollisionSolver::VFResolveNew(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3,
	BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, int debug, int iteration)
{
	float3 vp = vtxPrdP, p1p = p1PrdP, p2p = p2PrdP, p3p = p3PrdP;
	//if (debug == 7&& iteration >59)
	//{
	//	//if (i0 == 4608 || i1 == 4608 || i2 == 4608 || i3 == 4608)
	//	//{
	//		printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
	//		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
	//			vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
	//			vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	//	//}
	//}
	float3 crossColliFreeN = cross(p2Pos - p1Pos, p3Pos - p1Pos);
	float3 colliFreeN = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool colliFreeDir = RelativePos(vtxPos, p1Pos, colliFreeN);
	//bool colliFreeDir = contact.colliFreeDir;
	float3 crossPrdPN = cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP);
	float3 prdPN = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
	float3 p1ProjPrdP = p1PrdP + prdPN * (colliFreeDir ? 1 : -1) * 2.0 * m_thickness;
	float depth = dot(p1ProjPrdP - vtxPrdP, (prdPN * (colliFreeDir ? 1 : -1))) * 0.5; // must > 0
	float3 vtxFix = depth * (prdPN * (colliFreeDir ? 1 : -1));
	float3 p1Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	float3 p2Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	float3 p3Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	fixedBuffer.m_Data[i0] += vtxFix;
	fixedBuffer.m_Data[i1] += p1Fix;
	fixedBuffer.m_Data[i2] += p2Fix;
	fixedBuffer.m_Data[i3] += p3Fix;
	vFixedBuffer.m_Data[i0] += vtxFix;
	fFixedBuffer.m_Data[i1] += p1Fix;
	fFixedBuffer.m_Data[i2] += p2Fix;
	fFixedBuffer.m_Data[i3] += p3Fix;
	float3 weight = BarycentricCoord(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	float sw = weight.x * weight.x + weight.y * weight.y + weight.z * weight.z;
	vtxPrdP += vtxFix;
	p1PrdP += p1Fix * weight.x / sw;
	p2PrdP += p2Fix * weight.y / sw;
	p3PrdP += p3Fix * weight.z / sw; 
	if (isnan(vtxPrdP.x) || isnan(vtxPrdP.y) || isnan(vtxPrdP.z) || isnan(p1PrdP.x) || isnan(p1PrdP.y) || isnan(p1PrdP.z) ||
		isnan(p2PrdP.x) || isnan(p2PrdP.y) || isnan(p2PrdP.z) || isnan(p3PrdP.x) || isnan(p3PrdP.y) || isnan(p3PrdP.z))
	{
		printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
		printInfo("weight:", weight);
		printInfo("colliFreeN cross:", crossColliFreeN);
		printInfo("prdpN cross:", crossPrdPN);
		printf("------ depth:%f, sw: %f \n", depth, sw);
		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		vp.x, vp.y, vp.z, p1p.x, p1p.y, p1p.z, p2p.x, p2p.y, p2p.z, p3p.x, p3p.y, p3p.z);
	}
	//if ( debug == 7 && iteration >59)
	//{
		//if (i0 == 676 || i1 == 676 || i2 == 676 || i3 == 676)
		//{
		//	printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
		//	printInfo("weight:", weight);
		//	//printInfo("colliFreeN cross:", crossColliFreeN);
		//	//printInfo("prdpN cross:", crossPrdPN);
		//	printf("------ depth:%f, sw: %f \n", depth, sw);
		//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//		vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
		//}
	//}
}

void CollisionSolver::EEResolve(float3& e0PrdP, float3& e1PrdP, float3& e2PrdP, float3& e3PrdP, float* constrVal, Contact& contact)
{
	float sumW = contact.w[0] * contact.w[0] + contact.w[1] * contact.w[1] + contact.w[2] * contact.w[2] + contact.w[3] * contact.w[3];
	// contact.n should be the prdP Normal now
	float3 k = (-(*constrVal) * contact.n * (contact.colliFreeDir ? 1 : -1)) * (1.0 / sumW);
	e0PrdP += (contact.w[0] * k);
	e1PrdP += (contact.w[1] * k);
	e2PrdP += (contact.w[2] * k);
	e3PrdP += (contact.w[3] * k);

	if (isnan(e0PrdP.x) || isnan(e0PrdP.y) || isnan(e0PrdP.z) || isnan(e1PrdP.x) || isnan(e1PrdP.y) || isnan(e1PrdP.z) ||
		isnan(e2PrdP.x) || isnan(e2PrdP.y) || isnan(e2PrdP.z) || isnan(e3PrdP.x) || isnan(e3PrdP.y) || isnan(e3PrdP.z))
	{
		//printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
		printf("weight: %f %f %f %f", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
		printf("------ sw: %f \n", sumW);
		printInfo("k", k);
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		//vp.x, vp.y, vp.z, p1p.x, p1p.y, p1p.z, p2p.x, p2p.y, p2p.z, p3p.x, p3p.y, p3p.z);
	}
}

void CollisionSolver::CollisionResolveNew(
	BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, 
	int debug, int iteration, int& debugFrameId)
{	
	m_ccdSolverTimer->Tick();
	beforeColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer; // for collision debugging
	m_vfResolveTimes.m_Data.resize(m_pbdObj->meshTopol.posBuffer.GetSize(), 0);

	auto ctxIndices = &(contactData.ctxIndices);
	auto ctxStartNum = &(contactData.ctxStartNum);
	auto ctxList = &(contactData.ctxs);
	auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto m_prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
	{
		Contact contact;
		// For VF
		int start, i0, i1, i2, i3;
		float3 vtxPos, p1Pos, p2Pos, p3Pos;
		//// For EE
		//float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
		//int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges

		contact = ctxList->m_Data[ctxId];
		start = ctxStartNum->m_Data[ctxId].x;
		i0 = ctxIndices->m_Data[start];
		i1 = ctxIndices->m_Data[start + 1];
		i2 = ctxIndices->m_Data[start + 2];
		i3 = ctxIndices->m_Data[start + 3];
		vtxPos = posBuffer->m_Data[i0];
		p1Pos = posBuffer->m_Data[i1];
		p2Pos = posBuffer->m_Data[i2];
		p3Pos = posBuffer->m_Data[i3];
		if (contact.type == Contact::VF)
		{
			if (VFResolveTest(
				vtxPos, p1Pos, p2Pos, p3Pos, 
				m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3], 
				contact, 
				i0, i1, i2, i3))
			{
				VFResolveNew(vtxPos, p1Pos, p2Pos, p3Pos,
					m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
					contact, i0, i1, i2, i3, fixedBuffer, vFixedBuffer, fFixedBuffer, debug, iteration);
			}
			else
				continue;
		}
		//else if (contact.type == Contact::EE)
		//{
		//	start = ctxStartNum->m_Data[ctxId].x;
		//	e0 = ctxIndices->m_Data[start];
		//	e1 = ctxIndices->m_Data[start + 1];
		//	e2 = ctxIndices->m_Data[start + 2];
		//	e3 = ctxIndices->m_Data[start + 3];
		//	e0Pos = posBuffer->m_Data[e0];
		//	e1Pos = posBuffer->m_Data[e1];
		//	e2Pos = posBuffer->m_Data[e2];
		//	e3Pos = posBuffer->m_Data[e3];

		//	
		//	float constrVal = 0.0f;
		//	bool debugThis = false;

		//	if(e0 = 3686 && e1 == 1237 && e2 == 3060 && e3 == 3913)
		//	{
		//		debugThis = true;
		//	}
		//	if (EEResolveTest(
		//		m_prdPBuffer->m_Data[e0],
		//		m_prdPBuffer->m_Data[e1],
		//		m_prdPBuffer->m_Data[e2],
		//		m_prdPBuffer->m_Data[e3],
		//		&constrVal,
		//		contact,
		//		debugThis))
		//	{
		//		printf("This contact should be resolved!\n");
		//		printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
		//		printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
		//		printf("\tconstrVal: %f\n", constrVal);
		//		EEResolve(
		//			m_prdPBuffer->m_Data[e0],
		//			m_prdPBuffer->m_Data[e1],
		//			m_prdPBuffer->m_Data[e2],
		//			m_prdPBuffer->m_Data[e3],
		//			&constrVal,
		//			contact);
		//	}
		//}
		//else
		//{
		//	;   // impossible
		//}
	}
	for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
	{
		Contact contact;
		// For EE
		float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
		int start, e0, e1, e2, e3;  // e0-e1; e2-e3  two edges

		contact = ctxList->m_Data[ctxId];
		start = ctxStartNum->m_Data[ctxId].x;
		
		if (contact.type == Contact::EE)
		{
			start = ctxStartNum->m_Data[ctxId].x;
			e0 = ctxIndices->m_Data[start];
			e1 = ctxIndices->m_Data[start + 1];
			e2 = ctxIndices->m_Data[start + 2];
			e3 = ctxIndices->m_Data[start + 3];
			e0Pos = posBuffer->m_Data[e0];
			e1Pos = posBuffer->m_Data[e1];
			e2Pos = posBuffer->m_Data[e2];
			e3Pos = posBuffer->m_Data[e3];


			float constrVal = 0.0f;
			bool debugThis = true;

			//if (e0 = 3686 && e1 == 1237 && e2 == 3060 && e3 == 3913)
			//{
			//	debugThis = true;
			//}
			if (EEResolveTest(
				m_prdPBuffer->m_Data[e0],
				m_prdPBuffer->m_Data[e1],
				m_prdPBuffer->m_Data[e2],
				m_prdPBuffer->m_Data[e3],
				&constrVal,
				contact))
			{
				printf("This contact should be resolved!\n");
				printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
				printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
				printf("\tconstrVal: %f\n", constrVal);
				EEResolve(
					m_prdPBuffer->m_Data[e0],
					m_prdPBuffer->m_Data[e1],
					m_prdPBuffer->m_Data[e2],
					m_prdPBuffer->m_Data[e3],
					&constrVal,
					contact);
			}
		}
	}
	afterColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
	m_ccdSolverTimer->Tock();
	//PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}

bool ExtendedVFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3)
{
	if (i0 == i1 || i0 == i2 || i0 == i3)
		return false;
	const float3& x0 = vtxPos, v0 = vtxPrdP - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1Pos - x0, x2 = p2Pos - x0, x3 = p3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (p1PrdP - p1Pos) - v0, v2 = (p2PrdP - p2Pos) - v0, v3 = (p3PrdP - p3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[5];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	t[nsol + 1] = 0;
	std::sort(t, (t + nsol + 2));
	for (int i = 0; i < (nsol + 2); i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
		float3 colliPos0 = pos(vtxPos, vtxPrdP, t[i]), colliPos1 = pos(p1Pos, p1PrdP, t[i]), colliPos2 = pos(p2Pos, p2PrdP, t[i]), colliPos3 = pos(p3Pos, p3PrdP, t[i]);
		contact.n = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
		float3 colliN = normalize(cross(colliPos2 - colliPos1, colliPos3 - colliPos1));
		double dist;
		dist = dot(colliPos0 - colliPos1, colliN);
		float3 weight;
		weight = BarycentricCoord(colliPos0, colliPos1, colliPos2, colliPos3);
		contact.w[0] = 1;
		contact.w[1] = weight.x;
		contact.w[2] = weight.y;
		contact.w[3] = weight.z;
		bool inside;
		if (__USE_EPSILON)
		{
			float edge1 = IO::Distance(colliPos1, colliPos2), edge2 = IO::Distance(colliPos1, colliPos3), edge3 = IO::Distance(colliPos2, colliPos3); // Distance funtion declear many times
			float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
			bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
			bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
			bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
			if (inside1 && inside2 && inside3)
			{
				inside = true;
			}
			else
			{
				inside = false;
			}
		}
		else
		{
			inside = (min(weight.x, weight.y, weight.z) >= 1e-6);
		}
		//if (i0 == 264 && (i1 == 62 || i1 == 184 || i1 == 198) && (i2 == 62 || i2 == 184 || i2 == 198) && (i3 == 62 || i3 == 184 || i3 == 198))
	/*	if (i0 == 264 && i1 == 321 && i2 == 108 && i3 == 249)
		{
			printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
			printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
				vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
				contact.t);
			printf("w: (%f, %f, %f)\n", weight.x, weight.y, weight.z);
			cout << "inside: " << inside << endl;
			printf("contact t: %f\n", t[i]);
			printf("d :%f\n", abs(dist) - 2 * thickness);
		}*/
		if ((abs(dist) - 2.0f * thickness < -1e-6) && inside)
		{
			//if (i0 == 264 && i1 == 321 && i2 == 108 && i3 == 249)
			//	printf("--------t%f %f %f-----\n", t[i], dist, 2.0f * thickness);
			return true;
		}
	}
	return false;
}

bool NarrowVFCCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3)
{
	if (i0 == i1 || i0 == i2 || i0 == i3)
		return false;
	const float3& x0 = vtxPos, v0 = vtxPrdP - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1Pos - x0, x2 = p2Pos - x0, x3 = p3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (p1PrdP - p1Pos) - v0, v2 = (p2PrdP - p2Pos) - v0, v3 = (p3PrdP - p3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	for (int i = 0; i < nsol; i++)
	{
		if (t[i] < 0 || t[i] > 1)
		{
			continue;
		}
		contact.t = t[i];
		float3 colliPos0 = pos(vtxPos, vtxPrdP, t[i]), colliPos1 = pos(p1Pos, p1PrdP, t[i]), colliPos2 = pos(p2Pos, p2PrdP, t[i]), colliPos3 = pos(p3Pos, p3PrdP, t[i]);
		double d; // is vtx on tirangle
		bool inside;
		float3 weight;
		contact.n = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
		float3 cn = normalize(cross(colliPos2 - colliPos1, colliPos3 - colliPos1));
		d = dot(colliPos0 - colliPos1, cn);
		weight = BarycentricCoord(colliPos0, colliPos1, colliPos2, colliPos3);
		contact.w[0] = 1;
		contact.w[1] = weight.x;
		contact.w[2] = weight.y;
		contact.w[3] = weight.z;
		if (__USE_EPSILON)
		{
			float edge1 = IO::Distance(p1PrdP, p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
			float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
			bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
			bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
			bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
			if (inside1 && inside2 && inside3)
			{
				inside = true;
			}
			else
			{
				inside = false;
			}
		}
		else
		{
			inside = (min(weight.x, weight.y, weight.z) >= 1e-6);
		}
		if (abs(d) < 1e-6 && inside)
		{
			return true;
		}
	}
	return false;
}
>>>>>>> Stashed changes

bool CollisionSolver::VFTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
	float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
	Contact& contact, int i0, int i1, int i2, int i3)
{
	if (VFCCDTest(vtx_o, p1_o, p2_o, p3_o, vtx_p, p1_p, p2_p, p3_p, Contact::VF, contact, i0, i1, i2, i3))
	{
		//printf("vf contact\n");
		contact.type = Contact::VF;
		return true;
	}
	else if (VFDCDTest(vtx_o, p1_o, p2_o, p3_o, vtx_p, p1_p, p2_p, p3_p, contact))
	{
<<<<<<< Updated upstream
=======
		contact.type = Contact::VF;
		//// for collision debugging
		//vfIndices.m_Data.push_back(i0);
		//vfIndices.m_Data.push_back(i1);
		//vfIndices.m_Data.push_back(i2);
		//vfIndices.m_Data.push_back(i3);
		//
		//printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		//	vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
		//	contact.w[1], contact.w[2], contact.w[3], contact.t);
		//printf("------ dcd ------\n");
		//if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
		//{
		//	printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//		vtx_o.x, vtx_o.y, vtx_o.z, p1_o.x, p1_o.y, p1_o.z, p2_o.x, p1_o.y, p1_o.z, p3_o.x, p3_o.y, p3_o.z,
		//		vtx_p.x, vtx_p.y, vtx_p.z, p1_p.x, p1_p.y, p1_p.z, p2_p.x, p2_p.y, p2_p.z, p3_p.x, p3_p.y, p3_p.z,
		//		contact.w[1], contact.w[2], contact.w[3]);
		//	printf("------print after vf dcd test\n");
		//}
>>>>>>> Stashed changes
		//printf("vf thickness contact\n");
		contact.type = Contact::VFDCD;
		return true;
	}
	else
		return false;
}

bool CollisionSolver::VFDCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact)
{
<<<<<<< Updated upstream
	float d = Point2Plane(vtx_p, p1_p, p2_p, p3_p);
	if (d - 2 * m_thickness >= 1e-6)
		return false;
	contact.w[0] = 1;
	contact.w[1] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).x;
	contact.w[2] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).y;
	contact.w[3] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).z;
	contact.n = normalize(cross(normalize(p2_p - p1_p), normalize(p3_p - p1_p)));
	bool inside;
	inside = (min(contact.w[1], contact.w[2], contact.w[3]) >= 1e-6);
=======
	float d = Point2Plane(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	if (d - 2 * m_thickness >= -1e-6)
		return false;
	float3 weight = BarycentricCoord(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	contact.w[0] = 1;
	contact.w[1] = weight.x;
	contact.w[2] = weight.y;
	contact.w[3] = weight.z;
	contact.n = normalize(cross(normalize(p2PrdP - p1PrdP), normalize(p3PrdP - p1PrdP))); 
	contact.t = 0.0f;
	bool inside;
	float edge1 = IO::Distance(p1PrdP, p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
	float epsilon = m_thickness / ((edge1 + edge2 + edge3) / 3.0);
	bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
	bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
	bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
	if (inside1 && inside2 && inside3)
	{
		inside = true;
	}
	else
	{
		inside = false;
	}
	//inside = (min(contact.w[1], contact.w[2], contact.w[3]) >= 1e-6);
>>>>>>> Stashed changes
	if (!inside)
		return false;
	return true;
}

bool CollisionSolver::VFCCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact::Type type, Contact& contact, int i0, int i1, int i2, int i3)
{
	const float3& x0 = vtxPos, v0 = vtxPrdP - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1Pos - x0, x2 = p2Pos - x0, x3 = p3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (p1PrdP - p1Pos) - v0, v2 = (p2PrdP - p2Pos) - v0, v3 = (p3PrdP - p3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	for (int i = 0; i < nsol; i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
<<<<<<< Updated upstream
		float3 x0 = pos(vtx_o, vtx_p, t[i]), x1 = pos(p1_o, p1_p, t[i]), x2 = pos(p2_o, p2_p, t[i]), x3 = pos(p3_o, p3_p, t[i]);
=======
		float3 colliPos0 = pos(vtxPos, vtxPrdP, t[i]), colliPos1 = pos(p1Pos, p1PrdP, t[i]), colliPos2 = pos(p2Pos, p2PrdP, t[i]), colliPos3 = pos(p3Pos, p3PrdP, t[i]);
>>>>>>> Stashed changes
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
<<<<<<< Updated upstream
		contact.n = normalize(cross(p2_p - p1_p, p3_p - p1_p));
=======
		float3 weight;
		contact.n = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
>>>>>>> Stashed changes
		//compute weight and normal
		if (type == Contact::VF)
		{
			float3 cn = normalize(cross(x2 - x1, x3 - x1));
			d = dot(x0 - x1, cn);
			float3 weight = BarycentricCoord(x0, x1, x2, x3);
			contact.w[0] = 1;
			contact.w[1] = weight.x;
			contact.w[2] = weight.y;
			contact.w[3] = weight.z;
			inside = (min(w[1], w[2], w[3]) >= 1e-6);
			/*
			if (i0 == 497 )
			{
				cout << a;
				cout << " " << b << endl;
				printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					vtx_o.x, vtx_o.y, vtx_o.z, p1_o.x, p1_o.y, p1_o.z, p2_o.x, p2_o.y, p2_o.z, p3_o.x, p3_o.y, p3_o.z,
					vtx_p.x, vtx_p.y, vtx_p.z, p1_p.x, p1_p.y, p1_p.z, p2_p.x, p2_p.y, p2_p.z, p3_p.x, p3_p.y, p3_p.z);
				printf("-------(%d, %d, %d)-------\n", i1, i2, i3);
				printf("x0: (%f, %f, %f)\n", x0.x, x0.y, x0.z);
				printf("x1: (%f, %f, %f)\n", x1.x, x1.y, x1.z);
				printf("x2: (%f, %f, %f)\n", x2.x, x2.y, x2.z);
				printf("x3: (%f, %f, %f)\n", x3.x, x3.y, x3.z);
				printf("w: (%f, %f, %f)\n", w[1], w[2], w[3]);
				printf("normal: (%f, %f, %f)\n", n.x, n.y, n.z);
				printf("normal: (%f, %f, %f)\n", contact.n.x, contact.n.y, contact.n.z);
				cout << "inside: " << inside << endl;
				printf("contact t: %f\n", contact.t);
				printf("d :%f\n", abs(d));
			}
			*/
		}
		if (abs(d) < 1e-6 && inside)
		{
			return true;
		}
	}
	return false;
}

void CollisionSolver::VFResolve(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
	Contact contact, int i0, int i1, int i2, int i3)
{
	float3 fn = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool posRelativePos = RelativePos(vtxPos, p1Pos, fn);
	//printf("normal: (%f, %f, %f)\n", contact.n.x, contact.n.y, contact.n.z);
	//cout << "sameSide:" << sameSide << endl;
	//printf("------344 vtxPrd: (%f, %f, %f)--------\n", m_prdPBuffer.m_Data[344].x, m_prdPBuffer.m_Data[344].y, m_prdPBuffer.m_Data[344].z);
	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
	//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z, 
	//	vtxPrd.x, vtxPrd.y, vtxPrd.z, p1Prd.x, p1Prd.y, p1Prd.z, p2Prd.x, p2Prd.y, p2Prd.z, p3Prd.x, p3Prd.y, p3Prd.z);
	//printf("vtxPos: (%f, %f, %f)\n", vtxPos.x, vtxPos.y, vtxPos.z);
	//printf("p1Pos: (%f, %f, %f)\n", p1Pos.x, p1Pos.y, p1Pos.z);
	//printf("p2Pos: (%f, %f, %f)\n", p2Pos.x, p2Pos.y, p2Pos.z);
	//printf("p3Pos: (%f, %f, %f))\n", p3Pos.x, p3Pos.y, p3Pos.z);
	//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);

	//float3 n = contact.n;
	float3 n = normalize(cross(p2Prd - p1Prd, p3Prd - p1Prd));
	float depth = Point2Plane(vtxPrd, p1Prd, p2Prd, p3Prd);
	bool prdPRelativePos = RelativePos(vtxPrd, p1Prd, n);
<<<<<<< Updated upstream
	float dp;
	if ((posRelativePos == prdPRelativePos) && (depth < 2 * m_thickness)) // vf dcd
=======
	float dp = 0.0f;
	float3 weight = BarycentricCoord(vtxPrd, p1Prd, p2Prd, p3Prd);
	if (posRelativePos == prdPRelativePos && depth < 2 * m_thickness) // vf dcd
>>>>>>> Stashed changes
	{
		dp = (2 * m_thickness - depth) * 0.5;
		if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
			printf("---------vf dcd----------\n");
	}
	else if (posRelativePos != prdPRelativePos) // vf ccd
	{
		dp = (depth + 2 * m_thickness) * 0.5;
		if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
			printf("---------vf ccd----------\n");
	}
	// for collision debugging
	m_resolveDepths[i0].m_Data.push_back(depth);
	m_resolveDepths[i1].m_Data.push_back(depth);
	m_resolveDepths[i2].m_Data.push_back(depth);
	m_resolveDepths[i3].m_Data.push_back(depth);
	m_resolveTimes.m_Data[i0].y += 1;
	m_resolveTimes.m_Data[i1].y += 1;
	m_resolveTimes.m_Data[i2].y += 1;
	m_resolveTimes.m_Data[i3].y += 1;

	float sw = contact.w[1] * contact.w[1] + contact.w[2] * contact.w[2] + contact.w[3] * contact.w[3];

	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n", 
	//	vtxPrd.x, vtxPrd.y, vtxPrd.z, p1Prd.x, p1Prd.y, p1Prd.z, p2Prd.x, p2Prd.y, p2Prd.z, p3Prd.x, p3Prd.y, p3Prd.z, 
	//	depth, contact.w[1], contact.w[2], contact.w[3], Contact::VFDCD);

	if (posRelativePos)
	{
		vtxPrd += n * dp;
		p1Prd += -n * dp * contact.w[1] / sw;
		p2Prd += -n * dp * contact.w[2] / sw;
		p3Prd += -n * dp * contact.w[3] / sw;

	}
	else
	{
		vtxPrd += -n * dp;
		p1Prd += n * dp * contact.w[1] / sw;
		p2Prd += n * dp * contact.w[2] / sw;
		p3Prd += n * dp * contact.w[3] / sw;
	}
	//printf("normal: (%f, %f, %f) \n", n.x, n.y, n.z);
	//printf("dp: %f\n", dp);
	//printf("w1: (%f) ", (contact.w[1] / sw));
	//printf("w2: (%f) ", (contact.w[2] / sw));
	//printf("w3: (%f)\n", (contact.w[3] / sw));
	//printf("change1: (%f) ", (dp * contact.w[1] / sw));
	//printf("change2: (%f) ", (dp * contact.w[2] / sw));
	//printf("change3: (%f)\n", (dp * contact.w[3] / sw));
	//p1Prd += -n * dp * contact.w[1];
	//p2Prd += -n * dp * contact.w[2];
	//p3Prd += -n * dp * contact.w[3];
	//auto temp = n * dp;
	//m_resolveDepths[i0].m_Data.push_back((n * dp));

	//printf("normal: (%f, %f, %f) \n", n.x, n.y, n.z); 
	//printf("point id: %d ", i0);
	//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	//printf("point id: %d ", i1);
	//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	//printf("point id: %d ", i2);
	//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	//printf("point id: %d ", i3);
	//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
	//printf("depth: %f\n", depth);
	//printf("dp: %f\n", dp);
	//printf("------344 vtxPrd: (%f, %f, %f)--------\n", m_prdPBuffer.m_Data[344].x, m_prdPBuffer.m_Data[344].y, m_prdPBuffer.m_Data[344].z);
	//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
}

// 1 : resolve iteration ï¿½ï¿½ï¿½ï¿½
// 2 : resolve times
// 2 : condition ï¿½ï¿½
void CollisionSolver::CollisionResolve()
{
	m_ccdSolverTimer->Tick();
	beforeColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer; // for collision debugging
	for (int i = 0; i < 2; ++i)
	{
		//m_shs->UpdateSH(m_pbdObj->constrPBDBuffer.prdPBuffer);
		Contact contact;
		int start, i0, i1, i2, i3;
		float3 vtxPos, p1Pos, p2Pos, p3Pos;
		auto ctxIndices = &(contactData.ctxIndices);
		auto ctxStartNum = &(contactData.ctxStartNum);
		auto ctxList = &(contactData.ctxs);
		auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
		auto m_prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
		for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
		{
			contact = ctxList->m_Data[ctxId];
			start = ctxStartNum->m_Data[ctxId].x;
			i0 = ctxIndices->m_Data[start];
			i1 = ctxIndices->m_Data[start + 1];
			i2 = ctxIndices->m_Data[start + 2];
			i3 = ctxIndices->m_Data[start + 3];
			vtxPos = posBuffer->m_Data[i0];
			p1Pos = posBuffer->m_Data[i1];
			p2Pos = posBuffer->m_Data[i2];
			p3Pos = posBuffer->m_Data[i3];
			if (contact.type == Contact::VF)
			{
				// ï¿½ï¿½ï¿½ï¿½ï¿½Ç°Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½collision-freeï¿½ï¿½Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½Í¬ï¿½ï¿½ï¿½Ò¾ï¿½ï¿½ï¿½ï¿½ï¿½Ú¶ï¿½ï¿½ï¿½thicknessï¿½ï¿½ ï¿½Í²ï¿½ï¿½ï¿½ï¿½ï¿½
				float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
				bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
				float3 nn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
				bool cRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nn);
				///*
				if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
				{
					auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					if (cRelatedPos == fRelatedPos && distance >= 2 * m_thickness)
						continue;
					else
					{
						printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
						printf("distance: %f, 2*thickness: %f\n", distance, 2 * m_thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
						printf("------------------------------------------\n");

						VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
							m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
							contact, i0, i1, i2, i3);

						auto newDistance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
						float3 nnn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
						bool ncRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nnn);
						printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * m_thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;
						m_debugFrameID++;

						string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(m_debugFrameID) + ".cache";
						m_pbdObj->constrPBDBuffer.prdPBuffer.SetName("P");
						Topology temp;
						temp.indices = m_pbdObj->meshTopol.indices;
						temp.primList = m_pbdObj->meshTopol.primList;
						temp.posBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
						temp.indices.SetName("Indices");
						temp.primList.SetName("primList");
						temp.posBuffer.SetName("P");
						IO::SaveToplogy(temp, path);
						printf("------------------------------------------m_debugFrameID %dsave-------------------------------------\n", m_debugFrameID);
						//printf("vtxPrd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z);
						//printf("p1Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z);
						//printf("p2Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z);
						//printf("p3Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
					}
				}
				else
				{
					auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					if ((cRelatedPos == fRelatedPos) && (distance > 2 * m_thickness))
						continue;
					else
					{
						VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
							m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
							contact, i0, i1, i2, i3);
					}
				}
				//*/

				/*
				auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
				if (cRelatedPos == fRelatedPos && distance > 2 * m_thickness)
					continue;
				else
				{
					VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
						m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
						contact, i0, i1, i2, i3);
				}
				*/

			}
			//if (Contact::EE)
			//{
			//	EEResolve()
			//}
		}
	}
	afterColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}

void CollisionSolver::CollisionResolveTest()
{
	Contact contact;
	int start, i0, i1, i2, i3;
	float3 vtxPos, p1Pos, p2Pos, p3Pos;
	auto ctxIndices = &(contactData.ctxIndices);
	auto ctxStartNum = &(contactData.ctxStartNum);
	auto ctxList = &(contactData.ctxs);
	auto posBuffer = &(m_topol.posBuffer);
	for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
	{
		contact = ctxList->m_Data[ctxId];
		start = ctxStartNum->m_Data[ctxId].x;
		i0 = ctxIndices->m_Data[start];
		i1 = ctxIndices->m_Data[start + 1];
		i2 = ctxIndices->m_Data[start + 2];
		i3 = ctxIndices->m_Data[start + 3];
		vtxPos = posBuffer->m_Data[i0];
		p1Pos = posBuffer->m_Data[i1];
		p2Pos = posBuffer->m_Data[i2];
		p3Pos = posBuffer->m_Data[i3];
<<<<<<< Updated upstream

		if (Contact::VF || Contact::VFDCD)
=======
		if (Contact::VF)
>>>>>>> Stashed changes
		{
			// ï¿½ï¿½ï¿½ï¿½ï¿½Ç°Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½collision-freeï¿½ï¿½Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½Í¬ï¿½ï¿½ï¿½Ò¾ï¿½ï¿½ï¿½ï¿½ï¿½Ú¶ï¿½ï¿½ï¿½thicknessï¿½ï¿½ ï¿½Í²ï¿½ï¿½ï¿½ï¿½ï¿½
			float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
			bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
			//float3 nnew = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
			bool cRelatedPos = RelativePos(m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], contact.n);
			if (cRelatedPos == fRelatedPos &&
				(Point2Plane(m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], m_prdPBuffer.m_Data[i2], m_prdPBuffer.m_Data[i3]) > 2 * m_thickness))
				continue;
			else
			{
				VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
					m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], m_prdPBuffer.m_Data[i2], m_prdPBuffer.m_Data[i3],
					contact, i0, i1, i2, i3);
			}

			//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
			//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
			//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
			//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
		}
<<<<<<< Updated upstream
=======
	}
	afterColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}


void CollisionSolver::CollisionResolveTest()
{
	Contact contact;
	int start, i0, i1, i2, i3;
	float3 vtxPos, p1Pos, p2Pos, p3Pos;
	auto ctxIndices = &(contactData.ctxIndices);
	auto ctxStartNum = &(contactData.ctxStartNum);
	auto ctxList = &(contactData.ctxs);
	auto posBuffer = &(m_topol.posBuffer);
	for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
	{
		contact = ctxList->m_Data[ctxId];
		start = ctxStartNum->m_Data[ctxId].x;
		i0 = ctxIndices->m_Data[start];
		i1 = ctxIndices->m_Data[start + 1];
		i2 = ctxIndices->m_Data[start + 2];
		i3 = ctxIndices->m_Data[start + 3];
		vtxPos = posBuffer->m_Data[i0];
		p1Pos = posBuffer->m_Data[i1];
		p2Pos = posBuffer->m_Data[i2];
		p3Pos = posBuffer->m_Data[i3];

		if (contact.type == Contact::VF)
		{
			// ï¿½ï¿½ï¿½ï¿½ï¿½Ç°Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½collision-freeï¿½ï¿½Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½Í¬ï¿½ï¿½ï¿½Ò¾ï¿½ï¿½ï¿½ï¿½ï¿½Ú¶ï¿½ï¿½ï¿½thicknessï¿½ï¿½ ï¿½Í²ï¿½ï¿½ï¿½ï¿½ï¿½
			float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
			bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
			//float3 nnew = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
			bool cRelatedPos = RelativePos(m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], contact.n);
			if (cRelatedPos == fRelatedPos &&
				(Point2Plane(m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], m_prdPBuffer.m_Data[i2], m_prdPBuffer.m_Data[i3]) > 2 * m_thickness))
				continue;
			else
			{
				VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
					m_prdPBuffer.m_Data[i0], m_prdPBuffer.m_Data[i1], m_prdPBuffer.m_Data[i2], m_prdPBuffer.m_Data[i3],
					contact, i0, i1, i2, i3);
			}

			//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
			//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
			//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
			//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
		}
>>>>>>> Stashed changes
		//if (Contact::EE)
		//{
		//	EEResolve()
		//}
		//ColliWithShpGrd();
	}

}

void CollisionSolver::CCD_N2()
{
	cout << __FUNCTION__ << endl;
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	auto indices = &(m_pbdObj->meshTopol.indices);
	auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto triList = &(m_pbdObj->meshTopol.primList);
	auto m_prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		for (int vtxId = start; vtxId < start + num; ++vtxId)
		{
			i0 = indices->m_Data[vtxId];
			vtxPos = posBuffer->m_Data[i0];
			vtxPrdP = m_prdPBuffer->m_Data[i0];
			for (int nbtriId = 0; nbtriId < triList->GetSize(); ++nbtriId)
			{
				int start = triList->m_Data[nbtriId].x;
				i1 = indices->m_Data[start];
				i2 = indices->m_Data[start + 1];
				i3 = indices->m_Data[start + 2];
				triPos1 = posBuffer->m_Data[i1];
				triPos2 = posBuffer->m_Data[i2];
				triPos3 = posBuffer->m_Data[i3];
				triPrdP1 = m_prdPBuffer->m_Data[i1];
				triPrdP2 = m_prdPBuffer->m_Data[i2];
				triPrdP3 = m_prdPBuffer->m_Data[i3];
				if (nbtriId == triId)
					continue;
				Contact contact;


				if (VFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact, i0, i1, i2, i3))
				{
					contactData.ctxs.m_Data.push_back(contact);
					contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
					contactData.ctxIndices.m_Data.push_back(i0);
					contactData.ctxIndices.m_Data.push_back(i1);
					contactData.ctxIndices.m_Data.push_back(i2);
					contactData.ctxIndices.m_Data.push_back(i3);
				}
			}
		}
	}
}

void CollisionSolver::CCD_N2Test()
{
	cout << __FUNCTION__ << endl;
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	auto indices = &(m_topol.indices);
	auto posBuffer = &(m_topol.posBuffer);
	auto triList = &(m_topol.primList);
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		for (int vtxId = start; vtxId < start + num; ++vtxId)
		{
			i0 = indices->m_Data[vtxId];
			vtxPos = posBuffer->m_Data[i0];
			vtxPrdP = m_prdPBuffer.m_Data[i0];
			for (int nbtriId = 0; nbtriId < triList->GetSize(); ++nbtriId)
			{
				int start = triList->m_Data[nbtriId].x;
				i1 = indices->m_Data[start];
				i2 = indices->m_Data[start + 1];
				i3 = indices->m_Data[start + 2];
				triPos1 = posBuffer->m_Data[i1];
				triPos2 = posBuffer->m_Data[i2];
				triPos3 = posBuffer->m_Data[i3];
				triPrdP1 = m_prdPBuffer.m_Data[i1];
				triPrdP2 = m_prdPBuffer.m_Data[i2];
				triPrdP3 = m_prdPBuffer.m_Data[i3];
				if (nbtriId == triId)
					continue;
				Contact contact;
				if (VFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact, i0, i1, i2, i3))
				{
					contactData.ctxs.m_Data.push_back(contact);
					contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
					contactData.ctxIndices.m_Data.push_back(i0);
					contactData.ctxIndices.m_Data.push_back(i1);
					contactData.ctxIndices.m_Data.push_back(i2);
					contactData.ctxIndices.m_Data.push_back(i3);
				}
			}
		}
	}
}

static std::set<int> maskSet = { 17,21,22,23,26,34,41,47,52,55,57,58,87,92,103,106,107,134,136,137,152,153,162,164,168,170,171,174,192,214,217,228,233,248,251,252,259,274,280,285,307,308,310,317,318,319,320,321,323,373,384,411,412,414,415,419,420,421,422,424,425,432,447,448,454,455,456,466,478,479,488,494,525,541,542,544,596,599,601,606,609,610,611,631,643,644,646,650,658,662,664,666,683,728,729,730,732,742,753,758,759,761,763,787,809,827,836,869,871,874,875,877,878,885,887,902,904,905,908,909,913,914,915,919,921,929,930,936,937,938,939,979,1028,1115,1118,1119,1121,1127,1161,1169,1170,1171,1172,1174,1177,1187,1258,1297,1298,1299,1300,1301,1302,1303,1308,1309,1311,1313,1317,1318,1319,1335,1338,1340,1341,1342,1343,1345,1346,1347,1348,1349,1350,1351,1352,1354,1358,1359,1361,1362,1363,1365,1367,1370,1372,1373,1374,1393,1400,1401,1453,1468,1471,1475,1477,1479,1482,1483,1484,1485,1487,1510,1511,1512,1514,1515,1516,1518,1519,1520,1522,1556,1557,1560,1569,1600,1601,1606,1609,1611,1612,1613,1614,1632,1641,1644,1645,1648,1668,1670,1671,1689,1691,1692,1702,1707,1763,1787,1788,1789,1790,1791,1796,1798,1799,1813,1828,1835,1836,1837,1839,1842,1845,1846,1848,1849,1850,1851,1852,1853,1859,1861,1864,1896,1917,1976,1977,1978,1979,1981,1992,2014,2016,2017,2019,2021,2024,2030,2041,2043,2073,2076,2078,2089,2091,2092,2098,2099,2103,2106,2167,2168,2170,2291,2292,2294,2297,2303,2304,2330,2331,2333,2334,2335,2337,2338,2339,2340,2344,2345,2348,2359,2384,2386,2388,2425,2426,2428,2442,2445,2446,2448,2449,2450,2460,2461,2462,2464,2467,2468,2469,2483,2538,2562,2563,2566,2638,2639,2640,2642,2648,2649,2651,2653,2655,2662,2717,2718,2719,2731,2738,2752,2759,2760,2762,2771,2773,2800,2820,2823,2917,2919,2920,2921,2922,2928,2929,2930,2931,2937,2939,2940,2941,2942,2943,2944,2946,2947,2949,2950,2951,2952,2953,2955,2958,2978,2983,2986,2987,2988,2989,2990,2991,2996,2997,2998,3071,3072,3075,3078,3084,3085,3086,3087,3090,3097,3098,3099,3100,3101,3102,3103,3106,3118,3120,3122,3124,3125,3126,3127,3128,3129,3130,3131,3148,3149,3155,3157,3158,3159,3162,3165,3166,3188,3197,3199,3200,3201,3204,3205,3229,3233,3235,3236,3237,3238,3239,3240,3241,3242,3243,3261,3263,3264,3284,3286,3287,3303,3304,3305,3307,3308,3310,3331,3354,3364,3372,3387,3398,3407,3412,3428,3429,3430,3478,3504,3505,3508,3509,3511,3512,3516,3518,3519,3525,3526,3527,3548,3549,3555,3558,3560,3585,3586,3595,3641,3664,3665,3668,3706,3707,3721,3732,3752,3769,3775,3782,3783,3793,3795,3803,3805,3806,3807,3823,3824,3834,3836,3837,3838,3839,3840,3842,3845,3846,3848,3873,4056,4057,4058,4085,4086,4092,4093,4094,4097,4098,4099,4100,4101,4102,4104,4106,4127,4128,4129,4144,4221,4222,4223,4238,4239,4240,4244,4245,4246,4247,4248,4254,4255,4256,4492,4493,4494,4626,4627,4628,4738,4740,4747,4753,4758,4760,4770,4772,4780,4781,4827,4828,4831,4835,4838,4842,4844,4846,4847,4850,4851,4852,4859,4861,4863,4866,4867,4873,4875,4876,4878,4879,4922,4925,4928,4929,4931,4936,4938,4940,4942,4944,4963,5013,5022,5026,5030,5031,5047,5059,5063,5072,5073,5074,5077,5084,5090,5096,5098,5135,5136,5139,5141,5179,5181,5184,5186,5187,5190,5192,5199,5200,5212,5240,5247,5254,5259,5295,5305,5307,5308,5309,5312,5315,5317,5322,5323,5353,5355,5362,5363,5366,5371,5374,5406,5410,5411,5421,5466,5467,5472,5478,5489,5490,5499,5500,5503,5504,5505,5506,5507,5509,5515,5516,5517,5518,5519,5523,5524,5528,5531,5532,5533,5534,5536,5537,5538,5539,5540,5542,5547,5548,5551,5577,5578,5579,5584,5589,5590,5591,5598,5600,5603,5604,5606,5609,5638,5665,5680,5682,5683,5691,5740,5748,5750,5756,5765,5769,5770,5771,5776,5780,5781,5784,5786,5790,5798,5799,5801,5802,5819,5821,5823,5824,5844,5851,5854,5855,5856,5860,5861,5865,5867,5869,5871,5876,5877,5880,5882,5884,5888,5892,5896,5905,5908,5913,5917,5920,5987,6005,6028,6040,6042,6043,6049,6057,6058,6061,6066,6081,6082,6083,6114,6117,6124,6127,6128,6129,6130,6134,6143,6144,6162,6163,6166,6169,6172,6207,6213,6216,6251,6258,6261,6263,6264,6265,6270,6283,6285,6289,6300,6301,6304,6307,6318,6326,6327,6328,6330,6333,6335,6344,6350,6355,6361,6362,6376,6377,6411,6415,6434,6456,6457,6459,6461,6467,6470,6485,6495,6502,6504,6509,6511,6515,6517,6520,6521,6524,6525,6527,6528,6530,6535,6536,6539,6540,6542,6552,6574,6586,6598,6606,6661,6662,6663,6664,6665,6671,6673,6675,6677,6682,6699,6704,6706,6707,6708,6710,6713,6714,6715,6720,6724,6727,6729,6731,6732,6736,6738,6739,6740,6797,6801,6819,6828,6830,6838,6839,6841,6842,6843,6854,6855,6856,6857,6881,6907,6910,6918,6923,6925,6936,6939,6940,6941,6945,6951,6952,6990,7002,7003,7004,7008,7009,7017,7024,7025,7028,7029,7030,7032,7034,7035,7036,7073,7101,7108,7111,7124,7125,7126,7127,7128,7131,7132,7137,7148,7151,7165,7183,7184,7185,7191,7192,7202,7203,7208,7220,7221,7225,7229,7257,7264,7266,7270,7273,7276,7277,7278,7281,7284,7295,7305,7313,7319,7321,7322,7326,7351,7352,7353,7369,7378,7389,7392,7401,7416,7417,7418,7422,7425,7438,7441,7457,7459,7466,7472,7477,7486,7492,7498,7504,7510,7511,7513,7515,7519 };

void CollisionSolver::CCD_SH()
{
	m_ccdSolverTimer->Tick();
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	m_shs->UpdateSH(m_pbdObj->constrPBDBuffer.prdPBuffer);
	auto indices = &(m_pbdObj->meshTopol.indices);
	auto posBuffer = &(m_pbdObj->meshTopol.posBuffer);
	auto triList = &(m_pbdObj->meshTopol.primList);
	auto prdPBuffer = &(m_pbdObj->constrPBDBuffer.prdPBuffer);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	// For VF
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
<<<<<<< Updated upstream
	BufferInt vtxList;
=======
	// For EE
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	int tp0, tp1, tp2, np0, np1, np2;  // verticies of two triangles

>>>>>>> Stashed changes
	// for collision debugging
	m_nContact.m_Data.clear();
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection

		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		m_shs->FindNeighbors(neighborList, triId);
		//m_shs->FindNeighbors(neighborList, triId);
<<<<<<< Updated upstream
		/*printf("neighborlist size %d\n", neighborList.GetSize());
		for (int i = 0; i < neighborList.GetSize(); ++i)
			printf("\t %d\n", neighborList.m_Data[i]);*/
		for (int vtxId = start; vtxId < start + num; ++vtxId)
=======

		if (ENABLE_VF)
>>>>>>> Stashed changes
		{
			for (int vtxId = start; vtxId < start + num; ++vtxId)
			{
<<<<<<< Updated upstream
				int nbtriId = neighborList.m_Data[nbIdx];
				
				// For Mark Debug
				//printf("checking nbtriId: %d\n", nbtriId);
				auto mySetIterator = maskSet.find(nbtriId);
				int newNbId = std::distance(maskSet.begin(), mySetIterator);
				//printf("checking start: %d\n", triList->m_Data[newNbId].x);
				int start = triList->m_Data[newNbId].x;

				//int nbtriId = neighborList.m_Data[nbIdx];
				//int start = triList->m_Data[nbtriId].x;
				i1 = indices->m_Data[start];
				i2 = indices->m_Data[start + 1];
				i3 = indices->m_Data[start + 2];
				triPos1 = posBuffer->m_Data[i1];
				triPos2 = posBuffer->m_Data[i2];
				triPos3 = posBuffer->m_Data[i3];
				triPrdP1 = prdPBuffer->m_Data[i1];
				triPrdP2 = prdPBuffer->m_Data[i2];
				triPrdP3 = prdPBuffer->m_Data[i3];
				Contact contact;
=======
				i0 = indices->m_Data[vtxId];
				vtxPos = posBuffer->m_Data[i0];
				vtxPrdP = prdPBuffer->m_Data[i0];
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];
					int start = triList->m_Data[nbtriId].x;
					i1 = indices->m_Data[start];
					i2 = indices->m_Data[start + 1];
					i3 = indices->m_Data[start + 2];
					triPos1 = posBuffer->m_Data[i1];
					triPos2 = posBuffer->m_Data[i2];
					triPos3 = posBuffer->m_Data[i3];
					triPrdP1 = prdPBuffer->m_Data[i1];
					triPrdP2 = prdPBuffer->m_Data[i2];
					triPrdP3 = prdPBuffer->m_Data[i3];
					Contact contact;
					if (VFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact, i0, i1, i2, i3))
					{
						/*
						// for collision debugging
						m_nContact.m_Data[i0] += 1;
						m_nContact.m_Data[i1] += 1;
						m_nContact.m_Data[i2] += 1;
						m_nContact.m_Data[i3] += 1;
						*/
						contactData.ctxs.m_Data.push_back(contact);
						contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
						contactData.ctxIndices.m_Data.push_back(i0);
						contactData.ctxIndices.m_Data.push_back(i1);
						contactData.ctxIndices.m_Data.push_back(i2);
						contactData.ctxIndices.m_Data.push_back(i3);
					}
				}
			}
		}
		
>>>>>>> Stashed changes

		// EE: 9 edge-edge tests for each pair of triangles
		// for each edge of the current triangle
		// x0, x1; x1, x2; x2, x0
		if (ENABLE_EE)
		{

			tp0 = indices->m_Data[start];
			tp1 = indices->m_Data[start + 1];
			tp2 = indices->m_Data[start + 2];
			int2 currTriEdges[3];
			getEdgesOfTri(tp0, tp1, tp2, currTriEdges);

			for (int cEdgeIdx = 0; cEdgeIdx < 3; ++cEdgeIdx)
			{
				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				e0 = currTriEdges[cEdgeIdx].x;
				e0Pos = posBuffer->m_Data[e0];
				e0PrdP = prdPBuffer->m_Data[e0];

				e1 = currTriEdges[cEdgeIdx].y;
				e1Pos = posBuffer->m_Data[e1];
				e1PrdP = prdPBuffer->m_Data[e1];

				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
<<<<<<< Updated upstream
					/*
					bool isDetect = false;
					for (int i = 0; i < vtxList.GetSize(); ++i)
					{
						if (i0 == vtxList.m_Data[i])
						{
							isDetect = true;
							break;
						}
					}
					if (isDetect)
						continue;
					vtxList.m_Data.push_back(i0);
					*/
					/*
					// for collision debugging
					m_nContact.m_Data[i0] += 1;
					m_nContact.m_Data[i1] += 1;
					m_nContact.m_Data[i2] += 1;
					m_nContact.m_Data[i3] += 1;
					*/
					contactData.ctxs.m_Data.push_back(contact);
					contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
					contactData.ctxIndices.m_Data.push_back(i0);
					contactData.ctxIndices.m_Data.push_back(i1);
					contactData.ctxIndices.m_Data.push_back(i2);
					contactData.ctxIndices.m_Data.push_back(i3);
=======
					int nbtriId = neighborList.m_Data[nbIdx];

					// only compare curr triangle with triangles whose ids are larger than the current one
					if (nbtriId <= triId)  // actually shouldn't be equal
						continue;

					int nStart = triList->m_Data[nbtriId].x;
					np0 = indices->m_Data[nStart];
					np1 = indices->m_Data[nStart + 1];
					np2 = indices->m_Data[nStart + 2];
					int2 neighborTriEdges[3];
					getEdgesOfTri(np0, np1, np2, neighborTriEdges);

					for (int nEdgeIdx = 0; nEdgeIdx < 3; ++nEdgeIdx)
					{
						e2 = neighborTriEdges[nEdgeIdx].x;
						e2Pos = posBuffer->m_Data[e2];
						e2PrdP = prdPBuffer->m_Data[e2];

						e3 = neighborTriEdges[nEdgeIdx].y;
						e3Pos = posBuffer->m_Data[e3];
						e3PrdP = prdPBuffer->m_Data[e3];

						// call EETest
						Contact contact;
						if (EETest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, e0, e1, e2, e3))
						{
							//printf("Adding EE contact...\n");
							//printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(e0);
							contactData.ctxIndices.m_Data.push_back(e1);
							contactData.ctxIndices.m_Data.push_back(e2);
							contactData.ctxIndices.m_Data.push_back(e3);
						}
					}
>>>>>>> Stashed changes
				}
			}
		}
	}
	//std::cout << "Contact size is " << contactData.ctxs.GetSize() << std::endl;
	//for (int i = 0; i < contactData.ctxs.GetSize(); ++i)
	//{
	//	printf("edge0: %d-%d; edge1: %d-%d\n", 
	//		contactData.ctxIndices.m_Data[i*4 + 0],
	//		contactData.ctxIndices.m_Data[i*4 + 1],
	//		contactData.ctxIndices.m_Data[i*4 + 2],
	//		contactData.ctxIndices.m_Data[i*4 + 3]);
	//}

	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}

void CollisionSolver::getEdgesOfTri(int tp0, int tp1, int tp2, int2* rtnEdges)
{
	rtnEdges[0] = make_int2(tp0, tp1);
	rtnEdges[1] = make_int2(tp1, tp2);
	rtnEdges[2] = make_int2(tp2, tp0);
}

double CollisionSolver::signed_ee_distance(const float3& x0, const float3& x1,
	const float3& y0, const float3& y1,
	float3* n, double* w)
{
	// Does NOT consider parallel lines
	float3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	if (norm2(*n) < 1e-6) // parallel
		return infinity;  // not return infinity
	*n = normalize(*n);
	double h = dot(x0 - y0, *n);
	double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
		b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	w[0] = a0 / (a0 + a1);
	w[1] = a1 / (a0 + a1);
	w[2] = -b0 / (b0 + b1);
	w[3] = -b1 / (b0 + b1);
	// printf("w[0]: %f; w[1]: %f; w[2]: %f; w[3]: %f\n", w[0], w[1], w[2], w[3]);
	return h;

	// Deal with parallel lines
	//float3 _n; if (!n) n = &_n;
	//double _w[4]; if (!w) w = _w;
	//*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	//if (norm2(*n) < 1e-8) {
	//	// special case: parallel lines
	//	float3 e0 = normalize(x1 - x0), e1 = normalize(y1 - y0);

	//	double p0min = dot(x0, e0), p0max = dot(x1, e0), p1min = dot(y0, e0), p1max = dot(y1, e0);
	//	if (p1max < p1min) swap(p1max, p1min);

	//	double a = max(p0min, p1min), b = min(p0max, p1max), c = 0.5*(a + b);
	//	if (a > b) return infinity;

	//	float3 d = (y0 - x0) - dot(y0 - x0, e0)*e0;

	//	if (n) *n = normalize(-d);
	//	if (w) {
	//		w[1] = (c - dot(x0, e0)) / norm(x1 - x0);
	//		w[0] = 1.0 - w[1];
	//		w[3] = -(dot(e0, e1)*c - dot(y0, e1)) / norm(y1 - y0);
	//		w[2] = -1.0 - w[3];
	//	}
	//	return norm(d);
	//}
	//*n = normalize(*n);
	//double h = dot(x0 - y0, *n);
	//double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
	//	b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	//w[0] = a0 / (a0 + a1);
	//w[1] = a1 / (a0 + a1);
	//w[2] = -b0 / (b0 + b1);
	//w[3] = -b1 / (b0 + b1);
	//return h;
}

// EE Detection Test
bool CollisionSolver::EETest(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos,
	float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP,
	Contact& contact, int e0, int e1, int e2, int e3)
{
	if (e0 == e2 || e0 == e3 || e1 == e2 || e1 == e3)
		return false;
	if (EECCDTest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, e0, e1, e2, e3))
	{
		contact.type = Contact::EE;
		return true;
	}
	//else if (EEDCDTest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, e0, e1, e2, e3))
	//{
	//	printf("EE thickness contact\n");
	//	contact.type = Contact::EEDCD;
	//	return true;
	//}
	else
		return false;
}

bool CollisionSolver::EECCDTest(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos,
	float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP,
	Contact& contact, int e0, int e1, int e2, int e3)
{
	const float3& x0 = e0Pos, v0 = e0PrdP - x0;
	float3 x1 = e1Pos - x0, x2 = e2Pos - x0, x3 = e3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (e1PrdP - e1Pos) - v0, v2 = (e2PrdP - e2Pos) - v0, v3 = (e3PrdP - e3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3),
		a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3),
		a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3),
		a3 = stp(v1, v2, v3);

	if (abs(a0) < 1e-6 * norm(x1) * norm(x2) * norm(x3))
		return false; // initially coplanar

	double t[5];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	t[nsol + 1] = 0; // also check at start of timestep [prevent thickness DCD issue]
	std::sort(t, (t + nsol + 2));
	for (int i = 0; i < (nsol + 2); i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
		float3 colliPos0 = pos(e0Pos, e0PrdP, t[i]), colliPos1 = pos(e1Pos, e1PrdP, t[i]), colliPos2 = pos(e2Pos, e2PrdP, t[i]), colliPos3 = pos(e3Pos, e3PrdP, t[i]);
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
		d = signed_ee_distance(colliPos0, colliPos1, colliPos2, colliPos3, &n, w);
		if (__USE_EPSILON)
		{
			float edge1 = IO::Distance(e0PrdP, e1PrdP), edge2 = IO::Distance(e2PrdP, e3PrdP); // Distance funtion declear many times
			float epsilon = m_thickness / ((edge1 + edge2) / 2.0f);
			bool inside1 = (min(w[0], w[1]) >= -epsilon) && (max(w[0], w[1]) <= 1 + epsilon);
			bool inside2 = (min(-w[2], -w[3]) >= -epsilon) && (max(-w[2], -w[3]) <= 1 + epsilon);

			inside = inside1 && inside2;
		}
		else
			inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);

		// printf("EE: d is %f; inside is %d\n", d, inside);
		// TODO: check this again for EE resolve
		if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > -1e-6)  // calculate here for ee resolve
			n = -n;
		if (abs(d) < (2.0 * m_thickness) && inside)   // add the thickness for preventing EE DCD
		{
			float c = calcRestEEConstr(e0Pos, e1Pos, e2Pos, e3Pos, contact);

			//if (e0 == 74 && e1 == 36 && e2 == 62 && e3 == 40)
			//	printf("dubug ee info: d is %f\n", abs(d));
			// assert(c != 0.0);
			if (c == 0.0)
			{
				// printf("d and inside but constr is 0: %d, %d, %d, %d\n", e0, e1, e2, e3);
				return false;
			}

			//printf("Adding EE: %d,%d,%d,%d ----- c: %f\n", e0, e1,e2,e3, c);
			initContactFlag(c, contact);  // initialize contact colli free dir before return
			return true;
		}
	}
	return false;
}

//bool CollisionSolver::EEDCDTest(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos,
//	float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP,
//	Contact& contact, int e0, int e1, int e2, int e3)
//{
//	// TODO: fill in the DCD
//	return true;
//}

float CollisionSolver::calcRestEEConstr(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos, Contact& contact)
{
	float3 currNormal = cross(normalize(e1Pos - e0Pos), normalize(e3Pos - e2Pos));
	if (norm2(currNormal) <= 1e-6)
	{
		//printf("calcRestEEConstr get parallel\n");
		return 0.0f;
	}
	currNormal = normalize(currNormal);
	contact.n = currNormal;
	//printf("normal is %f %f %f\n", currNormal.x, currNormal.y, currNormal.z);
	double a0 = stp(e3Pos - e1Pos, e2Pos - e1Pos, currNormal), a1 = stp(e2Pos - e0Pos, e3Pos - e0Pos, currNormal),
		b0 = stp(e0Pos - e3Pos, e1Pos - e3Pos, currNormal), b1 = stp(e1Pos - e2Pos, e0Pos - e2Pos, currNormal);
	contact.w[0] = a0 / (a0 + a1);
	contact.w[1] = a1 / (a0 + a1);
	contact.w[2] = -b0 / (b0 + b1);
	contact.w[3] = -b1 / (b0 + b1);
	float3 dist = (contact.w[0] * e0Pos + contact.w[1] * e1Pos) - ((-contact.w[2]) * e2Pos + (-contact.w[3]) * e3Pos);
	//printf("dist is %f %f %f\n", dist.x, dist.y, dist.z);
	float c = dot(currNormal, dist);
	return c;
}

float CollisionSolver::calcCurrEEConstr(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, Contact& contact)
{
	float3 currNormal = cross(normalize(e1PrdP - e0PrdP), normalize(e3PrdP - e2PrdP));
	if (norm2(currNormal) <= 1e-6)
	{
		//printf("calcCurrEEConstr get parallel\n");
		return 0.0f;
	}
	currNormal = normalize(currNormal);
	contact.n = currNormal;
	//printf("normal is %f %f %f\n", currNormal.x, currNormal.y, currNormal.z);
	double a0 = stp(e3PrdP - e1PrdP, e2PrdP - e1PrdP, currNormal), a1 = stp(e2PrdP - e0PrdP, e3PrdP - e0PrdP, currNormal),
		b0 = stp(e0PrdP - e3PrdP, e1PrdP - e3PrdP, currNormal), b1 = stp(e1PrdP - e2PrdP, e0PrdP - e2PrdP, currNormal);
	contact.w[0] = a0 / (a0 + a1);
	contact.w[1] = a1 / (a0 + a1);
	contact.w[2] = -b0 / (b0 + b1);
	contact.w[3] = -b1 / (b0 + b1);

	float3 dist = (contact.w[0] * e0PrdP + contact.w[1] * e1PrdP) - ((-contact.w[2]) * e2PrdP + (-contact.w[3]) * e3PrdP);
	//printf("dist is %f %f %f\n", dist.x, dist.y, dist.z);
	float c = dot(currNormal * (contact.colliFreeDir ? 1 : -1), dist) - 2.0*m_thickness;

	//printf("\tEEResolveTest\n");
	//printf("colliFreeDir is %d\n", contact.colliFreeDir);
	//printf("\tprdPN: (%f, %f, %f)\n", contact.n.x, contact.n.y, contact.n.z);
	//printf("\te0PrdP: (%f, %f, %f); e1PrdP: (%f, %f, %f); e2PrdP: (%f, %f, %f); e3PrdP: (%f, %f, %f)\n",
	//	e0PrdP.x, e0PrdP.y, e0PrdP.z,
	//	e1PrdP.x, e1PrdP.y, e1PrdP.z,
	//	e2PrdP.x, e2PrdP.y, e2PrdP.z,
	//	e3PrdP.x, e3PrdP.y, e3PrdP.z);
	//printf("\ta0: %f; a1: %f; b0: %f; b1: %f\n", a0, a1, b0, b1);
	//printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);


	return c;
}

void CollisionSolver::initContactFlag(float constrVal, Contact& contact)
{
	contact.colliFreeDir = constrVal > 0.0 ? true : false;
}

bool CollisionSolver::EEResolveTest(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, float* constrVal, Contact& contact)
{
	*constrVal = calcCurrEEConstr(e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact);
	// printf("constrVal2: %f\n ", *constrVal);

	if ((*constrVal) >= 0.0f)
	{
		return false;
	}
	else
		return true;
}

void CCDTestMain()
{
	CCDTest ccdTest;
	ccdTest.PrepareTestData("D://0310ContinuousCollisionDectection//ccdTestData//3ClothStaticTestRestP.txt",
		"D://0310ContinuousCollisionDectection//ccdTestData//3ClothStaticTestPrdP.txt");
	//ccdTest.PrepareTestData("D://0310ContinuousCollisionDectection//ccdTestData//CCDResolvePos.txt",
	//										"D://0310ContinuousCollisionDectection//ccdTestData//CCDResolvePrdP.txt");
	/**/
	CollisionSolver colliSolver;

	//float3 cellSize = make_float3(3.0f, 3.0f, 3.0f);
	//float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	//uint3 gridSize = make_uint3(5, 1, 5);

	//// initialize SH
	//SpatialHashSystem shs(ccdTest.topol.posBuffer, ccdTest.topol.indices, CPU);
	//shs.SetGridCenter(gridCenter);
	//shs.SetGridSize(gridSize);
	//shs.SetDivision(cellSize);
	//shs.InitSH();
	//shs.UpdateSH(0.0f);

	colliSolver.SetTargetTest(ccdTest.topol, ccdTest.prdPBuffer);
	colliSolver.SetThickness(0.03f);

	//auto start = chrono::steady_clock::now();
	//clock_t tStart = clock();

	for (int i = 0; i < 4; ++i)
	{
		//// collision detection
		//colliSolver.CCD_SH();
		colliSolver.CCD_N2Test();
		//// collision resolves
		colliSolver.CollisionResolveTest();
	}

	//auto end = chrono::steady_clock::now();
	//auto diff = end - start;
	//cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;

	////printf("contact number: %d", colliSolver.contactData.ctxStartNum.GetSize());
	//set<int> idList;
	//for (int i = 0; i < colliSolver.contactData.ctxStartNum.GetSize(); ++i)
	//{
	//	int id = colliSolver.contactData.ctxStartNum.m_Data[i].x;
	//	idList.insert(colliSolver.contactData.ctxIndices.m_Data[id]);
	//}
	//set<int>::iterator it;
	//printf("id: ");
	//for (it = idList.begin(); it != idList.end(); it++)
	//{
	//	printf("%d -- ", *it);
	//}

	colliSolver.SaveResult();
}

void CollisionSolver::SaveResult()
{
	//printf("------229 vtxPrd: (%f, %f, %f)--------\n", m_prdPBuffer.m_Data[229].x, m_prdPBuffer.m_Data[229].y, m_prdPBuffer.m_Data[229].z);
	m_topol.indices.SetName("Indices");
	m_topol.posBuffer.SetName("P");
	m_prdPBuffer.SetName("P");
	m_topol.primList.SetName("primList");
	m_topol.posBuffer = m_prdPBuffer;
	IO::SaveToplogy(m_topol, "D://0310ContinuousCollisionDectection//ccdTestData//3clothCCDTestTopol.cache");
	IO::SaveBuffer(m_prdPBuffer, "D://0310ContinuousCollisionDectection//ccdTestData//prdP.cache");
	cout << "topol saved" << endl;
}

void CollisionSolver::SaveCollision(string path)
{
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;
	beforeColliPrdPBuffer.SetName("bPrdP");
	IO::SaveBuffer(beforeColliPrdPBuffer, ofs);
	afterColliPrdPBuffer.SetName("aPrdP");
	IO::SaveBuffer(afterColliPrdPBuffer, ofs);
	auto ctxPrim = &(contactData.ctxStartNum);
	auto ctxIndicies = &(contactData.ctxIndices);
	//m_nContact.m_Data.resize(m_pbdObj->meshTopol.posBuffer.GetSize(), 0);
	//for (int i = 0; i < ctxPrim->GetSize(); ++i)
	//{
	//	int i0 = ctxPrim->m_Data[i].x;
	//	int idx = ctxIndicies->m_Data[i0];
	//	m_nContact.m_Data[idx]+=1;
	//}

	m_nContact.SetName("nContact");
	IO::SaveBuffer(m_nContact, ofs);
	//contactData.Save(ofs);
}

<<<<<<< Updated upstream
//----------------------- CCDTest ----------------------------
=======
void CollisionSolver::SaveContact(string path)
{
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;
	printf("Saving Constact...\n");
	printf("contactVFMap size %d, contactVHitMap size %d\n", m_contactVFMap.size(), m_contactVHitMap.size());
	assert(m_contactVFMap.size() == m_contactVHitMap.size());
	map<int, std::set<int>>::iterator it;
	map<int, BufferVector3f>::iterator it2;
	for (it = m_contactVFMap.begin(), it2 = m_contactVHitMap.begin();
		it != m_contactVFMap.end() && it2 != m_contactVHitMap.end();
		it++, it2++)
	{
		auto vertexId = it->first;
		auto triIds = &(m_contactVFMap[vertexId]);
		auto hitPos = &(m_contactVHitMap[vertexId]);
		ofs << vertexId << ":";
		int i = 0;
		for (std::set<int>::iterator itSet = triIds->begin();
			itSet != triIds->end() && i < hitPos->GetSize();
			++itSet, ++i)
		{
			float3 currHitPos = hitPos->m_Data[i];
			ofs << (*itSet) << ",";
			ofs << currHitPos.x << "," << currHitPos.y << "," << currHitPos.z << ",";
		}
		ofs << ";";
	}
	ofs << std::endl;

	/*for (it2 = m_contactVHitMap.begin(); it2 != m_contactVHitMap.end(); it2++)
	{
		auto vertexId = it2->first;
		auto hitPos = &(m_contactVHitMap[vertexId]);
		ofs << vertexId << ":";
		for (int i = 0; i < hitPos->GetSize(); ++i)
		{
			float3 currHitPos = hitPos->m_Data[i];
			printf("%f %f %f ", currHitPos.x, currHitPos.y, currHitPos.z);
			ofs << currHitPos.x << "," << currHitPos.y << "," << currHitPos.z;
		}
		ofs << ";";
	}
	ofs << std::endl;*/
}

void CollisionSolver::SavePrdPBeforeCCD(string path)
{
	Topology temp;
	temp.indices = m_pbdObj->meshTopol.indices;
	temp.primList = m_pbdObj->meshTopol.primList;
	temp.posBuffer = beforeColliPrdPBuffer;
	temp.indices.SetName("Indices");
	temp.primList.SetName("primList");
	temp.posBuffer.SetName("P");
	IO::SaveToplogy(temp, path);
}

void CollisionSolver::SavePrdPAfterCCD(string path)
{
	Topology temp;
	temp.indices = m_pbdObj->meshTopol.indices;
	temp.primList = m_pbdObj->meshTopol.primList;
	temp.posBuffer = afterColliPrdPBuffer;
	temp.indices.SetName("Indices");
	temp.primList.SetName("primList");
	temp.posBuffer.SetName("P");
	IO::SaveToplogy(temp, path);
}
>>>>>>> Stashed changes

//----------------------- CCDTest ----------------------------
void CCDTest::PrepareTestData(string topolPath, string prdPath)
{
	readMeshFromTxt(topolPath);
	readBufferFromTxt(prdPath);
	createPrimList(3);
	cout << topol.indices.GetSize() << endl;
	cout << topol.primList.GetSize() << endl;
}

void CCDTest::readMeshFromTxt(string filename)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = splitString(buffer, firstSplit)[1];
		auto dataStringList = splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				// std::cout << dataStringList[i] << std::endl;
				float3 vertexPos = make_float3(std::stof(dataStringList[i]),
					std::stof(dataStringList[i + 1]),
					std::stof(dataStringList[i + 2]));
				// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				topol.posBuffer.m_Data.push_back(vertexPos);
			}
			std::cout << "topol.posBuffer: " << topol.posBuffer.GetSize() << endl;
			assert(topol.posBuffer.GetSize() == dataStringList.size() / 3);
		}
		else  // second line of the input file: vertices tet
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i]));
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i + 1]));
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i + 2]));
			}
		}
		++lineNum;
	}
	in.close();
	printf("Read File Done!\n");
}

void CCDTest::readBufferFromTxt(string filename)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = splitString(buffer, firstSplit)[1];
		auto dataStringList = splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				// std::cout << dataStringList[i] << std::endl;
				float3 vertexPos = make_float3(std::stof(dataStringList[i]),
					std::stof(dataStringList[i + 1]),
					std::stof(dataStringList[i + 2]));
				// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				prdPBuffer.m_Data.push_back(vertexPos);
			}
			std::cout << "topol.posBuffer: " << prdPBuffer.GetSize() << endl;
			assert(prdPBuffer.GetSize() == dataStringList.size() / 3);
		}
	}
	in.close();
	printf("Read File Done!\n");
}

void CCDTest::splitString(const std::string& src, std::vector<std::string>& v, const std::string& split)
{
	std::string::size_type pos1, pos2;
	pos2 = src.find(split);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(src.substr(pos1, pos2 - pos1));
		pos1 = pos2 + split.size();
		pos2 = src.find(split, pos1);
	}
	if (pos1 != src.length())
		v.push_back(src.substr(pos1));
}

std::vector<std::string> CCDTest::splitString(const std::string& src, const std::string& split)
{
	std::vector<std::string> _ret = std::vector<std::string>();
	splitString(src, _ret, split);
	return _ret;
}

void CCDTest::createPrimList(int count)
{
	for (int i = 0; i < topol.indices.GetSize() / count; ++i)
	{
		int2 p;
		p.x = i * count;
		p.y = count;
		topol.primList.m_Data.push_back(p);
	}
}
// -------------------------- contactData ---------------------------
void ContactData::Save(std::ofstream& ofs)
{
	ctxIndices.SetName("Indices");
	IO::SaveBuffer(ctxIndices, ofs);
	ctxStartNum.SetName("primList");
	IO::SaveBuffer(ctxStartNum, ofs);
}

<<<<<<< Updated upstream

// DoVF(dt)
/*{
	prdBuffer = posBuffer + t * v;    //calculate four points position after time t
	½âÈý´Î·½³Ì£¬È¡×îÐ¡·Ç¸ºÊµ¸ù
	float e = thickness/edgeAverage;    // calculate epsilon
	if( t < dt )
		½«t´øÈë·½³ÌÇó½âw1, w2
		if( -e < w1 < 1+e && -e < w2 < 1+e)
			w3 = 1 - w1 - w2
			if( -e < w1 < 1+e)
				colliPrimList.push_back(make_int2(colliIndices.GetSize(), 4))
				colliIndices.push_back(µ±Ç°µãidx£¬ÒÔ¼°µ±Ç°collideÈý½ÇÐÎÈý¸öµãidx)
				colliType.push_back(VF);
}*/

// DoEE(tolerance)
/*{
	¼ÆËãÁ½Ìõ±ßÏòÁ¿
	if(Á½Ìõ±ßÆ½ÐÐ)
		continue;
	½â·½³ÌÇóaºÍb£¬¼ÆËãÁ½¸ö½»µã
	if(a<1 && b < 1) //Á½¸öµã¶¼ÔÚÏß¶ÎÉÏ
		d = Distance(a1, a2)
		if( d < tolerance )
			colliPrimList.push_back(make_int2(colliIndices.GetSize(), 4))
			colliIndices.push_back(ËÄ¸öµãµÄidx)
			colliType.push_back(EE);
	else if()
}*/
=======
// ---------------- Data Oriented ----------------
// ------------- new method ------------
void CCD_SH_Extended(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol, // ï¿½ï¿½Ò»Ö¡ï¿½ï¿½collision freeï¿½ï¿½ï¿½ï¿½ï¿½ï¿½meshtopol
	BufferVector3f prdPBuffer,
	float thickness)
{
	PBD_DEBUG;
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	shs.UpdateSH(prdPBuffer);
	auto indices = &(meshTopol.indices);
	auto posBuffer = &(meshTopol.posBuffer);
	auto triList = &(meshTopol.primList);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	// For VF
	float3 vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP;
	int i0, i1, i2, i3;
	// For EE
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	int tp0, tp1, tp2, np0, np1, np2;  // verticies of two triangles

	int countEE = 0;
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection
		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		shs.FindNeighbors(neighborList, triId);

		if (ENABLE_VF)
		{
			for (int vtxId = start; vtxId < start + num; ++vtxId)
			{
				i0 = indices->m_Data[vtxId];
				vtxPos = posBuffer->m_Data[i0];
				vtxPrdP = prdPBuffer.m_Data[i0];
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];
					int start = triList->m_Data[nbtriId].x;
					i1 = indices->m_Data[start];
					i2 = indices->m_Data[start + 1];
					i3 = indices->m_Data[start + 2];
					p1Pos = posBuffer->m_Data[i1];
					p2Pos = posBuffer->m_Data[i2];
					p3Pos = posBuffer->m_Data[i3];
					p1PrdP = prdPBuffer.m_Data[i1];
					p2PrdP = prdPBuffer.m_Data[i2];
					p3PrdP = prdPBuffer.m_Data[i3];
					Contact contact;
					if (ExtendedVFTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, contact, thickness, i0, i1, i2, i3))
					{
						bool isAlreadyContact = IsDuplicated(contactData.ctxStartNum, contactData.ctxIndices, i0, i1, i2, i3);
						if (!isAlreadyContact)
						{
							contact.type = Contact::VF;
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(i0);
							contactData.ctxIndices.m_Data.push_back(i1);
							contactData.ctxIndices.m_Data.push_back(i2);
							contactData.ctxIndices.m_Data.push_back(i3);
						}
						else
							continue;
					}
				}
			}
		}

		// EE: 9 edge-edge tests for each pair of triangles
		// for each edge of the current triangle
		// x0, x1; x1, x2; x2, x0
		if (ENABLE_EE)
		{
			tp0 = indices->m_Data[start];
			tp1 = indices->m_Data[start + 1];
			tp2 = indices->m_Data[start + 2];
			int2 currTriEdges[3];
			getEdgesOfTri(tp0, tp1, tp2, currTriEdges);

			for (int cEdgeIdx = 0; cEdgeIdx < 3; ++cEdgeIdx)
			{
				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				e0 = currTriEdges[cEdgeIdx].x;
				e0Pos = posBuffer->m_Data[e0];
				e0PrdP = prdPBuffer.m_Data[e0];

				e1 = currTriEdges[cEdgeIdx].y;
				e1Pos = posBuffer->m_Data[e1];
				e1PrdP = prdPBuffer.m_Data[e1];

				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];

					// only compare curr triangle with triangles whose ids are larger than the current one
					if (nbtriId <= triId)  // actually shouldn't be equal
						continue;

					int nStart = triList->m_Data[nbtriId].x;
					np0 = indices->m_Data[nStart];
					np1 = indices->m_Data[nStart + 1];
					np2 = indices->m_Data[nStart + 2];
					int2 neighborTriEdges[3];
					getEdgesOfTri(np0, np1, np2, neighborTriEdges);

					for (int nEdgeIdx = 0; nEdgeIdx < 3; ++nEdgeIdx)
					{
						e2 = neighborTriEdges[nEdgeIdx].x;
						e2Pos = posBuffer->m_Data[e2];
						e2PrdP = prdPBuffer.m_Data[e2];

						e3 = neighborTriEdges[nEdgeIdx].y;
						e3Pos = posBuffer->m_Data[e3];
						e3PrdP = prdPBuffer.m_Data[e3];

						// call EETest
						Contact contact;
						if (EETest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness, e0, e1, e2, e3))
						{
							countEE++;
							//printf("Adding EE contact...\n");
							//printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(e0);
							contactData.ctxIndices.m_Data.push_back(e1);
							contactData.ctxIndices.m_Data.push_back(e2);
							contactData.ctxIndices.m_Data.push_back(e3);
						}
					}
				}
			}	
		}
	}
	printf("Find EE: %d\n", countEE);
}

void CCD_SH_Narrow(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol, // ï¿½ï¿½Ò»Ö¡ï¿½ï¿½collision freeï¿½ï¿½ï¿½ï¿½ï¿½ï¿½meshtopol
	BufferVector3f prdPBuffer,
	float thickness)
{
	PBD_DEBUG;
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	shs.UpdateSH(prdPBuffer);
	auto indices = &(meshTopol.indices);
	auto posBuffer = &(meshTopol.posBuffer);
	auto triList = &(meshTopol.primList);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	// For VF
	float3 vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP;
	int i0, i1, i2, i3;
	// For EE
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	int tp0, tp1, tp2, np0, np1, np2;  // verticies of two triangles

	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection
		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		shs.FindNeighbors(neighborList, triId);

		if (ENABLE_VF)
		{
			for (int vtxId = start; vtxId < start + num; ++vtxId)
			{
				i0 = indices->m_Data[vtxId];
				vtxPos = posBuffer->m_Data[i0];
				vtxPrdP = prdPBuffer.m_Data[i0];
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];
					int start = triList->m_Data[nbtriId].x;
					i1 = indices->m_Data[start];
					i2 = indices->m_Data[start + 1];
					i3 = indices->m_Data[start + 2];
					p1Pos = posBuffer->m_Data[i1];
					p2Pos = posBuffer->m_Data[i2];
					p3Pos = posBuffer->m_Data[i3];
					p1PrdP = prdPBuffer.m_Data[i1];
					p2PrdP = prdPBuffer.m_Data[i2];
					p3PrdP = prdPBuffer.m_Data[i3];
					Contact contact;
					if (NarrowVFCCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, contact, thickness, i0, i1, i2, i3))
					{
						bool isAlreadyContact = IsDuplicated(contactData.ctxStartNum, contactData.ctxIndices, i0, i1, i2, i3);
						if (!isAlreadyContact)
						{
							contact.type = Contact::VF;
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(i0);
							contactData.ctxIndices.m_Data.push_back(i1);
							contactData.ctxIndices.m_Data.push_back(i2);
							contactData.ctxIndices.m_Data.push_back(i3);
						}
						else
							continue;
					}
				}
			}
		}

		// EE: 9 edge-edge tests for each pair of triangles
		// for each edge of the current triangle
		// x0, x1; x1, x2; x2, x0
		if (ENABLE_EE)
		{
			tp0 = indices->m_Data[start];
			tp1 = indices->m_Data[start + 1];
			tp2 = indices->m_Data[start + 2];
			int2 currTriEdges[3];
			getEdgesOfTri(tp0, tp1, tp2, currTriEdges);

			for (int cEdgeIdx = 0; cEdgeIdx < 3; ++cEdgeIdx)
			{
				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				e0 = currTriEdges[cEdgeIdx].x;
				e0Pos = posBuffer->m_Data[e0];
				e0PrdP = prdPBuffer.m_Data[e0];

				e1 = currTriEdges[cEdgeIdx].y;
				e1Pos = posBuffer->m_Data[e1];
				e1PrdP = prdPBuffer.m_Data[e1];

				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];

					// only compare curr triangle with triangles whose ids are larger than the current one
					if (nbtriId <= triId)  // actually shouldn't be equal
						continue;

					int nStart = triList->m_Data[nbtriId].x;
					np0 = indices->m_Data[nStart];
					np1 = indices->m_Data[nStart + 1];
					np2 = indices->m_Data[nStart + 2];
					int2 neighborTriEdges[3];
					getEdgesOfTri(np0, np1, np2, neighborTriEdges);

					for (int nEdgeIdx = 0; nEdgeIdx < 3; ++nEdgeIdx)
					{
						e2 = neighborTriEdges[nEdgeIdx].x;
						e2Pos = posBuffer->m_Data[e2];
						e2PrdP = prdPBuffer.m_Data[e2];

						e3 = neighborTriEdges[nEdgeIdx].y;
						e3Pos = posBuffer->m_Data[e3];
						e3PrdP = prdPBuffer.m_Data[e3];

						// call EETest
						Contact contact;
						if (EETest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness, e0, e1, e2, e3))
						{
							//printf("Adding EE contact...\n");
							//printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(e0);
							contactData.ctxIndices.m_Data.push_back(e1);
							contactData.ctxIndices.m_Data.push_back(e2);
							contactData.ctxIndices.m_Data.push_back(e3);
						}
					}
				}
			}
		}
	}
}

void CCD_SH(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol, // ï¿½ï¿½Ò»Ö¡ï¿½ï¿½collision freeï¿½ï¿½ï¿½ï¿½ï¿½ï¿½meshtopol
	BufferVector3f prdPBuffer,
	float thickness) // ï¿½ï¿½Ç°Ö¡ï¿½ï¿½ï¿½ï¿½prdp
{
	PBD_DEBUG;
	contactData.ctxs.m_Data.clear();
	contactData.ctxStartNum.m_Data.clear();
	contactData.ctxIndices.m_Data.clear();
	shs.UpdateSH(prdPBuffer); 
	auto indices = &(meshTopol.indices);
	auto posBuffer = &(meshTopol.posBuffer);
	auto triList = &(meshTopol.primList);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	// For VF
	float3 vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP;
	int i0, i1, i2, i3;
	// For EE
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	int tp0, tp1, tp2, np0, np1, np2;  // verticies of two triangles

	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection
		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		shs.FindNeighbors(neighborList, triId);
		if (ENABLE_VF)
		{
			for (int vtxId = start; vtxId < start + num; ++vtxId)
			{
				i0 = indices->m_Data[vtxId];
				vtxPos = posBuffer->m_Data[i0];
				vtxPrdP = prdPBuffer.m_Data[i0];
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];
					int start = triList->m_Data[nbtriId].x;
					i1 = indices->m_Data[start];
					i2 = indices->m_Data[start + 1];
					i3 = indices->m_Data[start + 2];
					p1Pos = posBuffer->m_Data[i1];
					p2Pos = posBuffer->m_Data[i2];
					p3Pos = posBuffer->m_Data[i3];
					p1PrdP = prdPBuffer.m_Data[i1];
					p2PrdP = prdPBuffer.m_Data[i2];
					p3PrdP = prdPBuffer.m_Data[i3];
					Contact contact;

					if (VFTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, contact, thickness, i0, i1, i2, i3))
					{
						contactData.ctxs.m_Data.push_back(contact);
						contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
						contactData.ctxIndices.m_Data.push_back(i0);
						contactData.ctxIndices.m_Data.push_back(i1);
						contactData.ctxIndices.m_Data.push_back(i2);
						contactData.ctxIndices.m_Data.push_back(i3);
					}
				}
			}
		}
		
		// EE: 9 edge-edge tests for each pair of triangles
		// for each edge of the current triangle
		// x0, x1; x1, x2; x2, x0
		if (ENABLE_EE)
		{
			tp0 = indices->m_Data[start];
			tp1 = indices->m_Data[start + 1];
			tp2 = indices->m_Data[start + 2];
			int2 currTriEdges[3];
			getEdgesOfTri(tp0, tp1, tp2, currTriEdges);

			for (int cEdgeIdx = 0; cEdgeIdx < 3; ++cEdgeIdx)
			{
				// for each edge of the neighbors; do the EE Test on edges of its neighbor triangles
				e0 = currTriEdges[cEdgeIdx].x;
				e0Pos = posBuffer->m_Data[e0];
				e0PrdP = prdPBuffer.m_Data[e0];

				e1 = currTriEdges[cEdgeIdx].y;
				e1Pos = posBuffer->m_Data[e1];
				e1PrdP = prdPBuffer.m_Data[e1];
				//printf("neighborlist size %d\n", neighborList.GetSize());
				for (int nbIdx = 0; nbIdx < neighborList.GetSize(); ++nbIdx)
				{
					int nbtriId = neighborList.m_Data[nbIdx];

					// only compare curr triangle with triangles whose ids are larger than the current one
					if (nbtriId <= triId)  // actually shouldn't be equal
						continue;

					int nStart = triList->m_Data[nbtriId].x;
					np0 = indices->m_Data[nStart];
					np1 = indices->m_Data[nStart + 1];
					np2 = indices->m_Data[nStart + 2];
					int2 neighborTriEdges[3];
					getEdgesOfTri(np0, np1, np2, neighborTriEdges);

					for (int nEdgeIdx = 0; nEdgeIdx < 3; ++nEdgeIdx)
					{
						e2 = neighborTriEdges[nEdgeIdx].x;
						e2Pos = posBuffer->m_Data[e2];
						e2PrdP = prdPBuffer.m_Data[e2];

						e3 = neighborTriEdges[nEdgeIdx].y;
						e3Pos = posBuffer->m_Data[e3];
						e3PrdP = prdPBuffer.m_Data[e3];

						// call EETest
						Contact contact;
						//printf("Call EE Test\n");
						//printf("e0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
						if (EETest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness, e0, e1, e2, e3))
						{
							//printf("Adding EE contact...\n");
							printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
							contactData.ctxs.m_Data.push_back(contact);
							contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
							contactData.ctxIndices.m_Data.push_back(e0);
							contactData.ctxIndices.m_Data.push_back(e1);
							contactData.ctxIndices.m_Data.push_back(e2);
							contactData.ctxIndices.m_Data.push_back(e3);
						}
					}
				}
			}
		}
	}
	//std::cout << "Contact size is " << contactData.ctxs.GetSize() << std::endl;
	//for (int i = 0; i < contactData.ctxs.GetSize(); ++i)
	//{
	//	printf("edge0: %d-%d; edge1: %d-%d\n", 
	//		contactData.ctxIndices.m_Data[i*4 + 0],
	//		contactData.ctxIndices.m_Data[i*4 + 1],
	//		contactData.ctxIndices.m_Data[i*4 + 2],
	//		contactData.ctxIndices.m_Data[i*4 + 3]);
	//}
}

bool VFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3)
{
	if (VFCCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, Contact::VF,  contact, thickness, i0, i1, i2, i3))
	{
		contact.type = Contact::VF;
		//// for collision debugging
		//vfIndices.m_Data.push_back(i0);
		//vfIndices.m_Data.push_back(i1);
		//vfIndices.m_Data.push_back(i2);
		//vfIndices.m_Data.push_back(i3);
		//printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		//	vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
		//	contact.w[1], contact.w[2], contact.w[3], contact.t);
		//printf("------ ccd ------\n");
		//if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
		//{
		//	printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//	printf("------print after vf ccd test\n");
		//	
		//	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	//	contact.t,
		//	//	vtxPos.x, vtxPos.y, vtxPos.z, triPos1.x, triPos1.y, triPos1.z, triPos2.x, triPos2.y, triPos2.z, triPos3.x, triPos3.y, triPos3.z,
		//	//	vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, triPrdP1.x, triPrdP1.y, triPrdP1.z, triPrdP2.x, triPrdP2.y, triPrdP2.z, triPrdP3.x, triPrdP3.y, triPrdP3.z,
		//	//	contact.w[1], contact.w[2], contact.w[3]);
		//	//float3 x0 = pos(vtxPos, vtxPrdP, contact.t), x1 = pos(triPos1, triPrdP1, contact.t), x2 = pos(triPos2, triPrdP2, contact.t), x3 = pos(triPos3, triPrdP3, contact.t);
		//	//printf("(%f,%f,%f) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f)\n", x0.x, x0.y, x0.z, x1.x, x1.y, x1.z, x2.x, x2.y, x2.z, x3.x, x3.y, x3.z);
		//}
		return true;
	}
	else if (VFDCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, contact, thickness, i0, i1, i2, i3))
	{
		contact.type = Contact::VF;
		//// for collision debugging
		//vfIndices.m_Data.push_back(i0);
		//vfIndices.m_Data.push_back(i1);
		//vfIndices.m_Data.push_back(i2);
		//vfIndices.m_Data.push_back(i3);
		//
		//printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		//	vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
		//	contact.w[1], contact.w[2], contact.w[3], contact.t);
		// printf("------ dcd ------\n");
		//if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
		//{
		//	printf("v:%d, f:%d, %d, %d\n", i0, i1, i2, i3);
		//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//		vtx_o.x, vtx_o.y, vtx_o.z, p1_o.x, p1_o.y, p1_o.z, p2_o.x, p1_o.y, p1_o.z, p3_o.x, p3_o.y, p3_o.z,
		//		vtx_p.x, vtx_p.y, vtx_p.z, p1_p.x, p1_p.y, p1_p.z, p2_p.x, p2_p.y, p2_p.z, p3_p.x, p3_p.y, p3_p.z,
		//		contact.w[1], contact.w[2], contact.w[3]);
		//	printf("------print after vf dcd test\n");
		//}
		//printf("vf thickness contact\n");
		return true;
	}
	else
		return false;
}

bool VFDCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, float thickness, int i0, int i1, int i2, int i3)
{
	float d = Point2Plane(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	//if (i0 == 63 && i1 == 973 && i2 == 74 && i3 == 866)
	//{
	//	printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
	//	printf("dist: %f 2 * thickness: %f\n", d, 2 * thickness);
	//	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
	//	//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
	//	//	vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	//}
	if (d - 2 * thickness >= -1e-6)
	{
		return false;
	}
	float3 weight = BarycentricCoord(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	contact.w[0] = 1;
	contact.w[1] = weight.x;
	contact.w[2] = weight.y;
	contact.w[3] = weight.z;
	contact.n = normalize(cross(normalize(p2PrdP - p1PrdP), normalize(p3PrdP - p1PrdP)));
	contact.t = 0.0f;
	bool inside;
	float edge1 = IO::Distance(p1PrdP, p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
	float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
	bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
	bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
	bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
	if (inside1 && inside2 && inside3)
	{
		inside = true;
	}
	else
	{
		inside = false;
	}
	//inside = (min(contact.w[1], contact.w[2], contact.w[3]) >= 1e-6);
	if (!inside)
		return false;
	return true;
}

bool VFCCDTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact::Type type, Contact& contact, float thickness, int i0, int i1, int i2, int i3)
{
	//if (i0 == 729 && (i1 == 108|| i1 == 646 || i1 == 333) && (i2 == 108 || i2 == 646 || i2 == 333) && (i3 == 108 || i3 == 646 || i3 == 333))
	/*if (i0 == 729 && i1 == 1171 && i2 == 108 && i3 == 333)
	{
		printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		    vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	}*/
	const float3& x0 = vtxPos, v0 = vtxPrdP - x0;  
	float3 x1 = p1Pos - x0, x2 = p2Pos - x0, x3 = p3Pos - x0; 
	float3 v1 = (p1PrdP - p1Pos) - v0, v2 = (p2PrdP - p2Pos) - v0, v3 = (p3PrdP - p3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t);
	t[nsol] = 1; // also check at end of timestep
	t[nsol + 1] = 0; // also check at start of timestep [prevent thickness DCD issue]
	std::sort(t, (t + nsol + 2));
	for (int i = 0; i < (nsol + 2); i++)
	{
		if (t[i] < 0 || t[i] > 1)
		{
			continue;
		}
		contact.t = t[i];
		float3 colliPos0 = pos(vtxPos, vtxPrdP, t[i]), colliPos1 = pos(p1Pos, p1PrdP, t[i]), colliPos2 = pos(p2Pos, p2PrdP, t[i]), colliPos3 = pos(p3Pos, p3PrdP, t[i]);
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
		float3 weight;
		contact.n = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
		//compute weight and normal
		if (type == Contact::VF)
		{
			float3 cn = normalize(cross(colliPos2 - colliPos1, colliPos3 - colliPos1));
			d = dot(colliPos0 - colliPos1, cn);
			weight = BarycentricCoord(colliPos0, colliPos1, colliPos2, colliPos3);
			contact.w[0] = 1;
			contact.w[1] = weight.x;
			contact.w[2] = weight.y;
			contact.w[3] = weight.z;
			float edge1 = IO::Distance(p1PrdP,p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
			float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
			bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >=(-epsilon));
			bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >=(-epsilon));
			bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >=(-epsilon));
			if (inside1 && inside2 && inside3)
			{
				inside = true;
			}
			else
			{
				inside = false;
			}
			//inside = (min(w[1], w[2], w[3]) >= 1e-6);
			//if (i0 == 729 && (i1 == 108 || i1 == 646 || i1 == 333) && (i2 == 108 || i2 == 646 || i2 == 333) && (i3 == 108 || i3 == 646 || i3 == 333))
			//if (i0 == 729 && i1 == 1171 && i2 == 108 && i3 == 333)
			//{
			//	printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
			//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			//		vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
			//		vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
			//		contact.t);
			//	printf("w: (%f, %f, %f)\n", weight.x, weight.y, weight.z);
			//	cout << "inside: " << inside << endl;
			//	//printf("contact t: %f\n", contact.t);
			//	printf("d :%f\n", abs(d));
			//}			
		}
		if (abs(d) < 1e-6 && inside)
		{
			return true;
		}
	}
	return false;
}

void CollisionResolve( 
	Topology meshTopol, 
	BufferVector3f& prdPBuffer,
	ContactData contactData,
	int itereations, // for later resolve iteration times
	float thickness,
	int& debugFrameId,
	BufferInt& vfIndices, 
	BufferInt& resolveTimes)
{
	PBD_DEBUG;
	resolveTimes.m_Data.resize(meshTopol.posBuffer.GetSize(), 0);
	if (debugFrameId == 1)
	{
		string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(debugFrameId) + ".cache";
		Topology temp;
		temp.indices = meshTopol.indices;
		temp.primList = meshTopol.primList;
		temp.posBuffer = meshTopol.posBuffer;
		temp.indices.SetName("Indices");
		temp.primList.SetName("primList");
		temp.posBuffer.SetName("P");
		IO::SaveToplogy(temp, path);
		printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", debugFrameId);
	}
	for (int i = 0; i < 2; ++i) // change 2 to iterations
	{
		Contact contact;
		int start, i0, i1, i2, i3;
		float3 vtxPos, p1Pos, p2Pos, p3Pos;
		auto ctxIndices = &(contactData.ctxIndices);
		auto ctxStartNum = &(contactData.ctxStartNum);
		auto ctxList = &(contactData.ctxs);
		auto posBuffer = &(meshTopol.posBuffer);
		for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
		{
			contact = ctxList->m_Data[ctxId];
			start = ctxStartNum->m_Data[ctxId].x;
			i0 = ctxIndices->m_Data[start];
			i1 = ctxIndices->m_Data[start + 1];
			i2 = ctxIndices->m_Data[start + 2];
			i3 = ctxIndices->m_Data[start + 3];
			vtxPos = posBuffer->m_Data[i0];
			p1Pos = posBuffer->m_Data[i1];
			p2Pos = posBuffer->m_Data[i2];
			p3Pos = posBuffer->m_Data[i3];
			if (Contact::VF)
			{
				// ï¿½ï¿½ï¿½ï¿½ï¿½Ç°Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½collision-freeï¿½ï¿½Î»ï¿½Ã¹ï¿½Ïµï¿½ï¿½Í¬ï¿½ï¿½ï¿½Ò¾ï¿½ï¿½ï¿½ï¿½ï¿½Ú¶ï¿½ï¿½ï¿½thicknessï¿½ï¿½ ï¿½Í²ï¿½ï¿½ï¿½ï¿½ï¿½
				float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
				bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
				float3 nn = normalize(cross(prdPBuffer.m_Data[i2] - prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i3] - prdPBuffer.m_Data[i1]));
				bool cRelatedPos = RelativePos(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], nn);
				/*
				if (i0 == 1230 || i1 == 1230 || i2 == 1230 || i3 == 1230)
				{
					auto distance = Point2Plane(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3]);
					if (cRelatedPos == fRelatedPos && distance >= 2 * thickness)
						continue;
					else
					{
						vfIndices.m_Data.push_back(i0);
						vfIndices.m_Data.push_back(i1);
						vfIndices.m_Data.push_back(i2);
						vfIndices.m_Data.push_back(i3);
						printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
						printf("distance: %f, 2*thickness: %f\n", distance, 2 * thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
						printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
							vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
							prdPBuffer.m_Data[i0].x, prdPBuffer.m_Data[i0].y, prdPBuffer.m_Data[i0].z,
							prdPBuffer.m_Data[i1].x, prdPBuffer.m_Data[i1].y, prdPBuffer.m_Data[i1].z,
							prdPBuffer.m_Data[i2].x, prdPBuffer.m_Data[i2].y, prdPBuffer.m_Data[i2].z,
							prdPBuffer.m_Data[i3].x, prdPBuffer.m_Data[i3].y, prdPBuffer.m_Data[i3].z);
						printf("------------------------------------------\n");

						VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
											prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
											contact, thickness, i0, i1, i2, i3);

						auto newDistance = Point2Plane(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3]);
						float3 nnn = normalize(cross(prdPBuffer.m_Data[i2] - prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i3] - prdPBuffer.m_Data[i1]));
						bool ncRelatedPos = RelativePos(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], nnn);
						printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;

						// for collision debugging
						debugFrameId++; 

						string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(debugFrameId) + ".cache";

						prdPBuffer.SetName("P");
						Topology temp;
						temp.indices = meshTopol.indices;
						temp.primList = meshTopol.primList;
						temp.posBuffer = prdPBuffer;
						temp.indices.SetName("Indices");
						temp.primList.SetName("primList");
						temp.posBuffer.SetName("P");
						IO::SaveToplogy(temp, path);
						printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", debugFrameId);
					}
				}
				else
				{
					auto distance = Point2Plane(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3]);
					if (cRelatedPos == fRelatedPos && distance > 2 * thickness)
						continue;
					else
					{
						VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
										   prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
										   contact, thickness, i0, i1, i2, i3);
					}
				}
				*/
				///*
				auto distance = Point2Plane(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3]);
				if (cRelatedPos == fRelatedPos && distance >= 2 * thickness)
					continue;
				else
				{
					vfIndices.m_Data.push_back(i0);
					vfIndices.m_Data.push_back(i1);
					vfIndices.m_Data.push_back(i2);
					vfIndices.m_Data.push_back(i3);
					printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
					printf("distance: %f, 2*thickness: %f\n", distance, 2 * thickness);
					cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
					printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
						prdPBuffer.m_Data[i0].x, prdPBuffer.m_Data[i0].y, prdPBuffer.m_Data[i0].z,
						prdPBuffer.m_Data[i1].x, prdPBuffer.m_Data[i1].y, prdPBuffer.m_Data[i1].z,
						prdPBuffer.m_Data[i2].x, prdPBuffer.m_Data[i2].y, prdPBuffer.m_Data[i2].z,
						prdPBuffer.m_Data[i3].x, prdPBuffer.m_Data[i3].y, prdPBuffer.m_Data[i3].z);
					printf("------------------------------------------\n");

					VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
						prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
						contact, thickness, resolveTimes, i0, i1, i2, i3);

					auto newDistance = Point2Plane(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3]);
					float3 nnn = normalize(cross(prdPBuffer.m_Data[i2] - prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i3] - prdPBuffer.m_Data[i1]));
					bool ncRelatedPos = RelativePos(prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], nnn);
					printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * thickness);
					cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;

					// for collision debugging
					debugFrameId++;

					string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(debugFrameId) + ".cache";

					prdPBuffer.SetName("P");
					Topology temp;
					temp.indices = meshTopol.indices;
					temp.primList = meshTopol.primList;
					temp.posBuffer = prdPBuffer;
					temp.indices.SetName("Indices");
					temp.primList.SetName("primList");
					temp.posBuffer.SetName("P");
					IO::SaveToplogy(temp, path);
					printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", debugFrameId);
				}
				//*/

			}
			//if (Contact::EE)
			//{
			//	EEResolve()
			//}
		}
	}
}

void VFResolve(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
	Contact contact, float thickness, BufferInt& resolveTimes,
	int i0, int i1, int i2, int i3)
{
	float3 fn = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool posRelativePos = RelativePos(vtxPos, p1Pos, fn);
	//printf("normal: (%f, %f, %f)\n", contact.n.x, contact.n.y, contact.n.z);
	//cout << "sameSide:" << sameSide << endl;
	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
	//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z, 
	//	vtxPrd.x, vtxPrd.y, vtxPrd.z, p1Prd.x, p1Prd.y, p1Prd.z, p2Prd.x, p2Prd.y, p2Prd.z, p3Prd.x, p3Prd.y, p3Prd.z);
	//printf("vtxPos: (%f, %f, %f)\n", vtxPos.x, vtxPos.y, vtxPos.z);
	//printf("p1Pos: (%f, %f, %f)\n", p1Pos.x, p1Pos.y, p1Pos.z);
	//printf("p2Pos: (%f, %f, %f)\n", p2Pos.x, p2Pos.y, p2Pos.z);
	//printf("p3Pos: (%f, %f, %f))\n", p3Pos.x, p3Pos.y, p3Pos.z);
	//printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	//printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	//printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	//printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
	float3 n = normalize(cross(p2Prd - p1Prd, p3Prd - p1Prd));
	float depth = Point2Plane(vtxPrd, p1Prd, p2Prd, p3Prd);
	bool prdPRelativePos = RelativePos(vtxPrd, p1Prd, n);
	float dp = 0.0f;
	// ï¿½ï¿½prdpï¿½ï¿½ï¿½Â¼ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½weight
	float3 weight = BarycentricCoord(vtxPrd, p1Prd, p2Prd, p3Prd);
	if (posRelativePos == prdPRelativePos && depth < 2 * thickness) // vf dcd
	{
		dp = (2 * thickness - depth) * 0.5;
		printf("---------vf dcd----------\n");
	}
	else if (posRelativePos != prdPRelativePos) // relativePos different might be ccd
	{
		if (resolveTimes.m_Data[i0] == 0 || resolveTimes.m_Data[i1] == 0 || resolveTimes.m_Data[i2] == 0 || resolveTimes.m_Data[i3] == 0) // haven't been resolve before resove as ccd directly
		{
			dp = (depth + 2 * thickness) * 0.5;
			//if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
				//		printf("---------vf ccd----------\n");
		}
		else if (VFCCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrd, p1Prd, p2Prd, p3Prd, Contact::VF, contact, thickness, i0, i1, i2, i3)) // have been resolve before then test again
		{
			dp = (depth + 2 * thickness) * 0.5;
			//if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
				//		printf("---------vf ccd----------\n");
		}
	}
	if (dp != 0.0f)
	{
		resolveTimes.m_Data[i0] += 1;
		resolveTimes.m_Data[i1] += 1;
		resolveTimes.m_Data[i2] += 1;
		resolveTimes.m_Data[i3] += 1;
	}

	float sw = weight.x * weight.x + weight.y * weight.y + weight.z * weight.z;
	//if (i0 == 3806 || i1 == 3806 || i2 == 3806 || i3 == 3806)
	//{
	//	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z);
	//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\n",
	//		vtxPrd.x, vtxPrd.y, vtxPrd.z, p1Prd.x, p1Prd.y, p1Prd.z, p2Prd.x, p2Prd.y, p2Prd.z, p3Prd.x, p3Prd.y, p3Prd.z,
	//		depth, contact.w[1], contact.w[2], contact.w[3], Contact::VFDCD, posRelativePos);
	//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
	//		vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
	//		afterProjPrdpBuffer.m_Data[i0].x, afterProjPrdpBuffer.m_Data[i0].y, afterProjPrdpBuffer.m_Data[i0].z,
	//		afterProjPrdpBuffer.m_Data[i1].x, afterProjPrdpBuffer.m_Data[i1].y, afterProjPrdpBuffer.m_Data[i1].z,
	//		afterProjPrdpBuffer.m_Data[i2].x, afterProjPrdpBuffer.m_Data[i2].y, afterProjPrdpBuffer.m_Data[i2].z,
	//		afterProjPrdpBuffer.m_Data[i3].x, afterProjPrdpBuffer.m_Data[i3].y, afterProjPrdpBuffer.m_Data[i3].z,
	//		contact.t);
	//}
	if (posRelativePos)
	{
		vtxPrd += n * dp;
		p1Prd += -n * dp * weight.x / sw;
		p2Prd += -n * dp * weight.y / sw;
		p3Prd += -n * dp * weight.z / sw;

	}
	else
	{
		vtxPrd += -n * dp;
		p1Prd += n * dp * weight.x / sw;
		p2Prd += n * dp * weight.y / sw;
		p3Prd += n * dp * weight.z / sw;
	}
}
// ---------new method--------------
void CollisionResolveNew(
	Topology meshTopol,
	BufferVector3f& prdPBuffer,
	ContactData contactData,
	int itereations, // for later resolve iteration times
	float thickness,
	int& debugFrameId,
	BufferInt& vfIndices)
{
	PBD_DEBUG;

	for (int i = 0; i < itereations; ++i) // change 2 to iterations
	{
		auto ctxIndices = &(contactData.ctxIndices);
		auto ctxStartNum = &(contactData.ctxStartNum);
		auto ctxList = &(contactData.ctxs);
		auto posBuffer = &(meshTopol.posBuffer);
		for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
		{
			Contact contact;
			// For VF
			int start, i0, i1, i2, i3;
			float3 vtxPos, p1Pos, p2Pos, p3Pos;
			//// For EE
			//float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
			//int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges

			contact = ctxList->m_Data[ctxId];
			start = ctxStartNum->m_Data[ctxId].x;
			i0 = ctxIndices->m_Data[start];
			i1 = ctxIndices->m_Data[start + 1];
			i2 = ctxIndices->m_Data[start + 2];
			i3 = ctxIndices->m_Data[start + 3];
			vtxPos = posBuffer->m_Data[i0];
			p1Pos = posBuffer->m_Data[i1];
			p2Pos = posBuffer->m_Data[i2];
			p3Pos = posBuffer->m_Data[i3];
			if (contact.type == Contact::VF)
			{
				if (VFResolveTest(
					vtxPos, p1Pos, p2Pos, p3Pos,
					prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
					contact,
					i0, i1, i2, i3, thickness))
				{
					VFResolveNew(vtxPos, p1Pos, p2Pos, p3Pos,
						prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
						contact, thickness, i0, i1, i2, i3);
				}
				else
					continue;
			}

		}

		int resolveEECounter = 0;
		for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
		{
			Contact contact;
			// For EE
			float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
			int start, e0, e1, e2, e3;  // e0-e1; e2-e3  two edges

			contact = ctxList->m_Data[ctxId];
			start = ctxStartNum->m_Data[ctxId].x;
			
			if (contact.type == Contact::EE)
			{
				e0 = ctxIndices->m_Data[start];
				e1 = ctxIndices->m_Data[start + 1];
				e2 = ctxIndices->m_Data[start + 2];
				e3 = ctxIndices->m_Data[start + 3];
				e0Pos = posBuffer->m_Data[e0];
				e1Pos = posBuffer->m_Data[e1];
				e2Pos = posBuffer->m_Data[e2];
				e3Pos = posBuffer->m_Data[e3];

				float constrVal = 0.0f;
				//bool debugThis = true;

				if (EEResolveTest(
					prdPBuffer.m_Data[e0],
					prdPBuffer.m_Data[e1],
					prdPBuffer.m_Data[e2],
					prdPBuffer.m_Data[e3],
					&constrVal,
					contact,
					thickness))
				{
					resolveEECounter++;  // start from 1

					// save topol before resolve
					//string path = "D://0319CCDTest//VFEESimple//EEResolveCache//Before_ee." + to_string(resolveEECounter) + ".cache";
					//Topology temp;
					//temp.indices = meshTopol.indices;
					//temp.primList = meshTopol.primList;
					//temp.posBuffer = prdPBuffer;
					//temp.indices.SetName("Indices");
					//temp.primList.SetName("primList");
					//temp.posBuffer.SetName("P");
					//IO::SaveToplogy(temp, path);

					
					//string beforePath = "D://0319CCDTest//VFEESimple//before//beforeEE." + to_string(resolveEECounter) + ".cache";
					//SaveEEContact(contactData, ctxId, prdPBuffer, beforePath);
					//if (resolveEECounter == 1)
					//{
					//	beforePath = "D://0319CCDTest//VFEESimple//before//beforeEE.-1.cache";
					//	SaveEEContact(contactData, ctxId, prdPBuffer, beforePath);
					//	beforePath = "D://0319CCDTest//VFEESimple//before//beforeEE.0.cache";
					//	SaveEEContact(contactData, ctxId, prdPBuffer, beforePath);
					//}

					//if (debugThis)
					//{
					//	printf("This contact should be resolved!\n");
					//	printf("\te0-e1 e2-e3: : %d,%d,%d,%d\n", e0, e1, e2, e3);
					//	printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
					//	printf("\tconstrVal: %f\n", constrVal);
					//}

					// resolve
					EEResolve(
						prdPBuffer.m_Data[e0],
						prdPBuffer.m_Data[e1],
						prdPBuffer.m_Data[e2],
						prdPBuffer.m_Data[e3],
						&constrVal,
						contact);

					//float c = calcRestEEConstr(
					//	prdPBuffer.m_Data[e0],
					//	prdPBuffer.m_Data[e1],
					//	prdPBuffer.m_Data[e2],
					//	prdPBuffer.m_Data[e3], 
					//	contact);
					//printf("after resolve constr %f\n");

					//string afterPath = "D://0319CCDTest//VFEESimple//after//afterEE." + to_string(resolveEECounter) + ".cache";
					//SaveEEContact(contactData, ctxId, prdPBuffer, afterPath);
					//if (resolveEECounter == 1)
					//{
					//	afterPath = "D://0319CCDTest//VFEESimple//after//afterEE.-1.cache";
					//	SaveEEContact(contactData, ctxId, prdPBuffer, afterPath);
					//	afterPath = "D://0319CCDTest//VFEESimple//after//afterEE.0.cache";
					//	SaveEEContact(contactData, ctxId, prdPBuffer, afterPath);
					//}

					//printf("\n");

					//// save topol after resolve
					//path = "D://0319CCDTest//VFEESimple//EEResolveCache//After_ee." + to_string(resolveEECounter) + ".cache";
					//temp.indices = meshTopol.indices;
					//temp.primList = meshTopol.primList;
					//temp.posBuffer = prdPBuffer;
					//temp.indices.SetName("Indices");
					//temp.primList.SetName("primList");
					//temp.posBuffer.SetName("P");
					//IO::SaveToplogy(temp, path);
				}
			}
		}
		/*for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
		{
			Contact contact;
			// For VF
			int start, i0, i1, i2, i3;
			float3 vtxPos, p1Pos, p2Pos, p3Pos;
			// For EE
			float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
			int e0, e1, e2, e3;  // e0-e1; e2-e3  two edges

			contact = ctxList->m_Data[ctxId];
			start = ctxStartNum->m_Data[ctxId].x;
			i0 = ctxIndices->m_Data[start];
			i1 = ctxIndices->m_Data[start + 1];
			i2 = ctxIndices->m_Data[start + 2];
			i3 = ctxIndices->m_Data[start + 3];
			vtxPos = posBuffer->m_Data[i0];
			p1Pos = posBuffer->m_Data[i1];
			p2Pos = posBuffer->m_Data[i2];
			p3Pos = posBuffer->m_Data[i3];
			if (contact.type == Contact::VF)
			{
				if (VFResolveTest(
					vtxPos, p1Pos, p2Pos, p3Pos,
				    prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
					contact,
					i0, i1, i2, i3, thickness))
				{
					VFResolveNew(vtxPos, p1Pos, p2Pos, p3Pos,
						prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
						contact, thickness, i0, i1, i2, i3);
					vfIndices.m_Data.push_back(i0);
					vfIndices.m_Data.push_back(i1);
					vfIndices.m_Data.push_back(i2);
					vfIndices.m_Data.push_back(i3);
				}
				else
					continue;
			}
			else if (contact.type == Contact::EE)
			{
				start = ctxStartNum->m_Data[ctxId].x;
				e0 = ctxIndices->m_Data[start];
				e1 = ctxIndices->m_Data[start + 1];
				e2 = ctxIndices->m_Data[start + 2];
				e3 = ctxIndices->m_Data[start + 3];
				e0Pos = posBuffer->m_Data[e0];
				e1Pos = posBuffer->m_Data[e1];
				e2Pos = posBuffer->m_Data[e2];
				e3Pos = posBuffer->m_Data[e3];

				float constrVal = 0.0f;
				if (EEResolveTest(
					prdPBuffer.m_Data[e0],
					prdPBuffer.m_Data[e1],
					prdPBuffer.m_Data[e2],
					prdPBuffer.m_Data[e3],
					&constrVal,
					contact,
					thickness))
				{
					printf("This contact should be resolved!\n");
					printf("\te0-e1: %d-%d; e2-e3: %d-%d\n", e0, e1, e2, e3);
					printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);
					printf("\tconstrVal: %f\n", constrVal);
					EEResolve(
						prdPBuffer.m_Data[e0],
						prdPBuffer.m_Data[e1],
						prdPBuffer.m_Data[e2],
						prdPBuffer.m_Data[e3],
						&constrVal,
						contact);
				}
			}
			else
			{
				;   // impossible
			}
		}*/
	}
	//Topology temp;
	//temp.indices = meshTopol.indices;
	//temp.primList = meshTopol.primList;
	//temp.posBuffer = prdPBuffer;
	//temp.indices.SetName("Indices");
	//temp.primList.SetName("primList");
	//temp.posBuffer.SetName("P");
	//string path = "D://0319CCDTest//EETests//EETest3//EEResolve3.cache";
	//IO::SaveToplogy(temp, path);
}

bool VFResolveTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3, float thickness)
{
	//if(i0==1010 && i1 == 1084 && i2 == 212 && i3 == 345)
	float3 colliFreeN = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool colliFreeDir = RelativePos(vtxPos, p1Pos, colliFreeN);
	//bool colliFreeDir = contact.colliFreeDir;
	float3 prdPN = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
	float3 prdPn = cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP);
	if (dot(prdPn, prdPn) <= 1e-6) // ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ë»ï¿½ï¿½ï¿½Ò»ï¿½ï¿½ï¿½ß£ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Î£ï¿½ï¿½ï¿½ï¿½ï¿½
	{

		return false;
	}
	float dist = dot(vtxPrdP - p1PrdP, prdPN * (colliFreeDir ? 1 : -1)) - 2.0 * thickness;
	/*if (i0 == 50)
	{
		printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
		printf("dist: %f\n", dist);
		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
					vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	}*/
	if (dist >= 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void VFResolveNew(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrdP, float3& p1PrdP, float3& p2PrdP, float3& p3PrdP,
	Contact contact, float thickness,
	int i0, int i1, int i2, int i3)
{
	float3 vp = vtxPrdP, p1p = p1PrdP, p2p = p2PrdP, p3p = p3PrdP;
	//if (debug == 7&& iteration >59)
	//{
	//	//if (i0 == 4608 || i1 == 4608 || i2 == 4608 || i3 == 4608)
	//	//{
	//		printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
	//		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
	//			vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
	//			vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	//	//}
	//}
	float3 crossColliFreeN = cross(p2Pos - p1Pos, p3Pos - p1Pos);
	float3 colliFreeN = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool colliFreeDir = RelativePos(vtxPos, p1Pos, colliFreeN);
	//bool colliFreeDir = contact.colliFreeDir;
	float3 crossPrdPN = cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP);
	float3 prdPN = normalize(cross(p2PrdP - p1PrdP, p3PrdP - p1PrdP));
	float3 p1ProjPrdP = p1PrdP + prdPN * (colliFreeDir ? 1 : -1) * 2.0 * thickness;
	float depth = dot(p1ProjPrdP - vtxPrdP, (prdPN * (colliFreeDir ? 1 : -1))) * 0.5; // must > 0
	float3 vtxFix = depth * (prdPN * (colliFreeDir ? 1 : -1));
	float3 p1Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	float3 p2Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	float3 p3Fix = depth * (prdPN * (colliFreeDir ? -1 : 1));
	/*
	fixedBuffer.m_Data[i0] += vtxFix;
	fixedBuffer.m_Data[i1] += p1Fix;
	fixedBuffer.m_Data[i2] += p2Fix;
	fixedBuffer.m_Data[i3] += p3Fix;
	vFixedBuffer.m_Data[i0] += vtxFix;
	fFixedBuffer.m_Data[i1] += p1Fix;
	fFixedBuffer.m_Data[i2] += p2Fix;
	fFixedBuffer.m_Data[i3] += p3Fix;
	*/
	float3 weight = BarycentricCoord(vtxPrdP, p1PrdP, p2PrdP, p3PrdP);
	float sw = weight.x * weight.x + weight.y * weight.y + weight.z * weight.z;
	vtxPrdP += vtxFix;
	p1PrdP += p1Fix * weight.x / sw;
	p2PrdP += p2Fix * weight.y / sw;
	p3PrdP += p3Fix * weight.z / sw;
	if (isnan(vtxPrdP.x) || isnan(vtxPrdP.y) || isnan(vtxPrdP.z) || isnan(p1PrdP.x) || isnan(p1PrdP.y) || isnan(p1PrdP.z) ||
		isnan(p2PrdP.x) || isnan(p2PrdP.y) || isnan(p2PrdP.z) || isnan(p3PrdP.x) || isnan(p3PrdP.y) || isnan(p3PrdP.z))
	{
		printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
		printInfo("weight:", weight);
		printInfo("colliFreeN cross:", crossColliFreeN);
		printInfo("prdpN cross:", crossPrdPN);
		printf("------ depth:%f, sw: %f \n", depth, sw);
		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
			vp.x, vp.y, vp.z, p1p.x, p1p.y, p1p.z, p2p.x, p2p.y, p2p.z, p3p.x, p3p.y, p3p.z);
	}
}


void getEdgesOfTri(int tp0, int tp1, int tp2, int2* rtnEdges)
{
	rtnEdges[0] = make_int2(tp0, tp1);
	rtnEdges[1] = make_int2(tp1, tp2);
	rtnEdges[2] = make_int2(tp2, tp0);
}

double signed_ee_distance(const float3& x0, const float3& x1,
	const float3& y0, const float3& y1,
	float3* n, double* w)
{
	// Does NOT consider parallel lines
	float3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	if (norm2(*n) < 1e-6) // parallel
		return infinity;  // not return infinity
	*n = normalize(*n);
	double h = dot(x0 - y0, *n);
	double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
		b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	w[0] = a0 / (a0 + a1);
	w[1] = a1 / (a0 + a1);
	w[2] = -b0 / (b0 + b1);
	w[3] = -b1 / (b0 + b1);
	// printf("w[0]: %f; w[1]: %f; w[2]: %f; w[3]: %f\n", w[0], w[1], w[2], w[3]);
	return h;

	// Deal with parallel lines
	//float3 _n; if (!n) n = &_n;
	//double _w[4]; if (!w) w = _w;
	//*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	//if (norm2(*n) < 1e-8) {
	//	// special case: parallel lines
	//	float3 e0 = normalize(x1 - x0), e1 = normalize(y1 - y0);

	//	double p0min = dot(x0, e0), p0max = dot(x1, e0), p1min = dot(y0, e0), p1max = dot(y1, e0);
	//	if (p1max < p1min) swap(p1max, p1min);

	//	double a = max(p0min, p1min), b = min(p0max, p1max), c = 0.5*(a + b);
	//	if (a > b) return infinity;

	//	float3 d = (y0 - x0) - dot(y0 - x0, e0)*e0;

	//	if (n) *n = normalize(-d);
	//	if (w) {
	//		w[1] = (c - dot(x0, e0)) / norm(x1 - x0);
	//		w[0] = 1.0 - w[1];
	//		w[3] = -(dot(e0, e1)*c - dot(y0, e1)) / norm(y1 - y0);
	//		w[2] = -1.0 - w[3];
	//	}
	//	return norm(d);
	//}
	//*n = normalize(*n);
	//double h = dot(x0 - y0, *n);
	//double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
	//	b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	//w[0] = a0 / (a0 + a1);
	//w[1] = a1 / (a0 + a1);
	//w[2] = -b0 / (b0 + b1);
	//w[3] = -b1 / (b0 + b1);
	//return h;
}

float calcRestEEConstr(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos, Contact& contact)
{
	float3 currNormal = cross(normalize(e1Pos - e0Pos), normalize(e3Pos - e2Pos));
	if (norm2(currNormal) <= 1e-6)
	{
		//printf("calcRestEEConstr get parallel\n");
		return 0.0f;
	}
	currNormal = normalize(currNormal);
	contact.n = currNormal;
	//printf("normal is %f %f %f\n", currNormal.x, currNormal.y, currNormal.z);
	double a0 = stp(e3Pos - e1Pos, e2Pos - e1Pos, currNormal), a1 = stp(e2Pos - e0Pos, e3Pos - e0Pos, currNormal),
		b0 = stp(e0Pos - e3Pos, e1Pos - e3Pos, currNormal), b1 = stp(e1Pos - e2Pos, e0Pos - e2Pos, currNormal);
	contact.w[0] = a0 / (a0 + a1);
	contact.w[1] = a1 / (a0 + a1);
	contact.w[2] = -b0 / (b0 + b1);
	contact.w[3] = -b1 / (b0 + b1);
	float3 dist = (contact.w[0] * e0Pos + contact.w[1] * e1Pos) - ((-contact.w[2]) * e2Pos + (-contact.w[3]) * e3Pos);
	//printf("dist is %f %f %f\n", dist.x, dist.y, dist.z);
	float c = dot(currNormal, dist);
	
	//printf("\e0Pos: (%f, %f, %f); e1Pos: (%f, %f, %f); e2Pos: (%f, %f, %f); e3Pos: (%f, %f, %f)\n",
	//	e0Pos.x, e0Pos.y, e0Pos.z,
	//	e1Pos.x, e1Pos.y, e1Pos.z,
	//	e2Pos.x, e2Pos.y, e2Pos.z,
	//	e3Pos.x, e3Pos.y, e3Pos.z);
	return c;
}

float calcCurrEEConstr(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, Contact& contact, float thickness)
{
	float3 currNormal = cross(normalize(e1PrdP - e0PrdP), normalize(e3PrdP - e2PrdP));
	if (norm2(currNormal) <= 1e-6)
	{
		//printf("calcCurrEEConstr get parallel\n");
		return 0.0f;
	}
	currNormal = normalize(currNormal);
	contact.n = currNormal;
	//printf("normal is %f %f %f\n", currNormal.x, currNormal.y, currNormal.z);
	double a0 = stp(e3PrdP - e1PrdP, e2PrdP - e1PrdP, currNormal), a1 = stp(e2PrdP - e0PrdP, e3PrdP - e0PrdP, currNormal),
		b0 = stp(e0PrdP - e3PrdP, e1PrdP - e3PrdP, currNormal), b1 = stp(e1PrdP - e2PrdP, e0PrdP - e2PrdP, currNormal);
	contact.w[0] = a0 / (a0 + a1);
	contact.w[1] = a1 / (a0 + a1);
	contact.w[2] = -b0 / (b0 + b1);
	contact.w[3] = -b1 / (b0 + b1);

	float3 dist = (contact.w[0] * e0PrdP + contact.w[1] * e1PrdP) - ((-contact.w[2]) * e2PrdP + (-contact.w[3]) * e3PrdP);
	//printf("dist is %f %f %f\n", dist.x, dist.y, dist.z);
	float c = dot(currNormal * (contact.colliFreeDir ? 1 : -1), dist) - 2.0*thickness;

	//printf("\tEEResolveTest\n");
	//printf("colliFreeDir is %d\n", contact.colliFreeDir);
	//printf("\tprdPN: (%f, %f, %f)\n", contact.n.x, contact.n.y, contact.n.z);
	//printf("\te0PrdP: (%f, %f, %f); e1PrdP: (%f, %f, %f); e2PrdP: (%f, %f, %f); e3PrdP: (%f, %f, %f)\n",
	//	e0PrdP.x, e0PrdP.y, e0PrdP.z,
	//	e1PrdP.x, e1PrdP.y, e1PrdP.z,
	//	e2PrdP.x, e2PrdP.y, e2PrdP.z,
	//	e3PrdP.x, e3PrdP.y, e3PrdP.z);
	//printf("\ta0: %f; a1: %f; b0: %f; b1: %f\n", a0, a1, b0, b1);
	//printf("\tweight %f, %f, %f, %f\n", contact.w[0], contact.w[1], contact.w[2], contact.w[3]);


	return c;
}

void initContactFlag(float constrVal, Contact& contact)
{
	contact.colliFreeDir = constrVal > 0.0 ? true : false;
	// printf("colliFreeDir is %d\n", contact.colliFreeDir);
}

bool EETest(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos,
	float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP,
	Contact& contact, float thickness, int e0, int e1, int e2, int e3)
{

	if (e0 == e2 || e0 == e3 || e1 == e2 || e1 == e3)
		return false;
	if (EECCDTest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness, e0, e1, e2, e3))
	{
		contact.type = Contact::EE;
		return true;
	}
	//else if (EEDCDTest(e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, e0, e1, e2, e3))
	//{
	//	printf("EE thickness contact\n");
	//	contact.type = Contact::EEDCD;
	//	return true;
	//}
	else
		return false;
}

bool EECCDTest(float3 e0Pos, float3 e1Pos, float3 e2Pos, float3 e3Pos,
	float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP,
	Contact& contact, float thickness, int e0, int e1, int e2, int e3)
{
	const float3& x0 = e0Pos, v0 = e0PrdP - x0;
	float3 x1 = e1Pos - x0, x2 = e2Pos - x0, x3 = e3Pos - x0; // x: p's pos - V's pos
	float3 v1 = (e1PrdP - e1Pos) - v0, v2 = (e2PrdP - e2Pos) - v0, v3 = (e3PrdP - e3Pos) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3),
		a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3),
		a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3),
		a3 = stp(v1, v2, v3);

	if (abs(a0) < 1e-6 * norm(x1) * norm(x2) * norm(x3))
		return false; // initially coplanar

	double t[5];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	t[nsol + 1] = 0; // also check at start of timestep [prevent thickness DCD issue]
	std::sort(t, (t + nsol + 2));
	for (int i = 0; i < (nsol + 2); i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
		float3 colliPos0 = pos(e0Pos, e0PrdP, t[i]), colliPos1 = pos(e1Pos, e1PrdP, t[i]), colliPos2 = pos(e2Pos, e2PrdP, t[i]), colliPos3 = pos(e3Pos, e3PrdP, t[i]);
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
		d = signed_ee_distance(colliPos0, colliPos1, colliPos2, colliPos3, &n, w);
		float edge1 = IO::Distance(e0PrdP, e1PrdP), edge2 = IO::Distance(e2PrdP, e3PrdP); // Distance funtion declear many times
		if (__USE_EPSILON)
		{
			float edge1 = IO::Distance(e0PrdP, e1PrdP), edge2 = IO::Distance(e2PrdP, e3PrdP); // Distance funtion declear many times
			float epsilon = thickness / ((edge1 + edge2) / 2.0f);
			bool inside1 = (min(w[0], w[1]) >= -epsilon) && (max(w[0], w[1]) <= 1 + epsilon);
			bool inside2 = (min(-w[2], -w[3]) >= -epsilon) && (max(-w[2], -w[3]) <= 1 + epsilon);

			inside = inside1 && inside2;
		}
		else
			inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);
		// printf("EE: d is %f; inside is %d\n", d, inside);
		// TODO: check this again for EE resolve
		if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > -1e-6)  // calculate here for ee resolve
			n = -n;
		if (abs(d) < (2.0 * thickness) && inside)   // add the thickness for preventing EE DCD
		{
			float c = calcRestEEConstr(e0Pos, e1Pos, e2Pos, e3Pos, contact);

			//if (e0 == 74 && e1 == 36 && e2 == 62 && e3 == 40)
			//	printf("dubug ee info: d is %f\n", abs(d));
			// assert(c != 0.0);
			if (c == 0.0)
			{
				// printf("d and inside but constr is 0: %d, %d, %d, %d\n", e0, e1, e2, e3);
				return false;
			}

			//printf("Adding EE: %d,%d,%d,%d ----- c: %f\n", e0, e1,e2,e3, c);
			initContactFlag(c, contact);  // initialize contact colli free dir before return
			return true;
		}
	}
	return false;
}

bool EEResolveTest(float3 e0PrdP, float3 e1PrdP, float3 e2PrdP, float3 e3PrdP, float* constrVal, Contact& contact, float thickness)
{
	//float3 prdPN = cross(normalize(e1PrdP - e0PrdP), normalize(e3PrdP - e2PrdP));
	//if (norm2(prdPN) <= 1e-6)
	//	return false;
	//prdPN = normalize(prdPN);
	//contact.n = prdPN;  // update normal of the contact

	//double a0 = stp(e3PrdP - e1PrdP, e2PrdP - e1PrdP, prdPN), a1 = stp(e2PrdP - e0PrdP, e3PrdP - e0PrdP, prdPN),
	//	b0 = stp(e0PrdP - e3PrdP, e1PrdP - e3PrdP, prdPN), b1 = stp(e1PrdP - e2PrdP, e0PrdP - e2PrdP, prdPN);
	//contact.w[0] = a0 / (a0 + a1);
	//contact.w[1] = a1 / (a0 + a1);
	//contact.w[2] = -b0 / (b0 + b1);
	//contact.w[3] = -b1 / (b0 + b1);

	//// TODO: clamp the weights
	//// Recalculate current constr value
	////*constrVal = calcEEConstr(e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness);
	////printf("function constr %f\n", *constrVal);
	//*constrVal = dot(prdPN * (contact.colliFreeDir ? -1 : 1), (contact.w[0] * e0PrdP + contact.w[1] * e1PrdP) - ((-contact.w[2]) * e2PrdP + (-contact.w[3]) * e3PrdP)) - 2.0 * thickness;
	//printf("constrVal1: %f -------------- ", *constrVal);
	*constrVal = calcCurrEEConstr(e0PrdP, e1PrdP, e2PrdP, e3PrdP, contact, thickness);
	// printf("constrVal2: %f\n ", *constrVal);

	if ((*constrVal) >= 0.0f)
	{
		return false;
	}
	else
		return true;

}

void EEResolve(float3& e0PrdP, float3& e1PrdP, float3& e2PrdP, float3& e3PrdP, float* constrVal, Contact& contact)
{
	float sumW = contact.w[0] * contact.w[0] + contact.w[1] * contact.w[1] + contact.w[2] * contact.w[2] + contact.w[3] * contact.w[3];
	// contact.n should be the prdP Normal now
	float3 k = (-(*constrVal) * contact.n * (contact.colliFreeDir ? 1 : -1)) * (1.0 / sumW);
	e0PrdP += (contact.w[0] * k);
	e1PrdP += (contact.w[1] * k);
	e2PrdP += (contact.w[2] * k);
	e3PrdP += (contact.w[3] * k);
	//if (debugThis) 
	//{
	//	printf("\tEEResolve:\n");
	//	printf("\tsumW: %f; k: (%f, %f, %f)\n", sumW, k.x, k.y, k.z);
	//	printf("\te0PrdP: (%f, %f, %f); e1PrdP: (%f, %f, %f); e2PrdP: (%f, %f, %f); e3PrdP: (%f, %f, %f)\n",
	//		e0PrdP.x, e0PrdP.y, e0PrdP.z,
	//		e1PrdP.x, e1PrdP.y, e1PrdP.z,
	//		e2PrdP.x, e2PrdP.y, e2PrdP.z,
	//		e3PrdP.x, e3PrdP.y, e3PrdP.z);
	//}
	
}

void SaveEEContact(ContactData& contactData, int ctxId, BufferVector3f& prdPBuffer, string path)
{
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;
	float3 e0Pos, e1Pos, e2Pos, e3Pos, e0PrdP, e1PrdP, e2PrdP, e3PrdP;
	int start, e0, e1, e2, e3;  // e0-e1; e2-e3  two edges
	Contact contact;

	auto ctxIndices = &(contactData.ctxIndices);
	auto ctxStartNum = &(contactData.ctxStartNum);
	auto ctxList = &(contactData.ctxs);

	contact = ctxList->m_Data[ctxId];
	start = ctxStartNum->m_Data[ctxId].x;

	e0 = ctxIndices->m_Data[start];
	e1 = ctxIndices->m_Data[start + 1];
	e2 = ctxIndices->m_Data[start + 2];
	e3 = ctxIndices->m_Data[start + 3];
	e0Pos = prdPBuffer.m_Data[e0];
	e1Pos = prdPBuffer.m_Data[e1];
	e2Pos = prdPBuffer.m_Data[e2];
	e3Pos = prdPBuffer.m_Data[e3];

	float3 pOne0 = contact.w[0] * e0Pos + contact.w[1] * e1Pos;
	float3 pOne1 = -contact.w[2] * e2Pos + (-contact.w[3]) * e3Pos;

	float dist = dot(pOne0 - pOne1, contact.n * (contact.colliFreeDir ? -1 : 1));

	printf("Saving EE Constact...\n");

	// e0,e1:e2,e3-e0.x,e0.y,e0.z,e1.x,e1.y,e1.z,e2.x,e2.y,e2.z,e3.x,e3.y,e3.z;
	// pOne0.x,pOne0.y,pOne0.z,pOne1.x,pOne1.y,pOne1.z,
	// dist

	ofs << e0 << "," << e1 << ":" << e2 << "," << e3 << "/";
	ofs << e0Pos.x << "," << e0Pos.y << "," << e0Pos.z << ",";
	ofs << e1Pos.x << "," << e1Pos.y << "," << e1Pos.z << ",";
	ofs << e2Pos.x << "," << e2Pos.y << "," << e2Pos.z << ",";
	ofs << e3Pos.x << "," << e3Pos.y << "," << e3Pos.z << ";";
	
	ofs << std::endl;

	// adding hit point
	ofs << pOne0.x << "," << pOne0.y << "," << pOne0.z << ",";
	ofs << pOne1.x << "," << pOne1.y << "," << pOne1.z << ",";

	ofs << std::endl;

	ofs << dist;

	ofs << std::endl;
}

void SaveContact(
	string path,
	std::map<int, std::set<int>>& contactVFMap,
	std::map<int, BufferVector3f>& contactVHitMap)
{
	std::ofstream ofs(path);
	if (!ofs.is_open())
		return;
	printf("Saving Constact...\n");
	printf("contactVFMap size %d, contactVHitMap size %d\n", contactVFMap.size(), contactVHitMap.size());
	assert(contactVFMap.size() == contactVHitMap.size());
	map<int, std::set<int>>::iterator it;
	map<int, BufferVector3f>::iterator it2;
	for (it = contactVFMap.begin(), it2 = contactVHitMap.begin();
		it != contactVFMap.end() && it2 != contactVHitMap.end();
		it++, it2++)
	{
		auto vertexId = it->first;
		auto triIds = &(contactVFMap[vertexId]);
		auto hitPos = &(contactVHitMap[vertexId]);
		assert(triIds->size() == hitPos->GetSize());
		ofs << vertexId << ":";
		int i = 0;
		for (std::set<int>::iterator itSet = triIds->begin();
			itSet != triIds->end() && i < hitPos->GetSize();
			++itSet, ++i)
		{
			float3 currHitPos = hitPos->m_Data[i];
			ofs << (*itSet) << ",";
			ofs << currHitPos.x << "," << currHitPos.y << "," << currHitPos.z << ",";
		}
		ofs << ";";
	}
	ofs << std::endl;

	/*for (it2 = m_contactVHitMap.begin(); it2 != m_contactVHitMap.end(); it2++)
	{
		auto vertexId = it2->first;
		auto hitPos = &(m_contactVHitMap[vertexId]);
		ofs << vertexId << ":";
		for (int i = 0; i < hitPos->GetSize(); ++i)
		{
			float3 currHitPos = hitPos->m_Data[i];
			printf("%f %f %f ", currHitPos.x, currHitPos.y, currHitPos.z);
			ofs << currHitPos.x << "," << currHitPos.y << "," << currHitPos.z;
		}
		ofs << ";";
	}
	ofs << std::endl;*/
}

void SavePrdPBeforeCCD(string path, Topology meshTopol, BufferVector3f prdP)
{
	Topology temp;
	temp.indices = meshTopol.indices;
	temp.primList = meshTopol.primList;
	temp.posBuffer = prdP;
	temp.indices.SetName("Indices");
	temp.primList.SetName("primList");
	temp.posBuffer.SetName("P");
	IO::SaveToplogy(temp, path);
}

void SavePrdPAfterCCD(string path, Topology meshTopol, BufferVector3f prdP)
{
	Topology temp;
	temp.indices = meshTopol.indices;
	temp.primList = meshTopol.primList;
	temp.posBuffer = prdP;
	temp.indices.SetName("Indices");
	temp.primList.SetName("primList");
	temp.posBuffer.SetName("P");
	IO::SaveToplogy(temp, path);
}

// ------------- for ccdtest readfrom txt ----------------
void readMeshFromTxt(string filename, Topology& topol)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = UT::splitString(buffer, firstSplit)[1];
		auto dataStringList =UT::splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				// std::cout << dataStringList[i] << std::endl;
				float3 vertexPos = make_float3(std::stof(dataStringList[i]),
					std::stof(dataStringList[i + 1]),
					std::stof(dataStringList[i + 2]));
				// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				topol.posBuffer.m_Data.push_back(vertexPos);
			}
			std::cout << "topol.posBuffer: " << topol.posBuffer.GetSize() << endl;
			assert(topol.posBuffer.GetSize() == dataStringList.size() / 3);
		}
		else  // second line of the input file: vertices tet
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i]));
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i + 1]));
				topol.indices.m_Data.push_back(std::stoi(dataStringList[i + 2]));
			}
		}
		++lineNum;
	}
	in.close();
	printf("Read File Done!\n");
	for (int i = 0; i < topol.indices.GetSize() / 3; ++i)
	{
		int2 p;
		p.x = i * 3;
		p.y = 3;
		topol.primList.m_Data.push_back(p);
	}
	printf("Init PrimList Done!\n");
}

void readBufferFromTxt(string filename, BufferVector3f& prdPBuffer)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = UT::splitString(buffer, firstSplit)[1];
		auto dataStringList = UT::splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (uint i = 0; i < dataStringList.size(); i += 3)
			{
				// std::cout << dataStringList[i] << std::endl;
				float3 vertexPos = make_float3(std::stof(dataStringList[i]),
					std::stof(dataStringList[i + 1]),
					std::stof(dataStringList[i + 2]));
				// printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				prdPBuffer.m_Data.push_back(vertexPos);
			}
			std::cout << "topol.posBuffer: " << prdPBuffer.GetSize() << endl;
			assert(prdPBuffer.GetSize() == dataStringList.size() / 3);
		}
		lineNum++;
	}
	in.close();
	printf("Read File Done!\n");
}

>>>>>>> Stashed changes
