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
// -------------------------- Collision Solver--------------------------------------
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
	if (dot(prdPn, prdPn) <= 1e-6) // 如果三角形退化成一条线，就跳过这个三角形！！！
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

void CollisionSolver::CollisionResolveNew(
	BufferVector3f& fixedBuffer, BufferVector3f& vFixedBuffer, BufferVector3f& fFixedBuffer, 
	int debug, int iteration, int& debugFrameId)
{	
	m_ccdSolverTimer->Tick();
	beforeColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer; // for collision debugging
	m_vfResolveTimes.m_Data.resize(m_pbdObj->meshTopol.posBuffer.GetSize(), 0);
	//if (debugFrameId == 1)
	//{
	//	string path = "D://0324CCDTest//testData//test." + to_string(debugFrameId) + ".cache";
	//	Topology temp;
	//	temp.indices = m_pbdObj->meshTopol.indices;
	//	temp.primList = m_pbdObj->meshTopol.primList;
	//	temp.posBuffer = afterProjPrdpBuffer;
	//	temp.indices.SetName("Indices");
	//	temp.primList.SetName("primList");
	//	temp.posBuffer.SetName("P");
	//	IO::SaveToplogy(temp, path);
	//	//printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", debugFrameId);
	//}
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
		if (Contact::VF || Contact::VFDCD)
		{
			if (VFResolveTest(
				vtxPos, p1Pos, p2Pos, p3Pos, 
				m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3], 
				contact, 
				i0, i1, i2, i3))
			{
				/*
				if (debug == 50)
				{
					printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
					auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
					bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
					float3 nn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
					bool cRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nn);
					printf("distance: %f, 2*thickness: %f\n", distance, 2 * m_thickness);
					cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
					printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
						m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z,
						m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z,
						m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z,
						m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
					printf("-----------------------------------\n");
				}*/
				VFResolveNew(vtxPos, p1Pos, p2Pos, p3Pos,
					m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
					contact, i0, i1, i2, i3, fixedBuffer, vFixedBuffer, fFixedBuffer, debug, iteration);
				/*
				if (debug == 50)
				{
					vfIndices.m_Data.push_back(i0);
					vfIndices.m_Data.push_back(i1);
					vfIndices.m_Data.push_back(i2);
					vfIndices.m_Data.push_back(i3);
					string path = "D://0326Test//resoveContactTest//test." + to_string(debugFrameId) + ".cache";
					m_pbdObj->constrPBDBuffer.prdPBuffer.SetName("P");
					Topology temp;
					temp.indices = m_pbdObj->meshTopol.indices;
					temp.primList = m_pbdObj->meshTopol.primList;
					temp.posBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
					temp.indices.SetName("Indices");
					temp.primList.SetName("primList");
					temp.posBuffer.SetName("P");
					IO::SaveToplogy(temp, path); 
					debugFrameId++;

					auto newDistance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					float3 nnn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
					bool ncRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nnn);
					printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * m_thickness);
					float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
					bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
					cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;
					printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
						m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z,
						m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z,
						m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z,
						m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z); 
					printf("---------------------------------------------------------------------------------------------\n");
				}*/
			}
			else
				continue;
		}
		afterColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
		m_ccdSolverTimer->Tock();
		//PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
	}
}

// ------------- old method --------------
bool CollisionSolver::VFTest(
	float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3 vtxPrdP, float3 p1PrdP, float3 p2PrdP, float3 p3PrdP,
	Contact& contact, int i0, int i1, int i2, int i3)
{
	if (VFCCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, Contact::VF, contact, i0, i1, i2, i3))
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
	else if (VFDCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP, contact))
	{
		contact.type = Contact::VFDCD;
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
		//printf("vf thickness contact\n");
		
		return true;
	}
	else
		return false;
}

bool CollisionSolver::VFDCDTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
	float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
	Contact& contact)
{
	float d = Point2Plane(vtx_p, p1_p, p2_p, p3_p);
	if (d - 2 * m_thickness >= -1e-6)
		return false;
	float3 weight = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p);
	contact.w[0] = 1;
	contact.w[1] = weight.x;
	contact.w[2] = weight.y;
	contact.w[3] = weight.z;
	contact.n = normalize(cross(normalize(p2_p - p1_p), normalize(p3_p - p1_p))); 
	contact.t = 0.0f;
	bool inside;
	float edge1 = IO::Distance(p1_p, p2_p), edge2 = IO::Distance(p1_p, p3_p), edge3 = IO::Distance(p2_p, p3_p); // Distance funtion declear many times
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
	if (!inside)
		return false;
	return true;
}

bool CollisionSolver::VFCCDTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
	float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
	Contact::Type type, Contact& contact, int i0, int i1, int i2, int i3)
{
	const float3& x0 = vtx_o, v0 = vtx_p - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1_o - x0, x2 = p2_o - x0, x3 = p3_o - x0; // x: p's pos - V's pos
	float3 v1 = (p1_p - p1_o) - v0, v2 = (p2_p - p2_o) - v0, v3 = (p3_p - p3_o) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	for (int i = 0; i < nsol; i++)
	{
		if (t[i] < 0 || t[i] > 1)
			continue;
		contact.t = t[i];
		float3 colliPos0 = pos(vtx_o, vtx_p, t[i]), colliPos1 = pos(p1_o, p1_p, t[i]), colliPos2 = pos(p2_o, p2_p, t[i]), colliPos3 = pos(p3_o, p3_p, t[i]);
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
		float3 weight;
		contact.n = normalize(cross(p2_p - p1_p, p3_p - p1_p));
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
			float edge1 = IO::Distance(colliPos1, colliPos2), edge2 = IO::Distance(colliPos1, colliPos3), edge3 = IO::Distance(colliPos2, colliPos3); // Distance funtion declear many times
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
			//inside = (min(w[1], w[2], w[3]) >= 1e-6);
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

void CollisionSolver::VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd,
	Contact contact, int i0, int i1, int i2, int i3)
{
	float3 fn = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
	bool posRelativePos = RelativePos(vtxPos, p1Pos, fn);
	//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
	//	vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z, 
	//	vtxPrd.x, vtxPrd.y, vtxPrd.z, p1Prd.x, p1Prd.y, p1Prd.z, p2Prd.x, p2Prd.y, p2Prd.z, p3Prd.x, p3Prd.y, p3Prd.z);

	float3 n = normalize(cross(p2Prd - p1Prd, p3Prd - p1Prd));
	float depth = Point2Plane(vtxPrd, p1Prd, p2Prd, p3Prd);
	bool prdPRelativePos = RelativePos(vtxPrd, p1Prd, n);
	float dp = 0.0f;
	// 用prdp重新计算修正的weight
	float3 weight = BarycentricCoord(vtxPrd, p1Prd, p2Prd, p3Prd);
	if (posRelativePos == prdPRelativePos && depth < 2 * m_thickness) // vf dcd
	{
		dp = (2 * m_thickness - depth) * 0.5;
		//if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
			//printf("---------vf dcd----------\n");
	}
	else if (posRelativePos != prdPRelativePos) // relativePos different might be ccd
	{
		if (m_vfResolveTimes.m_Data[i0] == 0 || m_vfResolveTimes.m_Data[i1] == 0 || m_vfResolveTimes.m_Data[i2] == 0 || m_vfResolveTimes.m_Data[i3] == 0) // haven't been resolve before resove as ccd directly
		{
			dp = (depth + 2 * m_thickness) * 0.5;
			//if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
				//		printf("---------vf ccd----------\n");
		}
		else if (VFCCDTest(vtxPos, p1Pos, p2Pos, p3Pos, vtxPrd, p1Prd, p2Prd, p3Prd, Contact::VF, contact, i0, i1, i2, i3)) // have been resolve before then test again
		{
			dp = (depth + 2 * m_thickness) * 0.5;
			//if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
				//		printf("---------vf ccd----------\n");
		}		
	}
	//// for collision debugging
	//m_resolveDepths[i0].m_Data.push_back(depth);
	//m_resolveDepths[i1].m_Data.push_back(depth);
	//m_resolveDepths[i2].m_Data.push_back(depth);
	//m_resolveDepths[i3].m_Data.push_back(depth);
	if (dp != 0.0f)
	{
		m_vfResolveTimes.m_Data[i0] += 1;
		m_vfResolveTimes.m_Data[i1] += 1;
		m_vfResolveTimes.m_Data[i2] += 1;
		m_vfResolveTimes.m_Data[i3] += 1;
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
		p1Prd += -n * dp * weight.x/ sw;
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

// 1 : resolve iteration 次数
// 2 : resolve times
// 2 : condition 看
void CollisionSolver::CollisionResolve()
{
	m_ccdSolverTimer->Tick();
	beforeColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer; // for collision debugging
	m_vfResolveTimes.m_Data.resize(m_pbdObj->meshTopol.posBuffer.GetSize(), 0);
	//if (m_debugFrameID == 1)
	//{
	//	string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(m_debugFrameID) + ".cache";
	//	Topology temp;
	//	temp.indices = m_pbdObj->meshTopol.indices;
	//	temp.primList = m_pbdObj->meshTopol.primList;
	//	temp.posBuffer = afterProjPrdpBuffer;
	//	temp.indices.SetName("Indices");
	//	temp.primList.SetName("primList");
	//	temp.posBuffer.SetName("P");
	//	IO::SaveToplogy(temp, path);
	//	printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", m_debugFrameID);
	//}
	//for (int i = 0; i < 3; ++i)
	//{
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
			if (Contact::VF || Contact::VFDCD)
			{
				// 如果当前位置关系和collision-free的位置关系相同，且距离大于二倍thickness， 就不修了
				float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
				bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
				float3 nn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
				bool cRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nn);
				/*
				if (i0 == 1770 || i1 == 1770 || i2 ==1770 || i3 == 1770)
				{
					auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					if (cRelatedPos == fRelatedPos && distance >= 2 * m_thickness)
						continue;
					else
					{
						vfIndices.m_Data.push_back(i0);
						vfIndices.m_Data.push_back(i1);
						vfIndices.m_Data.push_back(i2);
						vfIndices.m_Data.push_back(i3);
						printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
						printf("distance: %f, 2*thickness: %f\n", distance, 2 * m_thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
						printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
								vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
							    m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z, 
							    m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z, 
							    m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z,
							    m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
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
						printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", m_debugFrameID);
						//printf("vtxPrd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z);
						//printf("p1Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z);
						//printf("p2Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z);
						//printf("p3Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
					}
				}
				else
				{
					auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					if (cRelatedPos == fRelatedPos && distance >= 2 * m_thickness)
						continue;
					else
					{
						VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
							m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
							contact, i0, i1, i2, i3);
					}
				}
				*/
				///*
				auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
				if (cRelatedPos == fRelatedPos && distance >= 2 * m_thickness)
					continue;
				else
				{
					/*vfIndices.m_Data.push_back(i0);
					vfIndices.m_Data.push_back(i1);
					vfIndices.m_Data.push_back(i2);
					vfIndices.m_Data.push_back(i3);*/
				/*	if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
					{
						printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
						printf("distance: %f, 2*thickness: %f\n", distance, 2 * m_thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
						printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
							vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
							m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z,
							m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z,
							m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z,
							m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
						printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
							vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
							afterProjPrdpBuffer.m_Data[i0].x, afterProjPrdpBuffer.m_Data[i0].y, afterProjPrdpBuffer.m_Data[i0].z,
							afterProjPrdpBuffer.m_Data[i1].x, afterProjPrdpBuffer.m_Data[i1].y, afterProjPrdpBuffer.m_Data[i1].z,
							afterProjPrdpBuffer.m_Data[i2].x, afterProjPrdpBuffer.m_Data[i2].y, afterProjPrdpBuffer.m_Data[i2].z,
							afterProjPrdpBuffer.m_Data[i3].x, afterProjPrdpBuffer.m_Data[i3].y, afterProjPrdpBuffer.m_Data[i3].z);
						printf("colli time %f\n", contact.t);
						printf("------------------------------------------\n");
					}*/

					VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
						m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
						contact, i0, i1, i2, i3);

					/*auto newDistance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					float3 nnn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
					bool ncRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nnn);					
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
					if (i0 == 1770 || i1 == 1770 || i2 == 1770 || i3 == 1770)
					{
						printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * m_thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;
						printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", m_debugFrameID);
					}*/
					//printf("vtxPrd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z);
					//printf("p1Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z);
					//printf("p2Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z);
					//printf("p3Prd: (%f, %f, %f)\n", m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
				}
				//*/

			}
			//if (Contact::EE)
			//{
			//	EEResolve()
			//}
		//}
	}
	afterColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
}

void CollisionSolver::CollisionResolveWithoutItreation()
{
	m_ccdSolverTimer->Tick();
	beforeColliPrdPBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer; // for collision debugging
	m_vfResolveTimes.m_Data.resize(m_pbdObj->meshTopol.posBuffer.GetSize(), 0);

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
		if (Contact::VF || Contact::VFDCD)
		{
			// 如果当前位置关系和collision-free的位置关系相同，且距离大于二倍thickness， 就不修了
			float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
			bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
			float3 nn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
			bool cRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nn);
			auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
			if (cRelatedPos == fRelatedPos && distance >= 2 * m_thickness)
				continue;
			else
			{
				VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
					m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3],
					contact, i0, i1, i2, i3);
			}

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

		if (Contact::VF || Contact::VFDCD)
		{
			// 如果当前位置关系和collision-free的位置关系相同，且距离大于二倍thickness， 就不修了
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
		//if (Contact::EE)
		//{
		//	EEResolve()
		//}
		//ColliWithShpGrd();
	}

}

void CollisionSolver::CCD_N2()
{
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
	// for t : triangles
		// for vtx : t.vertexs
			// for nbt : nbTriangle //记得去掉与当前三角形判断
				/*if( DoVF(vtx,float3 p0,float3 p1,float3 p2))
					ctxData.xxxx.push*/
					// for each edge
						// for each Triangle
							// for each edge
							// DoEE()
}

void CollisionSolver::CCD_N2Test()
{
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
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
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
	m_ccdSolverTimer->Tock();
	PBD_DEBUGTIME(m_ccdSolverTimer->GetFuncTime());
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


// ---------------- Data Oriented ----------------

void CCD_SH(
	ContactData& contactData,
	SpatialHashSystem& shs,
	Topology meshTopol, // 上一帧的collision free导出的meshtopol
	BufferVector3f prdPBuffer,
	float thickness) // 当前帧导出prdp
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
	float3 vtxPos, p1Pos, p2Pos, p3Pos, vtxPrdP, p1PrdP, p2PrdP, p3PrdP;
	int i0, i1, i2, i3;

	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		// for collision detection
		neighborList.m_Data.clear();
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		shs.FindNeighbors(neighborList, triId);

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
		//if (i0 == 50 && i1 == 316 && i2 == 67 && i3 == 276)
		//{
		//	printf("-----------------\n");
		//}	
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
	//float edge1 = IO::Distance(p1PrdP, p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
	//float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
	//bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >= (-epsilon));
	//bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >= (-epsilon));
	//bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >= (-epsilon));
	//if (inside1 && inside2 && inside3)
	//{
	//	inside = true;
	//}
	//else
	//{
	//	inside = false;
	//}
	inside = (min(contact.w[1], contact.w[2], contact.w[3]) >= 1e-6);
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
	if (i0 == 729 && i1 == 1171 && i2 == 108 && i3 == 333)
	{
		printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
		printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
		    vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z);
	}
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
			//float edge1 = IO::Distance(p1PrdP,p2PrdP), edge2 = IO::Distance(p1PrdP, p3PrdP), edge3 = IO::Distance(p2PrdP, p3PrdP); // Distance funtion declear many times
			//float epsilon = thickness / ((edge1 + edge2 + edge3) / 3.0);
			//bool inside1 = (weight.x <= (1 + epsilon)) && (weight.x >=(-epsilon));
			//bool inside2 = (weight.y <= (1 + epsilon)) && (weight.y >=(-epsilon));
			//bool inside3 = (weight.z <= (1 + epsilon)) && (weight.z >=(-epsilon));
			//if (inside1 && inside2 && inside3)
			//{
			//	inside = true;
			//}
			//else
			//{
			//	inside = false;
			//}
			inside = (min(w[1], w[2], w[3]) >= 1e-6);
			//if (i0 == 729 && (i1 == 108 || i1 == 646 || i1 == 333) && (i2 == 108 || i2 == 646 || i2 == 333) && (i3 == 108 || i3 == 646 || i3 == 333))
			if (i0 == 729 && i1 == 1171 && i2 == 108 && i3 == 333)
			{
				printf("v: %d, f: %d %d %d\n", i0, i1, i2, i3);
				printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
					vtxPrdP.x, vtxPrdP.y, vtxPrdP.z, p1PrdP.x, p1PrdP.y, p1PrdP.z, p2PrdP.x, p2PrdP.y, p2PrdP.z, p3PrdP.x, p3PrdP.y, p3PrdP.z,
					contact.t);
				printf("w: (%f, %f, %f)\n", weight.x, weight.y, weight.z);
				cout << "inside: " << inside << endl;
				//printf("contact t: %f\n", contact.t);
				printf("d :%f\n", abs(d));
			}			
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
			if (Contact::VF || Contact::VFDCD)
			{
				// 如果当前位置关系和collision-free的位置关系相同，且距离大于二倍thickness， 就不修了
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
	// 用prdp重新计算修正的weight
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
	//if (debugFrameId == 1)
	//{
	//	string path = "D://0319CCDTest//continueSimData//testData//test." + to_string(debugFrameId) + ".cache";
	//	Topology temp;
	//	temp.indices = meshTopol.indices;
	//	temp.primList = meshTopol.primList;
	//	temp.posBuffer = meshTopol.posBuffer;
	//	temp.indices.SetName("Indices");
	//	temp.primList.SetName("primList");
	//	temp.posBuffer.SetName("P");
	//	IO::SaveToplogy(temp, path);
	//	printf("------------------------------------------m_debugFrameID %d save-------------------------------------\n", debugFrameId);
	//}
	for (int i = 0; i < itereations; ++i) // change 2 to iterations
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
			if (Contact::VF || Contact::VFDCD)
			{
				if (VFResolveTest(
					vtxPos, p1Pos, p2Pos, p3Pos,
				    prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
					contact,
					i0, i1, i2, i3, thickness))
				{
				/*	if (i0 == 676 || i1 == 676 || i2 == 676 || i3 == 676)
					{
						printf("v: %d, f: %d, %d, %d\n", i0, i1, i2, i3);
						auto distance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
						float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
						bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
						float3 nn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
						bool cRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nn);
						printf("distance: %f, 2*thickness: %f\n", distance, 2 * thickness);
						cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << cRelatedPos << endl;
						printf("-----------------------------------\n");
					}*/
					VFResolveNew(vtxPos, p1Pos, p2Pos, p3Pos,
						prdPBuffer.m_Data[i0], prdPBuffer.m_Data[i1], prdPBuffer.m_Data[i2], prdPBuffer.m_Data[i3],
						contact, thickness, i0, i1, i2, i3);
					vfIndices.m_Data.push_back(i0);
					vfIndices.m_Data.push_back(i1);
					vfIndices.m_Data.push_back(i2);
					vfIndices.m_Data.push_back(i3);

					//if (i0 == 676 || i1 == 676 || i2 == 676 || i3 == 676)
					//{
					//	auto newDistance = Point2Plane(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i2], m_prdPBuffer->m_Data[i3]);
					//	float3 nnn = normalize(cross(m_prdPBuffer->m_Data[i2] - m_prdPBuffer->m_Data[i1], m_prdPBuffer->m_Data[i3] - m_prdPBuffer->m_Data[i1]));
					//	bool ncRelatedPos = RelativePos(m_prdPBuffer->m_Data[i0], m_prdPBuffer->m_Data[i1], nnn);
					//	printf("distance: %f, 2*thickness: %f\n", newDistance, 2 * m_thickness);
					//	float3 n = normalize(cross(p2Pos - p1Pos, p3Pos - p1Pos));
					//	bool fRelatedPos = RelativePos(vtxPos, p1Pos, n);
					//	cout << "fRealtedPos: " << fRelatedPos << "  " << "cRealtedPos:" << ncRelatedPos << endl;
					//	printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					//		vtxPos.x, vtxPos.y, vtxPos.z, p1Pos.x, p1Pos.y, p1Pos.z, p2Pos.x, p2Pos.y, p2Pos.z, p3Pos.x, p3Pos.y, p3Pos.z,
					//		m_prdPBuffer->m_Data[i0].x, m_prdPBuffer->m_Data[i0].y, m_prdPBuffer->m_Data[i0].z,
					//		m_prdPBuffer->m_Data[i1].x, m_prdPBuffer->m_Data[i1].y, m_prdPBuffer->m_Data[i1].z,
					//		m_prdPBuffer->m_Data[i2].x, m_prdPBuffer->m_Data[i2].y, m_prdPBuffer->m_Data[i2].z,
					//		m_prdPBuffer->m_Data[i3].x, m_prdPBuffer->m_Data[i3].y, m_prdPBuffer->m_Data[i3].z);
					//	printf("---------------------------------------------------------------------------------------------\n");
					//	/*debugFrameId++;
					//	string path = "D://0324CCDTest//testData//test." + to_string(debugFrameId) + ".cache";

					//	m_pbdObj->constrPBDBuffer.prdPBuffer.SetName("P");
					//	Topology temp;
					//	temp.indices = m_pbdObj->meshTopol.indices;
					//	temp.primList = m_pbdObj->meshTopol.primList;
					//	temp.posBuffer = m_pbdObj->constrPBDBuffer.prdPBuffer;
					//	temp.indices.SetName("Indices");
					//	temp.primList.SetName("primList");
					//	temp.posBuffer.SetName("P");
					//	IO::SaveToplogy(temp, path);*/
					//}
				}
				else
					continue;
			}
		}
	}
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
	if (dot(prdPn, prdPn) <= 1e-6) // 如果三角形退化成一条线，就跳过这个三角形！！！
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
	}
	in.close();
	printf("Read File Done!\n");
}

