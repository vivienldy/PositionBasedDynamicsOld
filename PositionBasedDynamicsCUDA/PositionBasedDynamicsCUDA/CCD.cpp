#include "CCD_Basic.h"


float3 pos(const float3 pos, const float3 prdPos, double t)
{
	return pos + t * (prdPos - pos);
}

bool IsSameSide(float3 vtx, float3 pvtx, float3 n)
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

bool CollisionSolver::VFTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
			float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,  
	        Contact& contact)
{		
	contact.type = Contact::VF;
	return(CCDCollisionTest(vtx_o, p1_o, p2_o, p3_o, vtx_p, p1_p, p2_p, p3_p, Contact::VF, contact));

	//if (CCDCollisionTest(vtx_o, p1_o, p2_o, p3_o, vtx_p, p1_p, p2_p, p3_p, Contact::VF, contact))
	//{
	//	contact.type = Contact::VF;
	//	return true;
	//}
	//else if(VFTestDCD(vtx_o, p1_o, p2_o, p3_o, vtx_p, p1_p, p2_p, p3_p, contact))
	//{
	//	contact.type = Contact::VFDCD;
	//	return true;
	//}
	//else
	//	return false;
}

bool CollisionSolver::VFTestDCD(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
	float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p,
	Contact& contact)
{
	
	if (Point2Plane(vtx_p, p1_p, p2_p, p3_p) > 2 * m_Thickness)
		return false;
	double w[4];
	w[0] = 1;
	w[1] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).x;
	w[2] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).x;
	w[3] = BarycentricCoord(vtx_p, p1_p, p2_p, p3_p).x;
	*contact.w = *w;
	float3* n  = &make_float3(0.0f,0.0f,0.0f);
	*n = cross(normalize(p2_p - p1_p), normalize(p3_p - p1_p));
	contact.n = *n;
	bool inside;
	inside = (min(-w[1], -w[2], -w[3]) >= 0);
	if (!inside)
		return false;
	return true;
}

bool CollisionSolver::CCDCollisionTest(float3 vtx_o, float3 p1_o, float3 p2_o, float3 p3_o,
	                         float3 vtx_p, float3 p1_p, float3 p2_p, float3 p3_p, 
	                         Contact::Type type, Contact& contact)
{
	const float3& x0 = vtx_o, v0 = vtx_p - x0;  // x0: V's pos   v0: V's prdPos - pos
	float3 x1 = p1_o - x0, x2 = p2_o - x0, x3 = p3_o - x0; // x: p's pos - V's pos
	float3 v1 = (p1_p - p1_o) - v0, v2 = (p2_p - p2_o) - v0, v3 = (p3_p - p3_o) - v0; // v: (p's prdPos - p's pos) - v0
	double a0 = stp(x1, x2, x3), a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3), a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3), a3 = stp(v1, v2, v3);
	double t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t); // number of solution
	t[nsol] = 1; // also check at end of timestep
	for (int i = 0; i < nsol; i++) {
		if (t[i] < 0 || t[i] > 1) 
			continue;
		contact.t = t[i];
		float3 x0 = pos(vtx_o, vtx_p, t[i]), x1 = pos(p1_o, p1_p, t[i]), x2 = pos(p2_o, p2_p, t[i]), x3 = pos(p3_o, p3_p, t[i]);
		float3& n = contact.n;
		double* w = contact.w;
		double d; // is vtx on tirangle
		bool inside;
		 //compute weight and normal
		if (type == Contact::VF) {
			d = signed_vf_distance(x0, x1, x2, x3, &n, w);
			inside = (min(-w[1], -w[2], -w[3]) >= -1e-6); 
			//inside = (min(w[1], w[2], w[3]) >= 1e-6);
		}
		else {// Impact::EE
			d = signed_ee_distance(x0, x1, x2, x3, &n, w);
			inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);
		}
		if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > 0)
			n = -n;
		if (abs(d) < 1e-6 && inside)
		{
			return true;
		}			
	}
	return false;
}

void CollisionSolver::VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
	                      float3& vtxPrd, float3& p1Prd, float3& p2Prd, float3& p3Prd, 
	                      Contact contact)
{
	cout << __FUNCTION__ << endl;
	printf("vtxPos: (%f, %f, %f)\n", vtxPos.x, vtxPos.y, vtxPos.z);
	printf("p1Pos: (%f, %f, %f)\n", p1Pos.x, p1Pos.y, p1Pos.z);
	printf("p2Pos: (%f, %f, %f)\n", p2Pos.x, p2Pos.y, p2Pos.z);
	printf("p3Pos: (%f, %f, %f))\n", p3Pos.x, p3Pos.y, p3Pos.z);
	printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
	float3 n = contact.n;
	float depth = Point2Plane(vtxPrd, p1Prd, p2Prd, p3Prd);
	printf("depth: %f\n", depth);
	bool sameSide = IsSameSide(vtxPos, p1Pos, contact.n);
	printf("normal: (%f, %f, %f)", n.x, n.y, n.z);
	cout << sameSide << endl;
	float dp;
	printf("w: (%f, %f, %f)\n", contact.w[1], contact.w[2], contact.w[3]);
	float sw = contact.w[1] * contact.w[1] + contact.w[2] * contact.w[2] + contact.w[3] * contact.w[3];
	printf("sw: %f", sw);
	if (sameSide)
	{
		if (contact.type == Contact::VF)
		{
			cout << "vf" << endl;
			dp = depth + 2 * m_Thickness;
			printf("dp: %f\n", dp);
		}
		else if(contact.type == Contact::VFDCD)
		{
			cout << "vfdcd" << endl;
			dp = 2 * m_Thickness - depth;
		}
		dp = 0.5f;
		float3 m = -n;
		vtxPrd += n * dp;
		p1Prd += -n * dp * contact.w[1] / sw;
		p2Prd += -n * dp * contact.w[2] / sw;
		p3Prd += -n * dp * contact.w[3] / sw;
	}
	else
	{
		if (contact.type == Contact::VF)
		{
			dp = depth + 2 * m_Thickness;
		}
		else if (contact.type == Contact::VFDCD)
		{
			dp = 2 * m_Thickness - depth;
		}
		vtxPrd += -n * dp;
		p1Prd += n * dp * contact.w[1] / sw;
		p2Prd += n * dp * contact.w[2] / sw;
		p3Prd += n * dp * contact.w[3] / sw;
	}
	printf("vtxPrd: (%f, %f, %f)\n", vtxPrd.x, vtxPrd.y, vtxPrd.z);
	printf("p1Prd: (%f, %f, %f)\n", p1Prd.x, p1Prd.y, p1Prd.z);
	printf("p2Prd: (%f, %f, %f)\n", p2Prd.x, p2Prd.y, p2Prd.z);
	printf("p3Prd: (%f, %f, %f)\n", p3Prd.x, p3Prd.y, p3Prd.z);
	// if 在同一侧
		// if VF
			// dp = depth + 2 * thickness;
		// if thickness穿插
			// dp = 2 * thickness - depth;
		// v prdP 沿着n方向移动 
		// p prdP 沿着n反方向移动 
	// if 不在同一侧
		// if (VF)
			// dp = depth + 2 * thickness;
		// if thickness穿插
			// dp = 2 * thickness - depth;
		// v prdP 沿着n反方向移动 
		// p prdP 沿着n方向移动 	
}

void CollisionSolver::CollisionResolve(ContactData& ctxData)
{
	// for ctx : ctxData.ctxs
		// if(Contact::VF)
				// VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
									// float3 vtxPrd, float3 p1Prd, float3 p2Prd, float3 p3Prd,
		                            // Contact ctx)
		// if(Contact::EE)
				// EEResolve()
}

void CollisionSolver::CollisionResolve()
{
	cout << __FUNCTION__ << endl;
	//for (int i = 0; i < m_Topol->posBuffer.GetSize(); ++i)
	//{
	//	printf(" (%f, %f, %f)", m_Topol->posBuffer.m_Data[i].x, m_Topol->posBuffer.m_Data[i].y, m_Topol->posBuffer.m_Data[i].z);
	//}
	Contact contact;
	int start, i0, i1, i2, i3;
	float3 vtxPos, p1Pos, p2Pos, p3Pos, vtxPrd, p1Prd, p2Prd, p3Prd;
	auto ctxIndices = &(contactData.ctxIndices);
	auto ctxStartNum = &(contactData.ctxStartNum);
	auto ctxList = &(contactData.ctxs);
	auto posBuffer = &(m_topol->posBuffer);
	for (int ctxId = 0; ctxId < ctxList->GetSize(); ++ctxId)
	{
		contact = ctxList->m_Data[ctxId];
		start = ctxStartNum->m_Data[ctxId].x;
		i0 = ctxIndices->m_Data[start];
		i1 = ctxIndices->m_Data[start+1];
		i2 = ctxIndices->m_Data[start+2];
		i3 = ctxIndices->m_Data[start+3];
		vtxPos = posBuffer->m_Data[i0];
		p1Pos = posBuffer->m_Data[i1];
		p2Pos = posBuffer->m_Data[i2];
		p3Pos = posBuffer->m_Data[i3];
		vtxPrd = m_prdPBuffer->m_Data[i0];
		p1Prd = m_prdPBuffer->m_Data[i1];
		p2Prd = m_prdPBuffer->m_Data[i2];
		p3Prd = m_prdPBuffer->m_Data[i3];
		if (Contact::VF || Contact::VFDCD)
		{
			VFResolve(vtxPos, p1Pos, p2Pos, p3Pos,
							   vtxPrd, p1Prd, p2Prd, p3Prd,
							   contact);
		}
		//if (Contact::EE)
		//{
		//	EEResolve()
		//}
	}
	// for ctx : ctxData.ctxs
		// if(Contact::VF)
				// VFResolve(float3 vtxPos, float3 p1Pos, float3 p2Pos, float3 p3Pos,
									// float3 vtxPrd, float3 p1Prd, float3 p2Prd, float3 p3Prd,
									// Contact ctx)
		// if(Contact::EE)
				// EEResolve()
}

// ArcSim VF
// EE
// collisionDetection 输出的数据结构?????
// 
//接口待确定
void CollisionSolver::CCD_N2()
{
	auto indices = &(m_topol->indices);
	auto posBuffer = &(m_topol->posBuffer);
	auto triList = &(m_topol->primList);
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
				if (VFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact))
				{
					contactData.ctxs.m_Data.push_back(contact);
					contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
					contactData.ctxIndices.m_Data.push_back(i0);
					contactData.ctxIndices.m_Data.push_back(i1);
					contactData.ctxIndices.m_Data.push_back(i2);
					contactData.ctxIndices.m_Data.push_back(i3);
					//printf("%d", ctxData.ctxIndices.GetSize());
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

void CollisionSolver::CCD_SH()
{
	auto indices = &(m_topol->indices);
	auto posBuffer = &(m_topol->posBuffer);
	auto triList = &(m_topol->primList);
	BufferInt neighborList;   // reserved for SH Find neighbor results
	float3 vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3;
	int i0, i1, i2, i3;
	//printf("--------- entering find neighbor ----------------\n");
	//m_shs->FindNeighbors(neighborList, 0);
	//printf("neighbors list: ");
	//for (int i = 0; i < neighborList.GetSize(); ++i)
		//printf("%d ", neighborList.m_Data[i]);
	for (int triId = 0; triId < triList->GetSize(); ++triId)
	{
		int start = triList->m_Data[triId].x;
		int num = triList->m_Data[triId].y;
		//printf("--------- entering find neighbor ----------------\n");
		m_shs->FindNeighbors(neighborList, triId);
		//printf("neighbors list: ");
		//for (int i = 0; i < neighborList.GetSize(); ++i)
			//printf("%d ", neighborList.m_Data[i]);
		for (int vtxId = start; vtxId < start + num; ++vtxId)
		{
			i0 = indices->m_Data[vtxId];
			vtxPos = posBuffer->m_Data[i0];
			vtxPrdP = m_prdPBuffer->m_Data[i0];
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
				triPrdP1 = m_prdPBuffer->m_Data[i1];
				triPrdP2 = m_prdPBuffer->m_Data[i2];
				triPrdP3 = m_prdPBuffer->m_Data[i3];
				Contact contact;
				if (VFTest(vtxPos, triPos1, triPos2, triPos3, vtxPrdP, triPrdP1, triPrdP2, triPrdP3, contact))
				{
					contactData.ctxs.m_Data.push_back(contact);
					contactData.ctxStartNum.m_Data.push_back(make_int2(contactData.ctxIndices.GetSize(), 4));
					contactData.ctxIndices.m_Data.push_back(i0);
					contactData.ctxIndices.m_Data.push_back(i1);
					contactData.ctxIndices.m_Data.push_back(i2);
					contactData.ctxIndices.m_Data.push_back(i3);
					//printf("%d", ctxData.ctxIndices.GetSize());
				}
			}
		}
	}
}

void CCDTestMain()
{
	CCDTest ccdTest;
	ccdTest.PrepareTestData("D://0310ContinuousCollisionDectection//ccdTestData//ccdClothRestData.txt",
											 "D://0310ContinuousCollisionDectection//ccdTestData//ccdClothData.txt");
	CollisionSolver colliSolver;

	float3 cellSize = make_float3(5.10f, 5.10f, 5.10f);
	float3 gridCenter = make_float3(0.0f, 0.0f, 0.0f);
	uint3 gridSize = make_uint3(3, 1, 3);

	// initialize SH
	SpatialHashSystem shs(ccdTest.topol.posBuffer, ccdTest.topol.indices, CPU);
	shs.SetGridCenter(gridCenter);
	shs.SetGridSize(gridSize);
	shs.SetDivision(cellSize);
	shs.InitSH();
	shs.UpdateSH(0.0f);

	colliSolver.SetTarget(&ccdTest.topol, &ccdTest.prdPBuffer);
	colliSolver.SetThickness(0.5f);
	colliSolver.SetAcceStruct(&shs);

	auto start = chrono::steady_clock::now();
	clock_t tStart = clock();

	// collision detection
	colliSolver.CCD_N2();
	//colliSolver.CCD_SH();
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;

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
	// collision resolves
	//colliSolver.CollisionResolve();
	//ccdTest.writePointsToFile();
}



double CollisionSolver::signed_vf_distance(const float3& x,
	const float3& y0, const float3& y1, const float3& y2,
	float3* n, double* w)
{
	//cout << __FUNCTION__ << endl;
	float3 _n;
	if (!n)
		n = &_n;
	double _w[4];
	if (!w)
		w = _w;
	*n = cross(normalize(y1 - y0), normalize(y2 - y0));
	//printf("vectors (%f, %f,  %f)\n", (y1 - y0).x, (y1 - y0).y, (y1 - y0).z);
	//printf("vectors (%f, %f,  %f)\n", (y2 - y0).x, (y2 - y0).y, (y2 - y0).z);
	if (norm2(*n) < 1e-6)
		return infinity;
	*n = normalize(*n);
	//printf("n (%f, %f,  %f)\n", n->x, n->y, n->z);
	double h = dot(x - y0, *n);
	double b0 = stp(y1 - x, y2 - x, *n),
		b1 = stp(y2 - x, y0 - x, *n),
		b2 = stp(y0 - x, y1 - x, *n);
	w[0] = 1;
	w[1] = -b0 / (b0 + b1 + b2);
	w[2] = -b1 / (b0 + b1 + b2);
	w[3] = -b2 / (b0 + b1 + b2);
	//printf("w (%f, %f,  %f)\n", w[1], w[2], w[3]);
	/*float3 weight = BarycentricCoord(x, y0, y1, y2);
	w[1] = weight.x;
	w[2] = weight.y;
	w[3] = weight.z;*/
	//printf("w (%f, %f,  %f)\n", BarycentricCoord(x, y0, y1, y2).x, BarycentricCoord(x, y0, y1, y2).y, BarycentricCoord(x, y0, y1, y2).z);
	
	return h;
}

double CollisionSolver::signed_ee_distance(const float3& x0, const float3& x1,
	const float3& y0, const float3& y1,
	float3* n, double* w)
{
	float3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	if (norm2(*n) < 1e-6)
		return infinity;
	*n = normalize(*n);
	double h = dot(x0 - y0, *n);
	double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
		b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	w[0] = a0 / (a0 + a1);
	w[1] = a1 / (a0 + a1);
	w[2] = -b0 / (b0 + b1);
	w[3] = -b1 / (b0 + b1);
	return h;
}
//----------------------- CCDTest ----------------------------

void CCDTest::PrepareTestData(string topolPath, string prdPath)
{
	readMeshFromTxt(topolPath);
	readBufferFromTxt(prdPath);
	createPrimList(3);
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


// DoVF(dt)
/*{
	prdBuffer = posBuffer + t * v;    //calculate four points position after time t
	解三次方程，取最小非负实根
	float e = thickness/edgeAverage;    // calculate epsilon
	if( t < dt )
		将t带入方程求解w1, w2
		if( -e < w1 < 1+e && -e < w2 < 1+e)
			w3 = 1 - w1 - w2
			if( -e < w1 < 1+e)
				colliPrimList.push_back(make_int2(colliIndices.GetSize(), 4))
				colliIndices.push_back(当前点idx，以及当前collide三角形三个点idx)
				colliType.push_back(VF);
}*/

// DoEE(tolerance)
/*{
	计算两条边向量
	if(两条边平行)
		continue;
	解方程求a和b，计算两个交点
	if(a<1 && b < 1) //两个点都在线段上
		d = Distance(a1, a2)
		if( d < tolerance )
			colliPrimList.push_back(make_int2(colliIndices.GetSize(), 4))
			colliIndices.push_back(四个点的idx)
			colliType.push_back(EE);
	else if()
}*/




	