#include "MyFunc.h"

double ComputeLength(Vector dest)
{
	double length;
	length = sqrt((dest).x * (dest).x + (dest).y * (dest).y + (dest).z * (dest).z);
	return length;
}

double DotProduct(Vector A, Vector B)
{
	double result = A.x*B.x + A.y*B.y + A.z*B.z;
	return result;
}

Vector Normalize(Vector src)
{
	double length;
	Vector dest;
	length = sqrt((src).x * (src).x + (src).y * (src).y + (src).z * (src).z);

	(dest).x = src.x / length;
	(dest).y = src.y / length;
	(dest).z = src.z / length;
	return dest;
}

Vector ForceInSpring(struct world * jello, double RestLength, struct point pos1, struct point pos2, struct point V1, struct point V2)
{
	// Hook Force
	double distanceL, RLength;
	Vector L, NorL;
	Vector FHook;
	RLength = RestLength;

	pDIFFERENCE(pos1, pos2, L);
	distanceL = ComputeLength(L);
	NorL = Normalize(L);
	FHook.x = -1 * jello->kElastic * (distanceL - RLength) * NorL.x;
	FHook.y = -1 * jello->kElastic * (distanceL - RLength) * NorL.y;
	FHook.z = -1 * jello->kElastic * (distanceL - RLength) * NorL.z;

	// Damping Froce
	Vector V1_V2;
	Vector FDamping;
	pDIFFERENCE(V1, V2, V1_V2);
	DotProduct(V1_V2, L);
	FDamping.x = -1 * jello->dElastic*(DotProduct(V1_V2, L) / distanceL)* NorL.x;
	FDamping.y = -1 * jello->dElastic*(DotProduct(V1_V2, L) / distanceL)* NorL.y;
	FDamping.z = -1 * jello->dElastic*(DotProduct(V1_V2, L) / distanceL)* NorL.z;

	Vector FinalForce;
	pSUM(FHook, FDamping, FinalForce);

	return FinalForce;
}

Vector InclinedForce(struct world * jello, int i, int j, int k)
{
	double state = jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y + jello->c * jello->p[i][j][k].z + jello->d;
	if (state * VerticesState[i][j][k] <= 0)
	{
		Vector L;
		Vector PlaneNormal;
		Vector NorL;
		double distance = std::abs(state) / (sqrt(jello->a * jello->a + jello->b * jello->b + jello->c * jello->c));
		PlaneNormal.x = -jello->a;
		PlaneNormal.y = -jello->b;
		PlaneNormal.z = -jello->c;
		NorL = Normalize(PlaneNormal);
		L.x = NorL.x * distance;
		L.y = NorL.y * distance;
		L.z = NorL.z * distance;

		// Hook Force
		Vector FHook;
		FHook.x = -1 * jello->kCollision * distance * NorL.x;
		FHook.y = -1 * jello->kCollision * distance * NorL.y;
		FHook.z = -1 * jello->kCollision * distance * NorL.z;

		// Damping Force
		Vector FDamping;
		Vector Va;
		Va.x = jello->v[i][j][k].x;
		Va.y = jello->v[i][j][k].y;
		Va.z = jello->v[i][j][k].z;

		if (distance == 0)
			distance = 1;

		FDamping.x = -1 * jello->dCollision * (DotProduct(Va, L)) / distance * NorL.x;
		FDamping.y = -1 * jello->dCollision * (DotProduct(Va, L)) / distance * NorL.y;
		FDamping.z = -1 * jello->dCollision * (DotProduct(Va, L)) / distance * NorL.z;

		Vector FinalInclinedForce;
		pSUM(FHook, FDamping, FinalInclinedForce);
		return FinalInclinedForce;
	}
	else
	{
		Vector zero;
		zero.x = 0;
		zero.y = 0;
		zero.z = 0;
		return zero;
	}
}

bool CollisionDetection(struct point pos, double boxsize, int *result)
{
	int i, j, k;

	if (pos.x > boxsize)
		*result = 1;
	else if (pos.x < -1 * boxsize)
		*result = -1;
	else if (pos.y > boxsize)
		*result = 2;
	else if (pos.y < -1 * boxsize)
		*result = -2;
	else if (pos.z > boxsize)
		*result = 3;
	else if (pos.z < -1 * boxsize)
		*result = -3;
	else
		return false;

	return true;
}

Vector PenaltyForce(struct world * jello, int i, int j, int k)
{
	int condition = 0;

	if (CollisionDetection(jello->p[i][j][k], 2.0, &condition))
	{
		Vector FHook, L;
		Vector NorL;
		switch (condition)
		{
		case 1:
			L.x = jello->p[i][j][k].x - 2;
			L.y = 0;
			L.z = 0;
			break;
		case -1:
			L.x = jello->p[i][j][k].x - (-2);
			L.y = 0;
			L.z = 0;
			break;
		case 2:
			L.x = 0;
			L.y = jello->p[i][j][k].y - 2;
			L.z = 0;
			break;
		case -2:
			L.x = 0;
			L.y = jello->p[i][j][k].y - (-2);
			L.z = 0;
			break;
		case 3:
			L.x = 0;
			L.y = 0;
			L.z = jello->p[i][j][k].z - 2;
			break;
		case -3:
			L.x = 0;
			L.y = 0;
			L.z = jello->p[i][j][k].z - (-2);
			break;
		}

		// Hook Force
		NorL = Normalize(L);
		FHook.x = -1 * jello->kCollision * (ComputeLength(L) - 0) * NorL.x;
		FHook.y = -1 * jello->kCollision * (ComputeLength(L) - 0) * NorL.y;
		FHook.z = -1 * jello->kCollision * (ComputeLength(L) - 0) * NorL.z;

		// Damping Force
		Vector FDamping;
		Vector Va;
		Va.x = jello->v[i][j][k].x;
		Va.y = jello->v[i][j][k].y;
		Va.z = jello->v[i][j][k].z;

		FDamping.x = -1 * jello->dCollision * (DotProduct(Va, L)) / ComputeLength(L) * NorL.x;
		FDamping.y = -1 * jello->dCollision * (DotProduct(Va, L)) / ComputeLength(L) * NorL.y;
		FDamping.z = -1 * jello->dCollision * (DotProduct(Va, L)) / ComputeLength(L) * NorL.z;

		Vector FinalCollisionForce;

		pSUM(FHook, FDamping, FinalCollisionForce);

		return FinalCollisionForce;
	}
	else
	{
		Vector zero;
		zero.x = 0;
		zero.y = 0;
		zero.z = 0;
		return zero;
	}
}

bool IsOtherForce(struct world * jello)
{
	if (jello->resolution > 0 && EXforce)
		return true;
	else
		return false;
}

Vector ExternalForce(struct world * jello, int i, int j, int k)
{
	point pos;
	//change to force field coordinate
	pos.x = ((jello->p[i][j][k].x + 2) / 4.0)*(jello->resolution - 1); 
	pos.y = ((jello->p[i][j][k].y + 2) / 4.0)*(jello->resolution - 1);
	pos.z = ((jello->p[i][j][k].z + 2) / 4.0)*(jello->resolution - 1);

	// boundary check
	if (pos.x > (jello->resolution - 1)) { pos.x = (jello->resolution - 1); }
	else if (pos.x < 0) { pos.x = 0; }
	if (pos.y > (jello->resolution - 1)) { pos.y = (jello->resolution - 1); }
	else if (pos.y < 0) { pos.y = 0; }
	if (pos.z > (jello->resolution - 1)) { pos.z = (jello->resolution - 1); }
	else if (pos.z < 0) { pos.z = 0; }

	int xf, xc, yf, yc, zf, zc;  // f = floor c = ceil
	double s, t, g; //distance to each axis

	xf = floor(pos.x); xc = ceil(pos.x);
	yf = floor(pos.y); yc = ceil(pos.y);
	zf = floor(pos.z); zc = ceil(pos.z);

	s = pos.x - xf; //x
	t = pos.y - yf; //y
	g = pos.z - zf; //z

	Vector ExternalF;

	ExternalF.x = (1 - s)*(1 - t)*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zf].x
		+ s * (1 - t)*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zf].x
		+ (1 - s)*t*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zf].x
		+ s * t*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zf].x
		+ s * (1 - t)*g*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zc].x
		+ (1 - s)*t*g*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zc].x
		+ (1 - s)*(1 - t)*g*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zc].x
		+ s * t*g*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zc].x;
	ExternalF.y = (1 - s)*(1 - t)*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zf].y
		+ s * (1 - t)*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zf].y
		+ (1 - s)*t*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zf].y
		+ s * t*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zf].y
		+ s * (1 - t)*g*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zc].y
		+ (1 - s)*t*g*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zc].y
		+ (1 - s)*(1 - t)*g*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zc].y
		+ s * t*g*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zc].y;

	ExternalF.z = (1 - s)*(1 - t)*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zf].z
		+ s * (1 - t)*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zf].z
		+ (1 - s)*t*(1 - g)*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zf].z
		+ s * t*(1 - g)*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zf].z
		+ s * (1 - t)*g*jello->forceField[xc * jello->resolution * jello->resolution + yf * jello->resolution + zc].z
		+ (1 - s)*t*g*jello->forceField[xf * jello->resolution * jello->resolution + yc * jello->resolution + zc].z
		+ (1 - s)*(1 - t)*g*jello->forceField[xf * jello->resolution * jello->resolution + yf * jello->resolution + zc].z
		+ s * t*g*jello->forceField[xc * jello->resolution * jello->resolution + yc * jello->resolution + zc].z;

	return ExternalF;
}

Vector MouseForce(struct world * jello)
{
	double distanceL;
	double magnitude = 0.8;
	Vector L, NorL;
	Vector MouseF;

	point pos1, pos2;
	pos1.x = 0; pos2.x = 0;
	pos1.y = UpPos[0]; pos2.y = DownPos[0];
	pos1.z = -UpPos[1]; pos2.z = -DownPos[1];

	pDIFFERENCE(pos1, pos2, L);

	distanceL = ComputeLength(L);
	if (distanceL != 0)
		NorL = Normalize(L);
	else
	{
		NorL.x = 0;
		NorL.y = 0;
		NorL.z = 0;
	}

	MouseF.x = magnitude * distanceL * NorL.x;
	MouseF.y = magnitude * distanceL * NorL.y;
	MouseF.z = magnitude * distanceL * NorL.z;

	return MouseF;
}