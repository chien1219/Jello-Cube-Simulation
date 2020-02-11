#pragma once
#include <stdio.h>
#include "jello.h"
#include <cmath>        // std::abs

extern point msForceUp, msForceHa;
extern double msforceScale;

double ComputeLength(Vector dest);
double DotProduct(Vector A, Vector B);
Vector Normalize(Vector dest);
bool CollisionDetection(struct point pos, double boxsize, int *result);
bool IsOtherForce(struct world * jello);
Vector PenaltyForce(struct world * jello, int i, int j, int k);
Vector InclinedForce(struct world * jello, int i, int j, int k);
Vector ExternalForce(struct world * jello, int i, int j, int k);
Vector ForceInSpring(struct world * jello, double RLength, struct point pos1, struct point pos2, struct point V1, struct point V2);
Vector MouseForce(struct world * jello);