/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include "MyFunc.h"

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
	/* for you to implement ... */
    Vector finalForce[8][8][8]={0};

    int i,j,k,ip,jp,kp;
    Vector oneDirectF;
	double structuralRL = JelloFixedR;
	double shearRL = sqrt(2) * JelloFixedR;
	double shearFRL = sqrt(3) * JelloFixedR;
	double bendRL = 2 * JelloFixedR;

	for (i=0; i<=7; i++)
		for (j=0; j<=7; j++)
			for (k=0; k<=7; k++)
			{
				// Structural
				NEIGHBOR_FORCE(1,0,0,structuralRL);
				NEIGHBOR_FORCE(0,1,0,structuralRL);
				NEIGHBOR_FORCE(0,0,1,structuralRL);
				NEIGHBOR_FORCE(-1,0,0,structuralRL);
				NEIGHBOR_FORCE(0,-1,0,structuralRL);
				NEIGHBOR_FORCE(0,0,-1,structuralRL);

				// Shear
				NEIGHBOR_FORCE(1,1,0,shearRL);
				NEIGHBOR_FORCE(-1,1,0,shearRL);
				NEIGHBOR_FORCE(-1,-1,0,shearRL);
				NEIGHBOR_FORCE(1,-1,0,shearRL);
				NEIGHBOR_FORCE(0,1,1,shearRL);
				NEIGHBOR_FORCE(0,-1,1,shearRL);
				NEIGHBOR_FORCE(0,-1,-1,shearRL);
				NEIGHBOR_FORCE(0,1,-1,shearRL);
				NEIGHBOR_FORCE(1,0,1,shearRL);
				NEIGHBOR_FORCE(-1,0,1,shearRL);
				NEIGHBOR_FORCE(-1,0,-1,shearRL);
				NEIGHBOR_FORCE(1,0,-1,shearRL);

				NEIGHBOR_FORCE(1,1,1,shearFRL);
				NEIGHBOR_FORCE(-1,1,1,shearFRL);
				NEIGHBOR_FORCE(-1,-1,1,shearFRL);
				NEIGHBOR_FORCE(1,-1,1,shearFRL);
				NEIGHBOR_FORCE(1,1,-1,shearFRL);
				NEIGHBOR_FORCE(-1,1,-1,shearFRL);
				NEIGHBOR_FORCE(-1,-1,-1,shearFRL);
				NEIGHBOR_FORCE(1,-1,-1,shearFRL);

				// Bend
				NEIGHBOR_FORCE(2,0,0,bendRL);
				NEIGHBOR_FORCE(0,2,0,bendRL);
				NEIGHBOR_FORCE(0,0,2,bendRL);
				NEIGHBOR_FORCE(-2,0,0,bendRL);
				NEIGHBOR_FORCE(0,-2,0,bendRL);
				NEIGHBOR_FORCE(0,0,-2,bendRL);

				// Compute Collision and Response
				Vector penaltyForce;
				penaltyForce = PenaltyForce(jello, i, j, k);
				pSUM(finalForce[i][j][k], penaltyForce, finalForce[i][j][k]);

				// Compute External Force
				if (IsExternalForce(jello))
				{
					Vector externalForce;
					externalForce = ExternalForce(jello, i, j, k);
					pSUM(finalForce[i][j][k], externalForce, finalForce[i][j][k]);
				}

				// Compute Inclined Collision and Response
				if (jello->incPlanePresent)
				{
					Vector inclineForce;
					inclineForce = InclinedForce( jello, i, j, k);
					pSUM( finalForce[i][j][k], inclineForce, finalForce[i][j][k]);
				}

				// Compute MouseForce
				if (isMouseForce)
				{
					Vector mouseForce;
					mouseForce = MouseForce(jello);
					pSUM(finalForce[i][j][k], mouseForce, finalForce[i][j][k]);
				}

				a[i][j][k].x=finalForce[i][j][k].x/jello->mass;
				a[i][j][k].y=finalForce[i][j][k].y/jello->mass;
				a[i][j][k].z=finalForce[i][j][k].z/jello->mass;
			}
	if (isMouseForce)
		isMouseForce = false;
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
	int i,j,k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i=0; i<=7; i++)
		for (j=0; j<=7; j++)
			for (k=0; k<=7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8],
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }
  return;
}

