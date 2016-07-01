/*
 * Physics.c
 *
 *  Created on: Feb 24, 2016
 *      Author: abauville
 */

#include "stokes.h"



void Physics_allocateMemory(Physics* Physics, Grid* Grid)
{


	int i;
	Physics->Vx 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
	Physics->P 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->eta 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->etaVisc		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->eta0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->n 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->rho 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );



#if (HEAT)
	Physics->k 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->DT 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
#endif
	Physics->psi 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dpsi 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->kD				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->SD				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );


	//Physics->psi 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	//Physics->Dpsi 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->G 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->cohesion 		= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->frictionAngle 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );





	// Initialize stuff
	//int i;
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0;
	}
	for (i = 0; i < Grid->nECTot; ++i) {
		//Physics->P[i]  = 0;

#if (HEAT)
		Physics->T[i]  = 0;
		Physics->DT[i] = 0;
#endif

		Physics->P[i] = 0;

		//Physics->eta[i] = 0;
		//Physics->rho[i] = 0;

		Physics->psi[i]  = 0;
		Physics->Dpsi[i] = 0;

		Physics->sigma_xx_0[i] = 0;
	}
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->sigma_xy_0[i] = 0;
	}










}


void Physics_freeMemory(Physics* Physics)
{
	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P );

	free( Physics->eta );
	//free( Physics->etaVisc );
	free( Physics->eta0 );

	free(Physics->etaVisc);

	free( Physics->n );
	free( Physics->rho );





#if (HEAT)
	free( Physics->k );
	free(Physics->T );
	free(Physics->DT );
#endif



	free(Physics->G );


	free(Physics->sigma_xx_0 );
	free(Physics->sigma_xy_0 );
	free(Physics->Dsigma_xx_0 );
	free(Physics->Dsigma_xy_0 );

	free(Physics->cohesion);
	free(Physics->frictionAngle);

	// Darcy
	free(Physics->psi );
	free(Physics->Dpsi );
	free(Physics->kD);
	free(Physics->SD);






}


void Physics_initPToLithostatic(Physics* Physics, Grid* Grid)
{
	int ix, iy, iCell;
	//printf("=== P ===\n");
	compute rho_g_h;
	// Set Temp to zero (the interpolation forced the ghost values to have a dirichlet value follow the dirichlet)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->T[iCell] = 0;
	}

	// Initialize the pressure at the lithostatic pressure
	// =========================

		// Initialize P at the lithostatic pressure
		// Contribution of y
		if (Physics->g[0]>0) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				rho_g_h = 0;
				for (iy = 0; iy < Grid->nyEC; --iy) {
					iCell = ix + iy*Grid->nxEC;
					rho_g_h += Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->dy;

					Physics->P[iCell] = 1*rho_g_h;


				}
				//printf("\n");
			}
		} else {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				rho_g_h = 0;
				for (iy = Grid->nyEC-1; iy >= 0; --iy) {
					iCell = ix + iy*Grid->nxEC;
					rho_g_h += Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->dy;

					Physics->P[iCell] = 1*rho_g_h;


				}
				//printf("\n");
			}
		}


		// Contribution of x // in case the gravity field is not vertical
		// be careful adding contribution from left to right. This assumes the model is diping right.
		// If the model dips in the other direction the loop should be from right to left
		if (Physics->g[0]>=0) {
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				rho_g_h = 0;
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					rho_g_h += Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->dx;
					//printf("%.2e  ", Physics->P[iCell]);
					Physics->P[iCell] += 1*rho_g_h;

				}
				//printf("\n");
			}
		} else {
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				rho_g_h = 0;
				for (ix = Grid->nxEC-1; ix >=0; --ix) {
					iCell = ix + iy*Grid->nxEC;
					rho_g_h += Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->dx;
					//printf("%.2e  ", Physics->P[iCell]);
					Physics->P[iCell] += 1*rho_g_h;

				}
				//printf("\n");
			}
		}
}



void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal, BC* BCThermal)
{

	// Declarations
	// =========================
	int iCell, i;
	int nNeighbours = 4;
	coord locX, locY;

	coord dx = Grid->dx;
	coord dy = Grid->dy;
	compute* sumOfWeights 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* eta0 			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* n    			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* rho  			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* G    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* cohesion 		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* frictionAngle 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* sigma_xx_0   	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* sigma_xy_0   	= (compute*) malloc(nNeighbours * Grid->nSTot * sizeof(compute));

	compute* psi   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* kD   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* SD   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#if (HEAT)
	compute* T    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* k    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#endif


	// Reinitialize Physics array
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		eta0[i] = 0;
		n[i] = 0;
		rho[i] = 0;


		G  [i] = 0;
		sigma_xx_0 [i] = 0;
		cohesion[i] = 0;
		frictionAngle[i] = 0;

		psi[i] = 0;
		kD[i] = 0;
		SD[i] = 0;

		sumOfWeights[i] = 0;
		#if (HEAT)
		k  [i] = 0;
		T  [i] = 0;
		#endif
	}

	for (i = 0; i < nNeighbours * Grid->nSTot; ++i) {
		sigma_xy_0 [i] = 0;
	}





	printf("Init ok\n");


	compute weight;


	//int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	//int IyMod[4] = {0,0,1,1};

	int phase;

	int nxEC = Grid->nxEC;
	compute xMod[4], yMod[4];
	int ix,  iy;

	int I, C;

	xMod[0] = -1; yMod[0] = -1;
	xMod[1] =  1; yMod[1] = -1;
	xMod[2] = -1; yMod[2] =  1;
	xMod[3] =  1; yMod[3] =  1;



	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right




	SingleParticle* thisParticle = NULL;
	// Loop through inner cells
	// ========================
	int iNode = 0;
	//int Count = 0;
	//#pragma omp parallel for private(ix, iNode, thisParticle, locX, locY, phase, Ix, Iy, i, locEta, iNode, dVxdy, dVydx, dVxdx, dVydy, locEps_II) schedule(static,32)







	//printf("=== Part Temp ===\n");
	// Loop through inner nodes
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 0; ix < Grid->nxS; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				phase = thisParticle->phase;
				//printf("phase = %i, locX= %.2f, locY=%.2f  \n", thisParticle->phase,locX,locY);
				// Loop through neighbours
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);

					weight = fabs((locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5));


					eta0			[iCell*4+i] += MatProps->eta0[phase] * weight;
					n				[iCell*4+i] += MatProps->n   [phase] * weight;


					G				[iCell*4+i] += 1/MatProps->G [phase] * weight; // harmonic average
					cohesion		[iCell*4+i] += MatProps->cohesion[phase] * weight;
					frictionAngle	[iCell*4+i] += MatProps->frictionAngle[phase] * weight;
#if (HEAT)
					rho				[iCell*4+i] += MatProps->rho0[phase] * weight;// * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
					k				[iCell*4+i] += MatProps->k   [phase] * weight;
					T 				[iCell*4+i] += thisParticle->T * weight;
#else
					rho				[iCell*4+i] += MatProps->rho0[phase]*weight;//* (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
#endif

					sigma_xx_0 		[iCell*4+i] += thisParticle->sigma_xx_0 * weight;
					sumOfWeights	[iCell*4+i] += weight;

					psi				[iCell*4+i] += thisParticle->psi * weight;

					SD				[iCell*4+i] += MatProps->SD  [phase] * weight;

					/*
					if (thisParticle->faulted == true) {
						kD				[iCell*4+i] += FAULT_MOD*MatProps->kD  [phase] * weight;
					} else {
						kD				[iCell*4+i] += MatProps->kD  [phase] * weight;
					}
					*/

				}
				thisParticle = thisParticle->next;
			}
		}
	}



	printf("Left Right contrib\n");
	// Add contribution from the other side in the case of periodic BC
	if(BCStokes->SetupType==SimpleShearPeriodic) {
		int IyNPeriod[2],IxNPeriod[2];
		IyNPeriod[0] = 0;
		IyNPeriod[1] = 1;
		compute xModPeriod[2], yModPeriod[2];
		yModPeriod[0] = -1.0;
		yModPeriod[1] =  1.0;

		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ix+=Grid->nxS-1) {
				if (ix == 0) {
					IxNPeriod[0]  = Grid->nxS;
					IxNPeriod[1]  = Grid->nxS;
					xModPeriod[0] = -1.0;
					xModPeriod[1] = -1.0;
				} else if (ix==Grid->nxS-1) {
					IxNPeriod[0]  = -(Grid->nxS-1);
					IxNPeriod[1]  = -(Grid->nxS-1);
					xModPeriod[0] =  1.0;
					xModPeriod[1] =  1.0;
				}

				//printf("ix = %i, ixN = %i, nxS = %i, nxEC = %i \n",ix, ix+IxNPeriod[0], Grid->nxS, Grid->nxEC);
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];

				// Loop through the particles in the shifted cell
				// ======================================
				while (thisParticle!=NULL) {

					locX = (thisParticle->x-Grid->xmin)/dx - ix;
					locY = (thisParticle->y-Grid->ymin)/dy - iy;


					phase = thisParticle->phase;
					//printf("phase = %i, locX= %.2f, locY=%.2f  \n", thisParticle->phase,locX,locY);
					// Loop through neighbours
					for (i=0; i<2; i++) {
						iCell = (ix+IxNPeriod[i] + (iy+IyNPeriod[i]) * nxEC);
						weight = fabs((locX + xModPeriod[i]*0.5)   *   (locY + yModPeriod[i]*0.5));


						eta0			[iCell*4+i] += MatProps->eta0[phase] * weight;
						n				[iCell*4+i] += MatProps->n   [phase] * weight;


						G				[iCell*4+i] += 1/MatProps->G [phase] * weight; // harmonic average
						cohesion		[iCell*4+i] += MatProps->cohesion[phase] * weight;
						frictionAngle	[iCell*4+i] += MatProps->frictionAngle[phase] * weight;
#if (HEAT)
						rho				[iCell*4+i] += MatProps->rho0[phase] * weight;// * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
						k				[iCell*4+i] += MatProps->k   [phase] * weight;
						T 				[iCell*4+i] += thisParticle->T * weight;
#else
						rho				[iCell*4+i] += MatProps->rho0[phase]*weight;//* (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
#endif

						sigma_xx_0 		[iCell*4+i] += thisParticle->sigma_xx_0 * weight;
						sumOfWeights	[iCell*4+i] += weight;

						psi				[iCell*4+i] += thisParticle->psi * weight;

						SD				[iCell*4+i] += MatProps->SD  [phase] * weight;



					}
					thisParticle = thisParticle->next;
				}
			}
		}
	}




	printf("Left Right contrib end\n");





	//printf("=== eta fill ===\n");

	compute sum;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			I = 4*iCell;
			sum = sumOfWeights[I+0] + sumOfWeights[I+1] + sumOfWeights[I+2] + sumOfWeights[I+3];
			//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
			if (sum==0) {
				printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
				exit(0);
			}
#if (HEAT)
			Physics->T  [iCell] =     (   T[I+0] +    T[I+1] +    T[I+2] +    T[I+3]) / sum;
			Physics->k  [iCell] =     (   k[I+0] +    k[I+1] +    k[I+2] +    k[I+3]) / sum;
#endif
			Physics->eta0[iCell]=     (eta0[I+0] + eta0[I+1] + eta0[I+2] + eta0[I+3]) / sum ; // harmonic average
			Physics->n  [iCell] =     (   n[I+0] +    n[I+1] +    n[I+2] +    n[I+3]) / sum;
			Physics->rho[iCell] =     ( rho[I+0] +  rho[I+1] +  rho[I+2] +  rho[I+3]) / sum;

			Physics->G  [iCell] = sum/(   G[I+0] +    G[I+1] +    G[I+2] +    G[I+3]) ; // harmonic average
			Physics->cohesion   [iCell] = ( cohesion  [I+0] + cohesion  [I+1] + cohesion  [I+2] + cohesion  [I+3])/sum;
			Physics->frictionAngle[iCell] = ( frictionAngle[I+0] + frictionAngle[I+1] + frictionAngle[I+2] + frictionAngle[I+3]) / sum;
			Physics->sigma_xx_0 [iCell] = ( sigma_xx_0[I+0] + sigma_xx_0[I+1] + sigma_xx_0[I+2] + sigma_xx_0[I+3]) / sum ; // harmonic average

			Physics->psi [iCell] =   ( psi[I+0] +  psi[I+1] +  psi[I+2] +  psi[I+3]) / sum;
			Physics->kD  [iCell] =   ( kD [I+0] +  kD [I+1] +  kD [I+2] +  kD [I+3]) / sum;
			Physics->SD  [iCell] =   ( SD [I+0] +  SD [I+1] +  SD [I+2] +  SD [I+3]) / sum;

			//printf("Physics->simga_xx_0[%i] %.2e, sigma_xx_0 = %.2f %.2f %.2f %.2f, sum = %.2f\n", iCell, Physics->sigma_xx_0[iCell], sigma_xx_0[I+0], sigma_xx_0[I+1], sigma_xx_0[I+2], sigma_xx_0[I+3], sum);


		}
	}
	//printf("=== end fill ===\n");




	// Replace boundary values by their neighbours
	int INeigh;
#if (HEAT)
	int IBC;
#endif
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (BCStokes->SetupType==SimpleShearPeriodic) {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		} else {
			if (ix==0) {
				INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
			} else if (ix==Grid->nxEC-1) {
				INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
			} else {
				INeigh =   ix + (iy+1)*Grid->nxEC  ;
			}
		}
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];

		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

		Physics->psi[I] = Physics->psi[INeigh]+Grid->dy;
		Physics->kD [I] = Physics->kD [INeigh];
		Physics->SD [I] = Physics->SD [INeigh];

#if (HEAT)
		Physics->k  [I] = Physics->k  [INeigh];
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}
		}
#endif

	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;


		if (BCStokes->SetupType==SimpleShearPeriodic) {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		} else {
			if (ix==0) {
				INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
			} else if (ix==Grid->nxEC-1) {
				INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
			} else {
				INeigh =   ix + (iy-1)*Grid->nxEC  ;
			}
		}
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];

		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

		Physics->psi[I] = Physics->psi[INeigh];
		Physics->kD [I] = Physics->kD [INeigh];
		Physics->SD [I] = Physics->SD [INeigh];
#if (HEAT)
		Physics->k  [I] = Physics->k  [INeigh];
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}

		}
#endif


	}





	if (BCStokes->SetupType!=SimpleShearPeriodic) {

		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {

			I = ix + iy*Grid->nxEC;

			INeigh =   ix+1 + (iy)*Grid->nxEC  ;
			Physics->eta0[I] = Physics->eta0[INeigh];
			Physics->n   [I] = Physics->n   [INeigh];
			Physics->rho[I] = Physics->rho[INeigh];

			Physics->G  [I] = Physics->G  [INeigh];
			Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
			Physics->cohesion     [I] = Physics->cohesion     [INeigh];
			Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

			Physics->psi[I] = Physics->psi[INeigh];
			Physics->kD [I] = Physics->kD [INeigh];
			Physics->SD [I] = Physics->SD [INeigh];

#if (HEAT)

			Physics->k  [I] = Physics->k  [INeigh];
			IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
			if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
				Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
			}
			else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
				if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
					Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
				}
				if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
					Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
				}

			}


#endif

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			Physics->eta0[I] = Physics->eta0[INeigh];
			Physics->n   [I] = Physics->n   [INeigh];
			Physics->rho[I] = Physics->rho[INeigh];

			Physics->G  [I] = Physics->G  [INeigh];
			Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
			Physics->cohesion     [I] = Physics->cohesion     [INeigh];
			Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

			Physics->psi[I] = Physics->psi[INeigh];
			Physics->kD [I] = Physics->kD [INeigh];
			Physics->SD [I] = Physics->SD [INeigh];
#if (HEAT)
			if (BCThermal->SetupType!=SimpleShearPeriodic) {
				Physics->k  [I] = Physics->k  [INeigh];
				IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n

				if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
					Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
				}
				else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
					if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
						Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
					}
					if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
						Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
					}
				}
			}
#endif
		}
	}
	printf("end neighbour stuff");



	// ==================================
	// Interpolate to nodes
	// ==================================

	// Reinitialize sum of weights
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		sumOfWeights[i] = 0;
	}

	int signX, signY, iNodeNeigh;
	int Counter;
	xMod[0] =  1; yMod[0] =  1;
	xMod[1] =  0; yMod[1] =  1;
	xMod[2] =  1; yMod[2] =  0;
	xMod[3] =  0; yMod[3] =  0;
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight, signX, signY, iNodeNeigh) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 0; ix < Grid->nxS; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			Counter = 0;

			while (thisParticle!=NULL) {
				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				phase = thisParticle->phase;

				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}




				for (i=0; i<4; i++) {
					iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

					if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
						printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
						printf("IX = %i, IY = %i, locX = %.3f, locY = %.3f, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
						printf("thisParticle->x = %.3f , y = %.3f \n", thisParticle->x, thisParticle->y);
						exit(0);
					}

					locX = fabs(locX);
					locY = fabs(locY);


					//printf("iNodeNeigh = %i, signX = %i, signY = %i\n", iNodeNeigh, signX, signY);
					weight = (locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5);

					sigma_xy_0 	[iNodeNeigh*4+i] += thisParticle->sigma_xy_0 * weight;
					sumOfWeights[iNodeNeigh*4+i] += weight; // using the same arrays



				}
				Counter++;
				thisParticle = thisParticle->next;
			}
		}
	}

	printf("end first loop for sigma_xy\n");

	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		I = 4*iNode;
		sum = sumOfWeights[I+0] + sumOfWeights[I+1] + sumOfWeights[I+2] + sumOfWeights[I+3];
		//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		Physics->sigma_xy_0 [iNode] = ( sigma_xy_0[I+0] +  sigma_xy_0[I+1] +  sigma_xy_0[I+2] +  sigma_xy_0[I+3]) / sum ; // harmonic average


	}

	printf("end filling loop for sigma_xy");









	if (DEBUG) {

		printf("=== Check eta0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta0[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== Check rho 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->rho[C]);
				C++;
			}
			printf("\n");
		}

#if (HEAT)
		printf("=== Check k 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->k[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check T 1 ===\n");
				C = 0;
				//int ix, iy;
				for (iy = 0; iy < Grid->nyEC; ++iy) {
					for (ix = 0; ix < Grid->nxEC; ++ix) {
						printf("%.3f  ", Physics->T[C]);
						C++;
					}
					printf("\n");
				}
#endif
		printf("=== Check G 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->G[C]);
				C++;
			}
			printf("\n");
		}



		printf("=== Check sigma_xx 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->sigma_xx_0[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check sigma_xy 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.3e  ", Physics->sigma_xy_0[C]);
				C++;
			}
			printf("\n");
		}

	}

	free(sumOfWeights);
	free(eta0);
	free(n);
	free(rho);


	free(G);
	free(sigma_xx_0);
	free(sigma_xy_0);
	free(cohesion);
	free(frictionAngle);

	free(psi);
	free(kD );
	free(SD );

#if (HEAT)
	free(T);
	free(k);
#endif



}


















#if (HEAT)
void Physics_interpTempFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;




	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;


				thisParticle->T  += ( .25*(1.0-locX)*(1.0-locY)*Physics->DT[ix  +(iy  )*Grid->nxEC]
									+ .25*(1.0-locX)*(1.0+locY)*Physics->DT[ix  +(iy+1)*Grid->nxEC]
									+ .25*(1.0+locX)*(1.0+locY)*Physics->DT[ix+1+(iy+1)*Grid->nxEC]
									+ .25*(1.0+locX)*(1.0-locY)*Physics->DT[ix+1+(iy  )*Grid->nxEC] );

				thisParticle = thisParticle->next;
			}
		}
	}


}
#endif


void Physics_interpPsiFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;





	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;



				thisParticle->psi+= ( .25*(1.0-locX)*(1.0-locY)*Physics->Dpsi[ix  +(iy  )*Grid->nxEC]
																			  + .25*(1.0-locX)*(1.0+locY)*Physics->Dpsi[ix  +(iy+1)*Grid->nxEC]
																														+ .25*(1.0+locX)*(1.0+locY)*Physics->Dpsi[ix+1+(iy+1)*Grid->nxEC]
																																								  + .25*(1.0+locX)*(1.0-locY)*Physics->Dpsi[ix+1+(iy  )*Grid->nxEC] );


				thisParticle = thisParticle->next;
			}
		}
	}


}








void Physics_interpStressesFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal, MatProps* MatProps)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	int signX, signY;


	// See Taras' book pp. 186-187
	compute sigma_xx_0_fromNodes;
	compute sigma_xy_0_fromNodes;

	compute d_ve = 0.95;
	compute dtm = Physics->dt;
	compute dtMaxwell;

	compute eta, G;

	compute* Dsigma_xy_sub_OnTheNodes = (compute*) malloc( 4*Grid->nSTot *sizeof(compute) );
	compute* sumOfWeights_OnTheNodes = (compute*) malloc( 4*Grid->nSTot *sizeof(compute) );
	compute* Dsigma_xx_sub_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* sumOfWeights_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );


	compute* Dsigma_xy_rem_OnTheNodes = (compute*) malloc( Grid->nSTot *sizeof(compute) );
	compute* Dsigma_xx_rem_OnTheCells = (compute*) malloc( Grid->nECTot *sizeof(compute) );

	int i;
	for (i=0;i<4*Grid->nSTot;++i) {
		Dsigma_xy_sub_OnTheNodes[i] = 0;
		sumOfWeights_OnTheNodes[i] = 0;
	}
	for (i=0;i<4*Grid->nECTot;++i) {
		Dsigma_xx_sub_OnTheCells[i] = 0;
		sumOfWeights_OnTheCells[i] = 0;
	}

	int iNodeNeigh;
	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right


	int iCell;
	int xModNode[4], yModNode[4], xModCell[4], yModCell[4];
	compute weight;
	xModNode[0] =  1; yModNode[0] =  1;
	xModNode[1] =  0; yModNode[1] =  1;
	xModNode[2] =  1; yModNode[2] =  0;
	xModNode[3] =  0; yModNode[3] =  0;


	xModCell[0] = -1; yModCell[0] = -1;
	xModCell[1] =  1; yModCell[1] = -1;
	xModCell[2] = -1; yModCell[2] =  1;
	xModCell[3] =  1; yModCell[3] =  1;

	compute Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart;


	// compute Dsigma_xx_0_sub on the particles and interpolate to the grid
	#pragma omp parallel for private(iy, ix, i, iNode, thisParticle, locX, locY, signX, signY, sigma_xx_0_fromNodes, sigma_xy_0_fromNodes, eta, G, dtMaxwell, Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart, iNodeNeigh, weight, iCell) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			//printf("iNode = %i\n",iNode);


			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;


				sigma_xx_0_fromNodes  = ( .25*(1.0-locX)*(1.0-locY)*Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC]
										+ .25*(1.0-locX)*(1.0+locY)*Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0+locY)*Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0-locY)*Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC] );


				eta  				  = ( .25*(1.0-locX)*(1.0-locY)*Physics->eta[ix  +(iy  )*Grid->nxEC]
										+ .25*(1.0-locX)*(1.0+locY)*Physics->eta[ix  +(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0+locY)*Physics->eta[ix+1+(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0-locY)*Physics->eta[ix+1+(iy  )*Grid->nxEC] );


				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from
				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}


				locX = fabs(locX)-1;
				locY = fabs(locY)-1;


				sigma_xy_0_fromNodes = ( .25*(1.0-locX)*(1.0-locY)*Physics->sigma_xy_0[ix      +(iy  )    *Grid->nxS]
											  + .25*(1.0-locX)*(1.0+locY)*Physics->sigma_xy_0[ix      +(iy+signY)*Grid->nxS]
											   + .25*(1.0+locX)*(1.0+locY)*Physics->sigma_xy_0[ix+signX+(iy+signY)*Grid->nxS]
												+ .25*(1.0+locX)*(1.0-locY)*Physics->sigma_xy_0[ix+signX+(iy  )    *Grid->nxS] );






				G = MatProps->G[thisParticle->phase];

				dtMaxwell = eta/G;

				// Compute Dsigma sub grid
				Dsigma_xx_sub_OnThisPart = ( sigma_xx_0_fromNodes - thisParticle->sigma_xx_0 ) * ( 1 - exp(-d_ve * dtm/dtMaxwell) );
				Dsigma_xy_sub_OnThisPart = ( sigma_xy_0_fromNodes - thisParticle->sigma_xy_0 ) * ( 1 - exp(-d_ve * dtm/dtMaxwell) );


				//printf("dve = %.2e, Dsigma_xx_sub_OnThisPart = %.2e, Dsigma_xy_sub_OnThisPart  = %.2e\n", d_ve, Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart );

				// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
				thisParticle->sigma_xx_0 += Dsigma_xx_sub_OnThisPart;
				thisParticle->sigma_xy_0 += Dsigma_xy_sub_OnThisPart;


				// redefine locX, locY (used to compute surface based weight, not used as weight directly)
				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;



				// Interp Dsigma_xy_sub from particles to Nodes
				for (i=0; i<4; i++) {

					iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

					if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
						printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
						//printf("IX = %i, IY = %i, locX = %.3f, locY = %.3f, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
						printf("thisParticle->x = %.3f , y = %.3f \n", thisParticle->x, thisParticle->y);
						exit(0);
					}


					locX = fabs(locX);
					locY = fabs(locY);


					//printf("iNodeNeigh = %i, signX = %i, signY = %i\n", iNodeNeigh, signX, signY);
					weight = (locX + xModNode[i]*0.5)   *   (locY + yModNode[i]*0.5);

					//printf("A, iNodeNeigh =%i\n", iNodeNeigh);
					Dsigma_xy_sub_OnTheNodes[iNodeNeigh*4+i] += Dsigma_xy_sub_OnThisPart * weight;
					sumOfWeights_OnTheNodes [iNodeNeigh*4+i] += weight; // using the same arrays



				}



				// Interp Dsigma_xx_sub from particles to Cells
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
					weight = fabs((locX + xModCell[i]*0.5)   *   (locY + yModCell[i]*0.5));

					Dsigma_xx_sub_OnTheCells[iCell*4+i] += Dsigma_xx_sub_OnThisPart * weight;
					sumOfWeights_OnTheCells	[iCell*4+i] += weight;

				}





				//printf("dve = %.2e, Dsigma_xy_sub_OnThisPart            = %.2e\n", d_ve, Dsigma_xy_sub_OnThisPart );
				//printf("dve = %.2e, Dsigma_xy_sub_OnTheNodes[iNode*4+0] = %.2e\n", d_ve, Dsigma_xy_sub_OnTheNodes[iNode*4+0] );

				thisParticle = thisParticle->next;

			}

			// Interpolate Dsigma_xy_sub to the nodes


		}
	}

	printf("A\n");

	int I;
	compute sum;
	compute Dsigma_xy_sub_OnThisNode;
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		I = 4*iNode;
		sum = sumOfWeights_OnTheNodes[I+0] + sumOfWeights_OnTheNodes[I+1] + sumOfWeights_OnTheNodes[I+2] + sumOfWeights_OnTheNodes[I+3];
		//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
		if (sum==0) {
			printf("error in Physics_interpStressesFromCellsToParticles: node #%i received no contribution from particles\n", iNode );
			exit(0);
		}

		Dsigma_xy_sub_OnThisNode = ( Dsigma_xy_sub_OnTheNodes[I+0] +  Dsigma_xy_sub_OnTheNodes[I+1] +  Dsigma_xy_sub_OnTheNodes[I+2] +  Dsigma_xy_sub_OnTheNodes[I+3]) / sum ; // harmonic average
		//printf("dve = %.2e, Dsigma_xy_sub_OnTheNodes[iCell*4+i] = %.2e\n", d_ve, Dsigma_xy_sub_OnTheNodes[iNode*4+i] );
		//printf("Dsigma_xy_sub_OnThisNode = %.2e\n", Dsigma_xy_sub_OnThisNode);

		Dsigma_xy_rem_OnTheNodes[iNode] = Physics->Dsigma_xy_0[iNode] - Dsigma_xy_sub_OnThisNode;


	}



	compute Dsigma_xx_sub_OnThisCell;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			I = 4*iCell;
			sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];
			//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
			if (sum==0) {
				printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
				exit(0);
			}

			Dsigma_xx_sub_OnThisCell = ( Dsigma_xx_sub_OnTheCells[I+0] + Dsigma_xx_sub_OnTheCells[I+1] + Dsigma_xx_sub_OnTheCells[I+2] + Dsigma_xx_sub_OnTheCells[I+3]) / sum ; // harmonic average
			Dsigma_xx_rem_OnTheCells[iCell] = Physics->Dsigma_xx_0[iCell] - Dsigma_xx_sub_OnThisCell;

			//printf("Physics->simga_xx_0[%i] %.2e, sigma_xx_0 = %.2f %.2f %.2f %.2f, sum = %.2f\n", iCell, Physics->sigma_xx_0[iCell], sigma_xx_0[I+0], sigma_xx_0[I+1], sigma_xx_0[I+2], sigma_xx_0[I+3], sum);


		}
	}

	// Replace boundary values by their neighbours
	int INeigh;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}

		Dsigma_xx_rem_OnTheCells[I] = Dsigma_xx_rem_OnTheCells[INeigh];

	}

	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}

		Dsigma_xx_rem_OnTheCells[I] = Dsigma_xx_rem_OnTheCells[INeigh];

	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;

		Dsigma_xx_rem_OnTheCells[I] = Dsigma_xx_rem_OnTheCells[INeigh];

	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;

		Dsigma_xx_rem_OnTheCells[I] = Dsigma_xx_rem_OnTheCells[INeigh];

	}







	// interpolate Dsigma sub back to the nodes and cell centers

















	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, signX, signY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;

				//compute locX0 = locX;
				//compute locY0 = locY;


				thisParticle->sigma_xx_0  += ( .25*(1.0-locX)*(1.0-locY)*Dsigma_xx_rem_OnTheCells[ix  +(iy  )*Grid->nxEC]
											+ .25*(1.0-locX)*(1.0+locY)*Dsigma_xx_rem_OnTheCells[ix  +(iy+1)*Grid->nxEC]
											+ .25*(1.0+locX)*(1.0+locY)*Dsigma_xx_rem_OnTheCells[ix+1+(iy+1)*Grid->nxEC]
 											+ .25*(1.0+locX)*(1.0-locY)*Dsigma_xx_rem_OnTheCells[ix+1+(iy  )*Grid->nxEC] );


				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from


				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}


				locX = fabs(locX)-1;
				locY = fabs(locY)-1;


				thisParticle->sigma_xy_0  += ( .25*(1.0-locX)*(1.0-locY)*Dsigma_xy_rem_OnTheNodes[ix      +(iy  )    *Grid->nxS]
											 + .25*(1.0-locX)*(1.0+locY)*Dsigma_xy_rem_OnTheNodes[ix      +(iy+signY)*Grid->nxS]
											 + .25*(1.0+locX)*(1.0+locY)*Dsigma_xy_rem_OnTheNodes[ix+signX+(iy+signY)*Grid->nxS]
											 + .25*(1.0+locX)*(1.0-locY)*Dsigma_xy_rem_OnTheNodes[ix+signX+(iy  )    *Grid->nxS] );


				thisParticle = thisParticle->next;
			}
		}
	}






	free(Dsigma_xy_sub_OnTheNodes);
	free(Dsigma_xx_sub_OnTheCells);

	free(sumOfWeights_OnTheNodes);
	free(sumOfWeights_OnTheCells);


	free(Dsigma_xy_rem_OnTheNodes);
	free(Dsigma_xx_rem_OnTheCells);




}










void Physics_get_VxVyP_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int InoDir, INeigh;
	compute maxVx = 0;
	compute maxVy = 0;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = -1;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = -1;
	}
	for (i = 0; i < Grid->nCTot; ++i) {
		Physics->P[i] = -1;
	}


	// Set Vx
	// =========================
	C = 0;
	int IBC;
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			I = ix + iy*Grid->nxVx;
			InoDir = Numbering->map[I];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vx[C] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1; // BC nodes are numbered -1 to -n
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vx[C] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (iy==0)  // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVx  ];
					if (iy==Grid->nyVx-1)  // lower boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVx  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vx[C] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (iy==0)  // lower boundary
							Physics->Vx[C] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dy;
						if (iy==Grid->nyVx-1)  // lower boundary
							Physics->Vx[C] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dy;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}

			// Get maxVx
			if (Physics->Vx[C]*Physics->Vx[C] > maxVx)
				maxVx = Physics->Vx[C]*Physics->Vx[C];


			C++;
		}
	}

	// Set Vy
	// =========================

	C = 0;
	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			I = ix + iy*Grid->nxVy + Grid->nVxTot;

			InoDir = Numbering->map[I];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vy[C] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vy[C] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (ix==0)  // lower boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVy + Grid->nVxTot ];
					if (ix==Grid->nxVy-1)  // lower boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVy + Grid->nVxTot  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vy[C] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (ix==0)  // lower boundary
							Physics->Vy[C] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dx;
						if (ix==Grid->nxVy-1)  // lower boundary
							Physics->Vy[C] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dx;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}

			if (Physics->Vy[C]*Physics->Vy[C] > maxVy)
				maxVy = Physics->Vy[C]*Physics->Vy[C];
			C++;
		}
	}


	int Iref;
	compute RefPressure;
	// Shift pressure, taking the pressure of the upper right node as reference (i.e. 0)
	Iref = Grid->nCTot-Grid->nxC + Grid->nVxTot + Grid->nVyTot;
	InoDir = Numbering->map[Iref];
	if (InoDir>=0) { // Not a Dirichlet node
		RefPressure = EqSystem->x[InoDir];
	} else {
		IBC = abs(InoDir)-1;
		RefPressure = BC->value[ IBC ];
	}
	int IE; // Index of embedded nodes
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			IE = ix + iy*Grid->nxEC; // Index embedded
			I = (ix-1) + (iy-1)*Grid->nxC + Grid->nVxTot + Grid->nVyTot; // Index in NumStokes
			InoDir = Numbering->map[I];
			if (InoDir>=0) { // Not a Dirichlet node
				Physics->P[IE] = EqSystem->x[InoDir];
			} else {
				IBC = abs(InoDir)-1;
				Physics->P[IE] = BC->value[ IBC ];
			}

			Physics->P[IE] -= RefPressure;
		}
	}
	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix-1 + (iy)*Grid->nxEC];
	}





	Physics->maxV = sqrt(maxVx+maxVy);





	if (DEBUG) {
		// Check Vx
		// =========================
		printf("=== Vx ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyVx; ++iy) {
			for (ix = 0; ix < Grid->nxVx; ++ix) {
				printf("%.3e  ", Physics->Vx[C]);
				C++;
			}
			printf("\n");
		}


		// Check Vy
		// =========================
		C = 0;
		printf("=== Vy ===\n");
		for (iy = 0; iy < Grid->nyVy; ++iy) {
			for (ix = 0; ix < Grid->nxVy; ++ix) {
				printf("%.3e  ", Physics->Vy[C]);
				C++;
			}
			printf("\n");
		}

		// Check P
		// =========================
		printf("=== P ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->P[C]);
				C++;
			}
			printf("\n");
		}
	}






}















#if (HEAT)
void Physics_get_T_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumThermal, EqSystem* EqThermal)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int INeigh, IBC;
	//int InoDir, INeigh;
	//compute maxVx = 0;
	//compute maxVy = 0;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================




	compute* Told = malloc(Grid->nECTot * sizeof(compute));
	for (i = 0; i < Grid->nECTot; ++i) {
		Told[i] = Physics->T[i];
	}


	printf("== Check T filling ==\n");
	C = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		for (ix = 0; ix<Grid->nxEC; ix++) {
			I = NumThermal->map[C];
			if (I>=0) {
				Physics->T[C] = EqThermal->x[I];
				Physics->DT[C] = Physics->T[C]  - Told[C];
			}
			else {

				IBC = abs(I)-1; // BC nodes are numbered -1 to -n

				// Get neighbours index
				if (iy==0) { // lower boundary
					if (ix==0) {
						INeigh = NumThermal->map[  ix+1 + (iy+1)*Grid->nxEC  ];
					} else if (ix==Grid->nxEC-1) {
						INeigh = NumThermal->map[  ix-1 + (iy+1)*Grid->nxEC  ];
					} else {
						INeigh = NumThermal->map[  ix + (iy+1)*Grid->nxEC  ];
					}

				} else if (iy==Grid->nyEC-1)  { //  upper boundary
					if (ix==0) {
						INeigh = NumThermal->map[  ix+1 + (iy-1)*Grid->nxEC  ];
					} else if (ix==Grid->nxEC-1) {
						INeigh = NumThermal->map[  ix-1 + (iy-1)*Grid->nxEC  ];
					} else {
						INeigh = NumThermal->map[  ix + (iy-1)*Grid->nxEC  ];
					}
				} else if (ix==0) { // left boundary
					INeigh = NumThermal->map[  ix+1 + (iy)*Grid->nxEC  ];
				} else if (ix==Grid->nxEC-1) { // right boundary
					INeigh = NumThermal->map[  ix-1 + (iy)*Grid->nxEC  ];
				}





				if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Physics->T[C] = 2.0*BC->value[IBC] - EqThermal->x[INeigh];
					Physics->DT[C] = Physics->T[C] - Told[C];
				}
				else if (BC->type[IBC]==NeumannGhost) { // Neumann
					if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
						Physics->T[C] = EqThermal->x[INeigh] - BC->value[IBC]*Grid->dx;
						Physics->DT[C] = Physics->T[C] - Told[C];
					}
					if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
						Physics->T[C] = EqThermal->x[INeigh] + BC->value[IBC]*Grid->dy;
						Physics->DT[C] = Physics->T[C] - Told[C];
					}
				}
				else {
					printf("error: unknown boundary type\n");
					exit(0);
				}
				//printf("C=%i, IBC=%i, Type=%i, value=%.3f, valueNeigh=%.3f, FinalValue=%.3f\n",C, IBC,BC->type[IBC], BC->value[IBC], EqThermal->x[INeigh], Physics->T[C]);


				//Physics->T[C] = BC->value[abs(I)];
			}
			C++;
		}
	}




	if (DEBUG) {
		// Check T
		// =========================
		printf("=== T ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->T[C]);
				C++;
			}
			printf("\n");
		}
	}

	free(Told);



}
#endif





void Physics_computeStressChanges(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes)
{

	// see Taras' book p. 186
	int ix, iy, iCell, iNode;
	compute Z;
	compute Eps_xx, Eps_xy;
	compute dVxdy, dVydx;
	compute GShear, etaShear;
	// compute stress
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;
			Eps_xx 	=  (Physics->Vx[ix + iy*Grid->nxVx] - Physics->Vx[ix-1 + iy*Grid->nxVx])/Grid->dx;


			Z 		= (Physics->G[iCell]*Physics->dt)  /  (Physics->eta[iCell] + Physics->G[iCell]*Physics->dt);

			Physics->Dsigma_xx_0[iCell] = ( 2*Physics->eta[iCell] * Eps_xx  -  Physics->sigma_xx_0[iCell] ) * Z;

		}
	}



	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix-1 + (iy)*Grid->nxEC];
	}






	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]
								  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]
								  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Eps_xy = 0.5*(dVxdy+dVydx);


			GShear 	 	= shearValue(Physics->G, ix, iy, Grid->nxEC);
			etaShear 	= shearValue(Physics->eta, ix, iy, Grid->nxEC);

			Z 			= (GShear*Physics->dt)  /  (etaShear + GShear*Physics->dt);

			Physics->Dsigma_xy_0[iNode] = ( 2*etaShear * Eps_xy   -   Physics->sigma_xy_0[iNode] ) * Z;

		}
	}






}








void Physics_computeStrainRateInvariant(Physics* Physics, Grid* Grid, compute* StrainRateInvariant)
{
	// Definition of second invariant: // E_II = sqrt( Eps_xx^2 + Eps_xy^2  );
	// Declarations
	// =========================
	int ix, iy, IE;
	compute dVxdy, dVydx, dVxdx, dVydy;
//	int iCell;

	// ix, iy modifiers
	//int iNode, Ix, Iy;
	//int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	//int IyMod[4] = {0,0,1,1};

//	int IxModC[4] = {0,1,-1,0, 0}; // lower left, lower right, upper right, upper left
//	int IyModC[4] = {0,0, 0,1,-1};

#pragma omp parallel for private(ix,iy, IE, dVxdx, dVydy, dVxdy, dVydx) schedule(static,32)
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			IE = ix+iy*Grid->nxEC;
			//I = (ix-1)+(iy-1)*Grid->nxC;


			/*
			// Method A: using the averageing of derivatives on the four nodes
			// Compute Eps_xy at the four nodes of the cell
			// 1. Sum contributions
			dVxdy = 0;
			dVydx = 0;
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];

				dVxdy += ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
						 - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


				dVydx += ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
						 - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;

			}
			// 2. Average
			dVxdy /= 4;
			dVydx /= 4;
			*/



			// Method A simplifies to:


			dVxdy = ( Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy-1)*Grid->nxVx] +
					  Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy-1)*Grid->nxVx] )/4./Grid->dy;


			dVydx = ( Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] - Physics->Vy[(ix-1)+(iy-1)*Grid->nxVy] +
					  Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix-1)+(iy  )*Grid->nxVy] )/4./Grid->dx;



			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
					   - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
				   - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			StrainRateInvariant[IE] = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  0.5*dVxdx*dVxdx  +  0.5*dVydy*dVydy);

		}
	}

	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix-1 + (iy)*Grid->nxEC];
	}
}

void Physics_computeStrainInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* EII)
{
	compute dVxdy, dVydx, dVxdx, dVydy;
	//int iNode, Ix, Iy;
	//int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	//int IyMod[4] = {0,0,1,1};


	dVxdy = ( Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy-1)*Grid->nxVx] +
			  Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy-1)*Grid->nxVx] )/4./Grid->dy;


	dVydx = ( Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] - Physics->Vy[(ix-1)+(iy-1)*Grid->nxVy] +
			  Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix-1)+(iy  )*Grid->nxVy] )/4./Grid->dx;


	dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

	dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;


	*EII = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  0.5*dVxdx*dVxdx  +  0.5*dVydy*dVydy);

}



void Physics_computeEta(Physics* Physics, Grid* Grid, Numerics* Numerics)
{
	int iCell, iy, ix;
	compute sigma_y, sigmaII;
	//compute* EIIGrid = (compute*) malloc(Grid->nECTot*sizeof(compute));
	//Physics_computeStrainRateInvariant(Physics, Grid, EIIGrid);



	compute EII_visc, eta_visc, EII;
	compute sigma_xx, sigma_xy;

	compute alpha, sigma_xxT;


	int C = 0;

	printf("timeStep = %i, itNonLin = %i\n", Numerics->timeStep, Numerics->itNonLin);


#pragma omp parallel for private(ix,iy, iCell, sigma_xy, sigma_xx, sigmaII, sigma_y, EII_visc, EII) schedule(static,32)
	//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			//printf("ix = %i, iy = %i, iCell = %i, Grid->nxEC = %i\n", ix, iy, iCell, Grid->nxEC);

			//  Compute new stresses:
			// xy from node to cell center
			sigma_xy  = Physics->sigma_xy_0[ix-1 + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix   + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix-1 + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix   + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
			sigma_xy /= 4.0;

			sigma_xx = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];

			// Rotation correction (might not be needed since no advection is not performed)
			/*
			alpha = - 0.5*Physics->dt*((Physics->Vy[ix+1+iy*Grid->nxVy]   - Physics->Vy[ix+(iy)*Grid->nxVy])/Grid->dx
					- (Physics->Vx[ix+(iy+1)*Grid->nxVx] - Physics->Vx[ix+(iy)*Grid->nxVx])/Grid->dy);
			sigma_xxT = sigma_xx* cos(alpha) * cos(alpha) - sigma_xy * sin(2*alpha);
			sigma_xy = sigma_xy * cos(2*alpha)  		 +  sigma_xx * sin(2*alpha);
			sigma_xx = sigma_xxT;
			*/

			if (Numerics->timeStep<=0 && Numerics->itNonLin <= 0){
				Physics->etaVisc[iCell] = Physics->eta0[iCell];
				Physics->eta[iCell] = Physics->etaVisc[iCell];

			} else {


				sigmaII = sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);
				EII_visc = sigmaII/(2*Physics->etaVisc[iCell]);

				// Compute powerlaw rheology
				Physics->etaVisc[iCell] = Physics->eta0[iCell] * pow(EII_visc/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);
				Physics->eta[iCell] = Physics->etaVisc[iCell];


				// Plasticity
				sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Physics->P[iCell] * sin(Physics->frictionAngle[iCell]);
				if (sigmaII>sigma_y) {
					Physics_computeStrainInvariantForOneCell(Physics, Grid, ix,iy, &EII);
					Physics->eta[iCell] = sigma_y / (2*EII);
					sigmaII = sigma_y;
				}




				if (Physics->eta[iCell]<Numerics->etaMin) {
					Physics->eta[iCell] = Numerics->etaMin;
				}
				else if (Physics->eta[iCell]>Numerics->etaMax) {
					Physics->eta[iCell] = Numerics->etaMax;
				}
			}

		}
	}





	/*
	//#pragma omp parallel for private(iCell, maxCorr, EII, EII_visc, eta_visc, sigma_y, sigmaII, corr, C) schedule(static,32)
	//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			//printf("ix = %i, iy = %i, iCell = %i, Grid->nxEC = %i\n", ix, iy, iCell, Grid->nxEC);

			// Initial guess


			alpha = - 0.5*Physics->dt*((Physics->Vy[ix+1+iy*Grid->nxVy]   - Physics->Vy[ix+(iy)*Grid->nxVy])/Grid->dx
					- (Physics->Vx[ix+(iy+1)*Grid->nxVx] - Physics->Vx[ix+(iy)*Grid->nxVx])/Grid->dy);

			//  Compute new stresses:
			sigma_xx = Physics->sigma_xx_0[iCell] * cos(alpha) * cos(alpha) - Physics->sigma_xy_0[iCell] * sin(2*alpha);
			sigma_xy = Physics->sigma_xy_0[iCell] * cos(2*alpha)  		 +  Physics->sigma_xx_0[iCell] * sin(2*alpha);

			if (Numerics->timeStep==0 && Numerics->itNonLin==0) {
				eta_visc = Physics->eta0[iCell] * pow(EII/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);
				//sigmaII = 2*Physics->eta0[iCell]*EII;
			} else {
				EII_visc = EII;
				sigmaII = sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);


			//sigmaII = 2*Physics->eta0[iCell]*EII;

			maxCorr = 1E20;

			EII = EIIGrid[iCell];





			C = 0;

			sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Physics->P[iCell] * sin(Physics->frictionAngle[iCell]);

			while (maxCorr>tol) {


			// Compute powerlaw rheology
			eta_visc = Physics->eta0[iCell] * pow(EII_visc/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);

			Physics->eta[iCell] = eta_visc;
			// Compute the yield stress

			//sigmaII = 2*Physics->eta0[iCell]*EII;

			if (sigmaII>sigma_y) {
				Physics->eta[iCell] = sigma_y / (2*EII);
				sigmaII = sigma_y;
			}


			corr = sigmaII/(2*eta_visc) - EII_visc;

			EII_visc += 0.9*corr;

			maxCorr = corr/EII;
			C++;
			}

			if (ix == 5 && iy == 5) {
				printf("******nonLinear local it: n = %.1f, C = %i, maxCorr = %.2e, corr = %.2e, sigmaII = %.2e\n", Physics->n[iCell], C, maxCorr, corr, sigmaII);
				printf("******sigma_xx = %.2e, sigma_xx_0 = %.2e\n",sigma_xx, Physics->sigma_xx_0[iCell]);
			}




			if (Physics->eta[iCell]<Numerics->etaMin) {
				Physics->eta[iCell] = Numerics->etaMin;
			}
			else if (Physics->eta[iCell]>Numerics->etaMax) {
				Physics->eta[iCell] = Numerics->etaMax;
			}
			}

		}
	}
	 */

	// Replace boundary values by their neighbours
	int INeigh, I;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Physics->eta[I] = Physics->eta[INeigh];
	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}
		Physics->eta[I] = Physics->eta[INeigh];
	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {

		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Physics->eta[I] = Physics->eta[INeigh];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Physics->eta[I] = Physics->eta[INeigh];

	}
















	if (DEBUG) {
		printf("===== Compute Eta =====\n");
		printf("=== Check eta0  ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta0[C]);
				C++;
			}
			printf("\n");
		}

		/*
		printf("=== Check  ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", EIIGrid[C]);
				C++;
			}
			printf("\n");
		}
		*/

		printf("=== Check eta ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
		}


		printf("=== Check n ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->n[C]);
				C++;
			}
			printf("\n");
		}
	}


	//free(EIIGrid);

}



void Physics_changePhaseOfFaults(Physics* Physics, Grid* Grid, MatProps* MatProps, Particles* Particles)
{

	compute* EII = (compute*) malloc(Grid->nECTot*sizeof(compute));
	Physics_computeStrainRateInvariant(Physics, Grid, EII);

	SingleParticle* thisParticle = NULL;

	int ix, iy, iNode;
	compute EII_node;


	for (iy = 0; iy < Grid->nyS; ++iy) {

		for (ix = 0; ix < Grid->nxS; ++ix) {


			iNode = ix+ iy*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			EII_node = (EII[ix+iy*Grid->nxEC] + EII[ix+1+iy*Grid->nxEC] + EII[ix+(iy+1)*Grid->nxEC] + EII[ix+1+(iy+1)*Grid->nxEC])/4;

			if (EII_node>10*Physics->epsRef) {

				while (thisParticle != NULL) {

					thisParticle->faulted = true;
					thisParticle = thisParticle->next;
				}
			}
		}
	}


	free(EII);

}



void Physics_updateDt(Physics* Physics, Numerics* Numerics)
{

	if (fabs(Physics->maxV)<1E-6)
		Physics->maxV = 1E-6;
	Physics->dt = Numerics->CFL_fac*Numerics->dLmin/(Physics->maxV); // note: the min(dx,dy) is the char length, so = 1
	//printf("maxV = %.3em, Physics.dt = %.3e, Physics.dt(SCALED)= %.3e yr, dtmin = %.2e, dtmax = %.2e, dtMax = %.2e\n",fabs(Physics.maxV), Physics.dt, Physics.dt*Char.time/3600/24/365, dtmin, dtmax, dtMax);


	if (Physics->dt<Numerics->dtmin) {
		Physics->dt = Numerics->dtmin;
	} else if (Physics->dt>Numerics->dtmax) {
		//	Physics.dt = dtmax;
	}

}





