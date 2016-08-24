/*
 * Physics->c
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
	Physics->rho0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );



#if (HEAT)
	Physics->k 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->T0 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->DT 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
#endif

	Physics->Plitho 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

#if (DARCY)

	Physics->Pc 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->divV0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->Pc0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->DPc 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->Pf 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->phi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->phi0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->Dphi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->perm0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // permeability
	Physics->perm 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // permeability
	Physics->eta_b 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // bulk viscosity
	Physics->B				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // elastic bulk modulus

#endif





	//Physics->Pc0 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	//Physics->DPc 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->G 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->cohesion 		= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->frictionAngle 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );


	Physics->phase 			= (int*) 	malloc( Grid->nECTot * sizeof(int) );


	// Initialize stuff
	//int i;
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0;
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0;
	}

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nECTot; ++i) {
		//Physics->P[i]  = 0;

#if (HEAT)
		Physics->T[i]  = 0;
		Physics->DT[i] = 0;
#endif

		Physics->P[i] = 0;

		//Physics->eta[i] = 0;
		//Physics->rho[i] = 0;
#if (DARCY)
		Physics->divV0[i] = 0;

		Physics->Pf  [i] = 0;
		Physics->Pc  [i] = 0;
		Physics->Pc0 [i] = 0;
		Physics->DPc [i] = 0;
		Physics->phi [i] = 0;
		Physics->phi0[i] = 0;
#endif

		Physics->sigma_xx_0[i] = 0;
		Physics->Dsigma_xx_0[i] = 0;
	}

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->sigma_xy_0[i] = 0;
		Physics->Dsigma_xy_0[i] = 0;
	}



	Physics->dtMaxwellMin = 1E-100;
	Physics->dtMaxwellMax = 1E100;





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

	free( Physics->rho0 );



#if (HEAT)
	free( Physics->k );
	free(Physics->T );
	free(Physics->DT );
#endif

	free(Physics->phase);

	free(Physics->G );


	free(Physics->sigma_xx_0 );
	free(Physics->sigma_xy_0 );
	free(Physics->Dsigma_xx_0 );
	free(Physics->Dsigma_xy_0 );

	free(Physics->cohesion);
	free(Physics->frictionAngle);

	free(Physics->Plitho);
	// Darcy
#if (DARCY)

	free(Physics->Pc);

	free(Physics->divV0);

	free(Physics->Pc0);
	free(Physics->DPc);
	free(Physics->Pf);
	free(Physics->phi);
	free(Physics->phi0);
	free(Physics->perm0);
	free(Physics->perm);
	free(Physics->eta_b);
	free(Physics->B);

#endif






}



void Physics_initPToLithostatic(Physics* Physics, Grid* Grid)
{
	int ix, iy, iCell;
	//printf("=== P ===\n");
	compute rho_g_h, rhof_g_h;
	// Set Temp to zero (the interpolation forced the ghost values to have a dirichlet value follow the dirichlet)

	/*
	// Initialize the pressure at the lithostatic pressure
	// =========================

	// Initialize P at the lithostatic pressure
	// Contribution of y
	if (Physics->g[0]>0) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {

			iy = 0;
			iCell = ix + iy*Grid->nxEC;


			rho_g_h  = -0.5* (Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->DYEC[0]);
#if (DARCY)
			rhof_g_h = -0.5* (Physics->rho_f * fabs(Physics->g[1]) * Grid->DYEC[0]);
#endif

			Physics->P[iCell] = 1*rho_g_h;
			for (iy = 1; iy < Grid->nyEC; ++iy) {
				iCell = ix + iy*Grid->nxEC;
				rho_g_h += Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->DYEC[iy-1];//Grid->dy;

				Physics->P[iCell] = 1*rho_g_h;
#if (DARCY)
				rhof_g_h += Physics->rho_f* fabs(Physics->g[1]) * Grid->DYEC[iy-1];//Grid->dy;
				Physics->Pf[iCell] = 1*rhof_g_h;
#endif
			}
			//printf("\n");
		}
	} else {
		for (ix = 0; ix < Grid->nxEC; ++ix) {

			iy = Grid->nyEC-1;
			iCell = ix + iy*Grid->nxEC;
			rho_g_h  = -0.5* (Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->DYEC[Grid->nyEC-2]);
#if (DARCY)
			rhof_g_h = -0.5* (Physics->rho_f * fabs(Physics->g[1]) * Grid->DYEC[Grid->nyEC-2]);
#endif
			Physics->P[iCell] = 1*rho_g_h;
			for (iy = Grid->nyEC-2; iy >= 0; --iy) {
				iCell = ix + iy*Grid->nxEC;
				rho_g_h += Physics->rho[iCell] * fabs(Physics->g[1]) * Grid->DYEC[iy];


				Physics->P[iCell] = 1*rho_g_h;
#if (DARCY)
				rhof_g_h += Physics->rho_f* fabs(Physics->g[1]) * Grid->DYEC[iy];//Grid->dy;
				Physics->Pf[iCell] = 1*rhof_g_h;
#endif

			}
			//printf("\n");
		}
	}


	// Contribution of x // in case the gravity field is not vertical
	// be careful adding contribution from left to right. This assumes the model is diping right.
	// If the model dips in the other direction the loop should be from right to left
	if (Physics->g[0]>=0) {
		for (iy = 0; iy < Grid->nyEC; ++iy) {

			ix = 0;
			iCell = ix + iy*Grid->nxEC;
			rho_g_h  = -0.5* (Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->DXEC[0]);
#if (DARCY)
			rhof_g_h = -0.5* (Physics->rho_f * fabs(Physics->g[0]) * Grid->DXEC[0]);
#endif
			//Physics->P[iCell] = 1*rho_g_h;
			for (ix = 1; ix < Grid->nxEC; ++ix) {
				iCell = ix + iy*Grid->nxEC;
				rho_g_h += Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->DXEC[ix-1];

				//printf("%.2e  ", Physics->P[iCell]);
				Physics->P[iCell] += 1*rho_g_h;
#if (DARCY)
				rhof_g_h += Physics->rho_f * fabs(Physics->g[0]) * Grid->DXEC[ix-1];
				Physics->Pf[iCell] += 1*rhof_g_h;
#endif
			}
			//printf("\n");
		}
	} else {
		for (iy = 0; iy < Grid->nyEC; ++iy) {

			ix = Grid->nxEC-1;
			iCell = ix + iy*Grid->nxEC;
			rho_g_h  = -0.5* (Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->DXEC[Grid->nxEC-2]);
#if (DARCY)
			rhof_g_h = -0.5* (Physics->rho_f * fabs(Physics->g[0]) * Grid->DXEC[Grid->nxEC-2]);
#endif
			//Physics->P[iCell] = 1*rho_g_h;
			for (ix = Grid->nxEC-2; ix >=0; --ix) {
				iCell = ix + iy*Grid->nxEC;
				rho_g_h += Physics->rho[iCell] * fabs(Physics->g[0]) * Grid->DXEC[ix];

				//printf("%.2e  ", Physics->P[iCell]);
				Physics->P[iCell] += 1*rho_g_h;
#if (DARCY)
				rhof_g_h += Physics->rho_f * fabs(Physics->g[0]) * Grid->DXEC[ix];
				Physics->Pf[iCell] += 1*rhof_g_h;
#endif
			}
			//printf("\n");
		}
	}


#if (DARCY)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->Pc[iCell] = Physics->P[iCell] - Physics->Pf[iCell];
		Physics->Pc0[iCell] = Physics->P[iCell] - Physics->Pf[iCell];
		Physics->DPc[iCell] = Physics->P[iCell] - Physics->Pf[iCell];
	}


#endif




	if (DEBUG) {
		printf("=== init P to Litho P here ===\n");
		int C;
		// Check P
		// =========================
		printf("=== P here ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->P[C]);
				C++;
			}
			printf("\n");
		}
#if (DARCY)
		// Check Pf
		// =========================
		printf("=== Pf ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Pf[C]);
				C++;
			}
			printf("\n");
		}

		// Check Pc
		// =========================
		printf("=== Pc ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Pc[C]);
				C++;
			}
			printf("\n");
		}

#endif
	}


*/

		Physics_computePlitho(Physics, Grid);
		for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
			Physics->P[iCell] = Physics->Plitho[iCell];
#if (DARCY)
			Physics->Pf[iCell] = Physics->Plitho[iCell];
			Physics->Pc[iCell] = 0.0;
			Physics->Pc0[iCell] = 0.0;
			Physics->DPc[iCell] = 0.0;
#endif
		}

		if (DEBUG) {
			// Check P
			// =========================
			printf("=== P here ===\n");
			int C = 0;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					printf("%.3e  ", Physics->P[C]);
					C++;
				}
				printf("\n");
			}
#if (DARCY)
			// Check Pf
			// =========================
			printf("=== Pf ===\n");
			C = 0;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					printf("%.3e  ", Physics->Pf[C]);
					C++;
				}
				printf("\n");
			}

			// Check Pc
			// =========================
			printf("=== Pc ===\n");
			C = 0;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					printf("%.3e  ", Physics->Pc[C]);
					C++;
				}
				printf("\n");
			}

#endif
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


	compute* sumOfWeights 	= (compute*) calloc(nNeighbours * Grid->nECTot , sizeof(compute));


	compute* eta0 			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* n    			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* rho0  			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* G    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* cohesion 		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* frictionAngle 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* sigma_xx_0   	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));


#if (HEAT)
	compute* T    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* k    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#endif

#if (DARCY)
	compute* Pc0   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* phi0   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* perm0  		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* eta_b  		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* B  			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#endif


	compute* sigma_xy_0   	= (compute*) malloc(nNeighbours * Grid->nSTot * sizeof(compute));


	// Reinitialize Physics array

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		eta0[i] = 0;
		n[i] = 0;
		rho0[i] = 0;


		G  [i] = 0;
		sigma_xx_0 [i] = 0;
		cohesion[i] = 0;
		frictionAngle[i] = 0;


		sumOfWeights[i] = 0;
#if (HEAT)
		k  [i] = 0;
		T  [i] = 0;
#endif
#if (DARCY)
		Pc0[i] 		= 0;
		phi0[i] 		= 0;
		perm0[i] 	= 0;
		eta_b[i] 	= 0;
		B[i] 		= 0;
#endif
	}


#pragma omp parallel for private(i) schedule(static,32)
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





	printf("Main loop\n");
	int iArr;
	//printf("=== Part Temp ===\n");
	// Loop through inner nodes
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 0; ix < Grid->nxS; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			//printf("ix = %i, iy = %i\n", ix, iy);
			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix  ]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy  ]);
				}


				phase = thisParticle->phase;
				//printf("phase = %i, locX= %.2f, locY=%.2f  \n", thisParticle->phase,locX,locY);
				// Loop through neighbours
				//printf("B\n");
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);

					weight = fabs((locX + xMod[i]*1.0)   *   (locY + yMod[i]*1.0));

					eta0			[iCell*4+i] += MatProps->eta0[phase] * weight;
					n				[iCell*4+i] += MatProps->n   [phase] * weight;


					G				[iCell*4+i] += 1/MatProps->G [phase] * weight; // harmonic average
					cohesion		[iCell*4+i] += MatProps->cohesion[phase] * weight;
					frictionAngle	[iCell*4+i] += MatProps->frictionAngle[phase] * weight;


					sigma_xx_0 		[iCell*4+i] += thisParticle->sigma_xx_0 * weight;




					rho0			[iCell*4+i] += MatProps->rho0[phase]*weight;//* (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;

#if (HEAT)
					rho0			[iCell*4+i] += MatProps->rho0[phase] * weight * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell]);
					k				[iCell*4+i] += MatProps->k   [phase] * weight;
					T 				[iCell*4+i] += thisParticle->T * weight;
#else

#endif


#if (DARCY)
					Pc0				[iCell*4+i] += thisParticle->Pc0 * weight;
					phi0				[iCell*4+i] += thisParticle->phi * weight;

					perm0			[iCell*4+i] += MatProps->perm0   [phase] * weight;
					eta_b			[iCell*4+i] += MatProps->eta_b   [phase] * weight;
					B				[iCell*4+i] += MatProps->B   	 [phase] * weight;

#endif

					sumOfWeights	[iCell*4+i] += weight;

				}
				thisParticle = thisParticle->next;
			}
		}
	}



	//printf("Left Right contrib\n");
	// Add contribution from the other side in the case of periodic BC

	int nPointers = 7;
	int nPointersArithm = 6;
#if (HEAT)
	nPointers += 2;
	nPointersArithm += 2;
#endif
#if (DARCY)
	nPointers += 5;
	nPointersArithm += 5;
#endif
	compute** ArrayOfPointers;
	ArrayOfPointers = malloc((nPointers+1) * sizeof(compute*)); // the last one is sumofWeights
	compute** ArrayOfPointersPhysics;
	ArrayOfPointersPhysics = malloc((nPointers) * sizeof(compute*));

	// In this section goes the arrays for which arithmetic averaging is used
	// =======================
	ArrayOfPointers			[ 0] = eta0;
	ArrayOfPointersPhysics	[ 0] = Physics->eta0;
	ArrayOfPointers			[ 1] = n;
	ArrayOfPointersPhysics	[ 1] = Physics->n;

	ArrayOfPointers			[ 2] = cohesion;
	ArrayOfPointersPhysics	[ 2] = Physics->cohesion;
	ArrayOfPointers			[ 3] = frictionAngle;
	ArrayOfPointersPhysics	[ 3] = Physics->frictionAngle;

	ArrayOfPointers			[ 4] = rho0;
	ArrayOfPointersPhysics	[ 4] = Physics->rho0;
	ArrayOfPointers			[ 5] = sigma_xx_0;
	ArrayOfPointersPhysics	[ 5] = Physics->sigma_xx_0;



	i = 5;
#if (HEAT)
	i++;
	ArrayOfPointers			[ i] = k;
	ArrayOfPointersPhysics	[ i] = Physics->k;
	i++;
	ArrayOfPointers			[ i] = T;
	ArrayOfPointersPhysics	[ i] = Physics->T;
#endif



#if (DARCY)
	i++;
	ArrayOfPointers			[ i] = Pc0;
	ArrayOfPointersPhysics	[ i] = Physics->Pc0;
	i++;
	ArrayOfPointers			[ i] = phi0;
	ArrayOfPointersPhysics	[ i] = Physics->phi0;
	i++;
	ArrayOfPointers			[ i] = perm0;
	ArrayOfPointersPhysics	[ i] = Physics->perm0;
	i++;
	ArrayOfPointers			[ i] = eta_b;
	ArrayOfPointersPhysics	[ i] = Physics->eta_b;
	i++;
	ArrayOfPointers			[ i] = B;
	ArrayOfPointersPhysics	[ i] = Physics->B;

#endif
	// =======================

	// In this section goes the arrays for which harmonic averaging is used
	// =======================
	i++;
	ArrayOfPointers			[ i] = G;
	ArrayOfPointersPhysics	[ i] = Physics->G;

	// =======================


	i++;
	ArrayOfPointers			[ i] = sumOfWeights;



	int iPtr;
	if(BCStokes->SetupType==SimpleShearPeriodic) {
		int iCellS, iCellD, j;
#pragma omp parallel for private(iy, j, iCellS, iCellD,i, iPtr) schedule(static,32)
		for (iy = 0; iy < Grid->nyEC; ++iy) {

			for (j = 0; j<2; ++j) {
				if (j==0) {
					iCellS = 0 + iy*Grid->nxEC;
					iCellD = Grid->nxEC-2 + iy*Grid->nxEC;
				} else {
					iCellD = 1 + iy*Grid->nxEC;
					iCellS = Grid->nxEC-1 + iy*Grid->nxEC;
				}


				for (iPtr = 0; iPtr < nPointers+1 ; ++iPtr) {
					for (i = 0; i < 4; ++i) {
						ArrayOfPointers[iPtr][iCellD*4+i] += ArrayOfPointers[iPtr][iCellS*4+i];
						ArrayOfPointers[iPtr][iCellS*4+i]  = ArrayOfPointers[iPtr][iCellD*4+i];
					}
				}
			}
		}
	}




	//printf("Left Right contrib end\n");


	compute sum;
#pragma omp parallel for private(iCell, sum, I, iPtr) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {

		I = 4*iCell;
		sum = sumOfWeights[I+0] + sumOfWeights[I+1] + sumOfWeights[I+2] + sumOfWeights[I+3];
		//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		// Arithmetic averaging
		for (iPtr = 0; iPtr < nPointersArithm; ++iPtr) {
			ArrayOfPointersPhysics[iPtr][iCell] = (ArrayOfPointers[iPtr][I+0] + ArrayOfPointers[iPtr][I+1] + ArrayOfPointers[iPtr][I+2] + ArrayOfPointers[iPtr][I+3]) /sum;
		}
		// Harmonic averaging
		for (iPtr = nPointersArithm; iPtr < nPointers; ++iPtr) {
			ArrayOfPointersPhysics[iPtr][iCell] = sum/ (ArrayOfPointers[iPtr][I+0] + ArrayOfPointers[iPtr][I+1] + ArrayOfPointers[iPtr][I+2] + ArrayOfPointers[iPtr][I+3]);
		}

	}



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
		for (iPtr = 0; iPtr < nPointers; ++iPtr) {
			ArrayOfPointersPhysics[iPtr][I] = ArrayOfPointersPhysics[iPtr][INeigh];
		}

#if (HEAT)
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->DYEC[iy];
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
		for (iPtr = 0; iPtr < nPointers; ++iPtr) {
			ArrayOfPointersPhysics[iPtr][I] = ArrayOfPointersPhysics[iPtr][INeigh];
		}
#if (HEAT)
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->DYEC[iy];
		}
#endif
	}


	if (BCStokes->SetupType!=SimpleShearPeriodic) {
		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;
			INeigh =   ix+1 + (iy)*Grid->nxEC  ;
			for (iPtr = 0; iPtr < nPointers; ++iPtr) {
				ArrayOfPointersPhysics[iPtr][I] = ArrayOfPointersPhysics[iPtr][INeigh];
			}
#if (HEAT)
			IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
			if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
				Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
			}
			else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->DXEC[ix];
			}
#endif
		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			for (iPtr = 0; iPtr < nPointers; ++iPtr) {
				ArrayOfPointersPhysics[iPtr][I] = ArrayOfPointersPhysics[iPtr][INeigh];
			}

#if (HEAT)
			if (BCThermal->SetupType!=SimpleShearPeriodic) {
				IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
				if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
					Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
				}
				else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
					Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->DXEC[ix-1];
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
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, weight, signX, signY, iNodeNeigh) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 0; ix < Grid->nxS; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================

			while (thisParticle!=NULL) {
				//locX = (thisParticle->x-Grid->xmin)/dx - ix;
				//locY = (thisParticle->y-Grid->ymin)/dy - iy;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					signX = -1;
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					signX = 1;
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					signY = -1;
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					signY = 1;
					locY = 2.0*(locY/Grid->DYS[iy]);
				}


				phase = thisParticle->phase;



				for (i=0; i<4; i++) {
					iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

					if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
						printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
						//printf("IX = %i, IY = %i, locX = %.2e, locY = %.2e, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
						printf("thisParticle->x = %.2e , y = %.2e \n", thisParticle->x, thisParticle->y);
						exit(0);
					}

					locX = fabs(locX);
					locY = fabs(locY);


					//printf("iNodeNeigh = %i, signX = %i, signY = %i\n", iNodeNeigh, signX, signY);
					weight = (locX + xMod[i]*1.0)   *   (locY + yMod[i]*1.0);

					sigma_xy_0 	[iNodeNeigh*4+i] += thisParticle->sigma_xy_0 * weight;
					sumOfWeights[iNodeNeigh*4+i] += weight; // using the same arrays



				}
				thisParticle = thisParticle->next;
			}
		}
	}

	if(BCStokes->SetupType==SimpleShearPeriodic) {
		int iCellS, iCellD;
#pragma omp parallel for private(iy, iCellS, iCellD,i) schedule(static,32)
		for (iy = 0; iy < Grid->nyS; ++iy) {

			iCellS = 0 + iy*Grid->nxS; // Source
			iCellD = Grid->nxS-1 + iy*Grid->nxS; // destination

			//printf("before: sigma_xy_0[iCellS*4] = %.2e\n",sigma_xy_0[iCellS*4]);
			for (iPtr = 0; iPtr < nPointers+1 ; ++iPtr) {
				for (i = 0; i < 4; ++i) {
					sigma_xy_0 	[iCellD*4+i] += sigma_xy_0  [iCellS*4+i];
					sigma_xy_0 	[iCellS*4+i]  = sigma_xy_0  [iCellD*4+i];
					sumOfWeights[iCellD*4+i] += sumOfWeights[iCellS*4+i];
					sumOfWeights[iCellS*4+i]  = sumOfWeights[iCellD*4+i];
				}
			}
			//printf("after: sigma_xy_0[iCellS*4] = %.2e\n",sigma_xy_0[iCellS*4]);

		}
	}






	printf("end first loop for sigma_xy\n");
#pragma omp parallel for private(iNode, I, sum) schedule(static,32)
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

	printf("end filling loop for sigma_xy\n");









	if (DEBUG) {

		printf("=== Check eta0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->eta0[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== Check rho0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->rho0[C]);
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
				printf("%.2e  ", Physics->k[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check T 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->T[C]);
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
				printf("%.2e  ", Physics->G[C]);
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

#if (DARCY)
		printf("=== Check perm0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->perm0[C]);
				C++;
			}
			printf("\n");
		}


#endif

	}



	free(sumOfWeights);
	free(eta0);
	free(n);
	free(rho0);


	free(G);
	free(sigma_xx_0);
	free(sigma_xy_0);
	free(cohesion);
	free(frictionAngle);



#if (HEAT)
	free(T);
	free(k);
#endif

#if (DARCY)
	free(Pc0);
	free(phi0);

	free(perm0);
	free(eta_b);
	free(B);

#endif


	free(ArrayOfPointers);
	free(ArrayOfPointersPhysics);



}


















#if (HEAT)
void Physics_interpTempFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	compute* DT_sub_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* sumOfWeights_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* DT_rem_OnTheCells = (compute*) malloc( Grid->nECTot *sizeof(compute) );

	int i, iCell;

#pragma omp parallel for private(i) schedule(static,32)
	for (i=0;i<4*Grid->nECTot;++i) {
		DT_sub_OnTheCells[i] = 0;
		sumOfWeights_OnTheCells[i] = 0;
	}

	compute TFromNodes, DT_sub_OnThisPart, PFromNodes;
	compute rhoParticle;
	compute dtDiff;
	compute d = 0.98;

	int phase;

	compute weight;
	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right

	int xModCell[4], yModCell[4];
	xModCell[0] = -1; yModCell[0] = -1;
	xModCell[1] =  1; yModCell[1] = -1;
	xModCell[2] = -1; yModCell[2] =  1;
	xModCell[3] =  1; yModCell[3] =  1;

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, phase, locX, locY, TFromNodes, PFromNodes, rhoParticle, dtDiff, DT_sub_OnThisPart, i, iCell, weight) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			phase = thisParticle->phase;
			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				//locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				//locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}


				TFromNodes  = 		( .25*(1.0-locX)*(1.0-locY)*Physics->T[ix  +(iy  )*Grid->nxEC]
																		   + .25*(1.0-locX)*(1.0+locY)*Physics->T[ix  +(iy+1)*Grid->nxEC]
																												  + .25*(1.0+locX)*(1.0+locY)*Physics->T[ix+1+(iy+1)*Grid->nxEC]
																																						 + .25*(1.0+locX)*(1.0-locY)*Physics->T[ix+1+(iy  )*Grid->nxEC] );
				PFromNodes  = 		( .25*(1.0-locX)*(1.0-locY)*Physics->P[ix  +(iy  )*Grid->nxEC]
																		   + .25*(1.0-locX)*(1.0+locY)*Physics->P[ix  +(iy+1)*Grid->nxEC]
																												  + .25*(1.0+locX)*(1.0+locY)*Physics->P[ix+1+(iy+1)*Grid->nxEC]
																																						 + .25*(1.0+locX)*(1.0-locY)*Physics->P[ix+1+(iy  )*Grid->nxEC] );


				rhoParticle = MatProps->rho0[phase] * (1+MatProps->beta[phase]*PFromNodes) * (1-MatProps->alpha[phase]*thisParticle->T);

				dtDiff = (Physics->Cp*rhoParticle)/(  MatProps->k[phase]*( 2/(Grid->dx*Grid->dx) + 2/(Grid->dy*Grid->dy) )  );

				DT_sub_OnThisPart = ( TFromNodes - thisParticle->T ) * ( 1 - exp(-d * Physics->dtT/dtDiff) );


				// redefine locX, locY (used to compute surface based weight, not used as weight directly)
				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				// Interp Dsigma_xx_sub from particles to Cells
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
					weight = fabs((locX + xModCell[i]*0.5)   *   (locY + yModCell[i]*0.5));

					DT_sub_OnTheCells[iCell*4+i] += DT_sub_OnThisPart * weight;
					sumOfWeights_OnTheCells	[iCell*4+i] += weight;

				}

				thisParticle->T += DT_sub_OnThisPart;

				thisParticle = thisParticle->next;
			}
		}
	}


	compute DT_sub_OnThisCell;
	compute sum;
	int I;
#pragma omp parallel for private(iCell, I, sum, DT_sub_OnThisCell) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		I = 4*iCell;
		sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];
		//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		DT_sub_OnThisCell = ( DT_sub_OnTheCells[I+0] + DT_sub_OnTheCells[I+1] + DT_sub_OnTheCells[I+2] + DT_sub_OnTheCells[I+3]) / sum ; // harmonic average
		DT_rem_OnTheCells[iCell] = Physics->DT[iCell] - DT_sub_OnThisCell;

		//printf("Physics->simga_xx_0[%i] %.2e, sigma_xx_0 = %.2f %.2f %.2f %.2f, sum = %.2f\n", iCell, Physics->sigma_xx_0[iCell], sigma_xx_0[I+0], sigma_xx_0[I+1], sigma_xx_0[I+2], sigma_xx_0[I+3], sum);

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

		DT_rem_OnTheCells[I] = DT_rem_OnTheCells[INeigh];

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

		DT_rem_OnTheCells[I] = DT_rem_OnTheCells[INeigh];

	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;

		DT_rem_OnTheCells[I] = DT_rem_OnTheCells[INeigh];

	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;

		DT_rem_OnTheCells[I] = DT_rem_OnTheCells[INeigh];

	}


	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				//locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				//locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}

				//compute locX0 = locX;
				//compute locY0 = locY;


				thisParticle->T  += ( .25*(1.0-locX)*(1.0-locY)*DT_rem_OnTheCells[ix  +(iy  )*Grid->nxEC]
																				  + .25*(1.0-locX)*(1.0+locY)*DT_rem_OnTheCells[ix  +(iy+1)*Grid->nxEC]
																																+ .25*(1.0+locX)*(1.0+locY)*DT_rem_OnTheCells[ix+1+(iy+1)*Grid->nxEC]
																																											  + .25*(1.0+locX)*(1.0-locY)*DT_rem_OnTheCells[ix+1+(iy  )*Grid->nxEC] );



				thisParticle = thisParticle->next;
			}
		}
	}



	free(DT_sub_OnTheCells);
	free(sumOfWeights_OnTheCells);
	free(DT_rem_OnTheCells);

}
#endif


void Physics_interpPhiFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics)
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

				//locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				//locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}

#if (DARCY)


				thisParticle->Pc0 += ( .25*(1.0-locX)*(1.0-locY)*Physics->DPc[ix  +(iy  )*Grid->nxEC]
																			   + .25*(1.0-locX)*(1.0+locY)*Physics->DPc[ix  +(iy+1)*Grid->nxEC]
																														 + .25*(1.0+locX)*(1.0+locY)*Physics->DPc[ix+1+(iy+1)*Grid->nxEC]
																																								   + .25*(1.0+locX)*(1.0-locY)*Physics->DPc[ix+1+(iy  )*Grid->nxEC] );
				thisParticle->phi += ( .25*(1.0-locX)*(1.0-locY)*Physics->Dphi[ix  +(iy  )*Grid->nxEC]
																			   + .25*(1.0-locX)*(1.0+locY)*Physics->Dphi[ix  +(iy+1)*Grid->nxEC]
																														 + .25*(1.0+locX)*(1.0+locY)*Physics->Dphi[ix+1+(iy+1)*Grid->nxEC]
																																								   + .25*(1.0+locX)*(1.0-locY)*Physics->Dphi[ix+1+(iy  )*Grid->nxEC] );

				#endif

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

	compute d_ve = 0.98;
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
#pragma omp parallel for private(i) schedule(static,32)
	for (i=0;i<4*Grid->nSTot;++i) {
		Dsigma_xy_sub_OnTheNodes[i] = 0;
		sumOfWeights_OnTheNodes[i] = 0;
	}
#pragma omp parallel for private(i) schedule(static,32)
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

				//locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				//locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}



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
						//printf("IX = %i, IY = %i, locX = %.2e, locY = %.2e, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
						printf("thisParticle->x = %.2e , y = %.2e \n", thisParticle->x, thisParticle->y);
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


	int I;
	compute sum;
	compute Dsigma_xy_sub_OnThisNode;
#pragma omp parallel for private(iNode, I, sum, Dsigma_xy_sub_OnThisNode) schedule(static,32)
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
#pragma omp parallel for private(iCell, I, sum, Dsigma_xx_sub_OnThisCell) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
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

				//locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				//locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;
				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}



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










void Physics_get_VxVy_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int InoDir, INeigh;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = -1;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = -1;
	}
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P[i] = -1;
	}


	// Set Vx
	// =========================
	int IBC;

#pragma omp parallel for private(iy, ix, I, InoDir, IBC, INeigh) schedule(static,32) // maxVx would conflict
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			I = ix + iy*Grid->nxVx;
			InoDir = Numbering->map[I];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vx[I] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1; // BC nodes are numbered -1 to -n
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vx[I] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (iy==0)  // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVx  ];
					if (iy==Grid->nyVx-1)  // lower boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVx  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vx[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (iy==0)  // lower boundary
							Physics->Vx[I] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dy;
						if (iy==Grid->nyVx-1)  // lower boundary
							Physics->Vx[I] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dy;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}



		}
	}

	// Set Vy
	// =========================

	int IMap;
#pragma omp parallel for private(iy, ix, I, IMap, InoDir, IBC, INeigh) schedule(static,32) // maxVx would conflict
	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			IMap = ix + iy*Grid->nxVy + Grid->nVxTot;
			I = ix + iy*Grid->nxVy;

			InoDir = Numbering->map[IMap];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vy[I] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vy[I] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (ix==0)  // lower boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVy + Grid->nVxTot ];
					if (ix==Grid->nxVy-1)  // lower boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVy + Grid->nVxTot  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vy[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (ix==0)  // lower boundary
							Physics->Vy[I] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dx;
						if (ix==Grid->nxVy-1)  // lower boundary
							Physics->Vy[I] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dx;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}


		}
	}


	compute maxV;
	compute Vx, Vy;
	maxV = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			Vx = (Physics->Vx[ix-1+  iy   *Grid->nxVx]+Physics->Vx[ix+  iy   *Grid->nxVx])/2.0;
			Vy = (Physics->Vy[ix  + (iy-1)*Grid->nxVy]+Physics->Vx[ix+ (iy-1)*Grid->nxVy])/2.0;
			if (Vx*Vx+Vy*Vy>maxV) {
				maxV = Vx*Vx + Vy*Vy;
			}
		}
	}

	Physics->maxV = sqrt(maxV);






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

	}






}






void Physics_get_P_FromSolution(Physics* Physics, Grid* Grid, BC* BCStokes, Numbering* NumStokes, EqSystem* EqStokes, Numerics* Numerics)
{
	int iy, ix, I, InoDir, IBC, iCell;

	compute * thisP;
	int eq0;








#if (!DARCY)
	// /!\ For visuit's better if all sides are Neumann
	Physics_get_ECVal_FromSolution (Physics->P, 2, Grid, BCStokes, NumStokes, EqStokes);

	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	compute RefPressure = Physics->P[1 + (Grid->nyEC-1)*Grid->nxEC];
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure;
	}

#else

	int i;


	// save the value from the previous time step
	/*
	if (Numerics->itNonLin == -1) {
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->Pc0[i] = Physics->Pc[i];
		}
	}
	*/
	printf("Pf\n");
	Physics_get_ECVal_FromSolution (Physics->Pf, 2, Grid, BCStokes, NumStokes, EqStokes);
	printf("Pc\n");
	Physics_get_ECVal_FromSolution (Physics->Pc, 3, Grid, BCStokes, NumStokes, EqStokes);




	// Fill P, the total pressure
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = Physics->Pc[iCell] + Physics->Pf[iCell];
	}

	/*
	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	compute RefPressure = Physics->P[1 + (Grid->nyEC-1)*Grid->nxEC];
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure;
		//Physics->Pc [iCell] = Physics->Pc [iCell] - RefPressure;
		Physics->Pf [iCell] = Physics->Pf [iCell] - RefPressure;
	}
	*/


	// get the increment from the previous time step DT
	//if (Numerics->itNonLin == -1) {
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DPc[i] = Physics->Pc[i] - Physics->Pc0[i];
		}
	//}






#endif




	if (DEBUG) {
		int C;
		// Check P
		// =========================
		printf("=== P here ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->P[C]);
				C++;
			}
			printf("\n");
		}
#if (DARCY)
		// Check Pf
		// =========================
		printf("=== Pf ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Pf[C]);
				C++;
			}
			printf("\n");
		}

		// Check Pc
		// =========================
		printf("=== Pc ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Pc[C]);
				C++;
			}
			printf("\n");
		}
		// Check Pc0
		// =========================
		printf("=== Pc0 ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Pc0[C]);
				C++;
			}
			printf("\n");
		}
		// Check DPc
		// =========================
		printf("=== DPc ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->DPc[C]);
				C++;
			}
			printf("\n");
		}

#endif
	}









}


#if (HEAT)
void Physics_get_T_FromSolution(Physics* Physics, Grid* Grid, BC* BCThermal, Numbering* NumThermal, EqSystem* EqThermal, Numerics* Numerics)
{
	// Declarations
	// =========================
	int i;
	//int InoDir, INeigh;
	//compute maxVx = 0;
	//compute maxVy = 0;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================


	// save the value from the previous time step
	if (Numerics->itNonLin == -1) {
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->T0[i] = Physics->T[i];
		}
	}
	Physics_get_ECVal_FromSolution (Physics->T, 0, Grid, BCThermal, NumThermal, EqThermal);

	// get the increment from the previous time step DT
	if (Numerics->itNonLin == -1) {
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DT[i] = Physics->T[i] - Physics->T0[i];
		}
	}


	if (DEBUG) {
		// Check T
		// =========================
		printf("=== T ===\n");

		int C = 0;
		int iy, ix;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->T[C]);
				C++;
			}
			printf("\n");
		}
	}



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
#pragma omp parallel for private(iy, ix, iCell, Eps_xx, Z) schedule(static,32)
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





#pragma omp parallel for private(iy, ix, iNode, dVxdy, dVydx, Eps_xy, GShear, etaShear, Z) schedule(static,32)
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



void Physics_computeEta(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes,MatProps* MatProps)
{
	int iCell, iy, ix;
	compute sigma_y, sigmaII;
	//compute* EIIGrid = (compute*) malloc(Grid->nECTot*sizeof(compute));
	//Physics_computeStrainRateInvariant(Physics, Grid, EIIGrid);

	int C = 0;

	compute EII_visc, eta_visc, EII;
	compute sigma_xx, sigma_xy;

	compute alpha, sigma_xxT;






	//printf("timeStep = %i, itNonLin = %i\n", Numerics->timeStep, Numerics->itNonLin);

	Physics->dtMaxwellMin = 1E100;
	Physics->dtMaxwellMax = 0;

	compute dtMaxwell;
	compute corr, etaViscNew;
	compute tolerance = 1e-8;
	compute etaVisc0;
#pragma omp parallel for private(ix,iy, iCell, sigma_xy, sigma_xx, sigmaII, etaVisc0, corr, etaViscNew, sigma_y, EII_visc, EII) schedule(static,32)
	//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;

			//printf("Physics->eta_b[iCell] = %.2e, phi = %.2e\n", Physics->eta_b[iCell], Physics->phi[iCell]);

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

			if (Numerics->timeStep<=0 && Numerics->itNonLin == -1){
				Physics->etaVisc[iCell] = Physics->eta0[iCell];
				Physics->eta[iCell] = Physics->etaVisc[iCell];
#if (DARCY)
				Physics->eta_b[iCell] 	=  	Physics->eta0[iCell]/Physics->phi[iCell];
#endif

			} else {
				sigmaII = sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);



#if (DARCY)
				if (MatProps->isWater[Physics->phase[iCell]] || MatProps->isAir[Physics->phase[iCell]]){
					Physics->eta_b[iCell] 	=  	MatProps->eta_b[Physics->phase[iCell]];
				} else {
					Physics->eta_b[iCell] 	=  	Physics->eta0[iCell]/Physics->phi[iCell];
				}

				Physics->eta [iCell] 	= 	Physics->eta0[iCell] * exp(27.0*Physics->phi[iCell]);
#else


				etaVisc0 = Physics->etaVisc[iCell];
				corr = 2*etaVisc0; // dummy initial value, just needs to be higher than etaVisc0
				//C = 0;
				while (fabs(corr/etaVisc0)>tolerance) {

					EII_visc = sigmaII/(2*Physics->etaVisc[iCell]);
					etaViscNew = Physics->eta0[iCell] * pow(EII_visc/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);
					corr = etaViscNew-Physics->etaVisc[iCell];

					Physics->etaVisc[iCell] += 1.0*corr;

					//C++;
				}
				Physics->etaVisc[iCell] = (Physics->etaVisc[iCell] + etaVisc0)/2;
				/*
				if (ix==10 && iy==10) {
					printf("C = %i\n",C);
				}
				 */














				// Compute powerlaw rheology
				//Physics->etaVisc[iCell] = Physics->eta0[iCell] * pow(EII_visc/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);
				Physics->eta[iCell] = Physics->etaVisc[iCell];
#endif


				// Plasticity
				sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Physics->P[iCell] * sin(Physics->frictionAngle[iCell]);
				if (sigmaII>sigma_y) {
					Physics_computeStrainInvariantForOneCell(Physics, Grid, ix,iy, &EII);
					Physics->eta[iCell] = sigma_y / (2*EII);
					//sigmaII = sigma_y;
				}


			}

			dtMaxwell = Physics->eta[iCell]/Physics->G[iCell];
			if (dtMaxwell<Physics->dtMaxwellMin) {
				Physics->dtMaxwellMin = dtMaxwell;
			}
			if (dtMaxwell>Physics->dtMaxwellMax) {
				Physics->dtMaxwellMax = dtMaxwell;
			}

			if (Physics->eta[iCell]<Numerics->etaMin) {

				Physics->eta[iCell] = Numerics->etaMin;
			}
			else if (Physics->eta[iCell]>Numerics->etaMax) {
				Physics->eta[iCell] = Numerics->etaMax;
			}

#if (DARCY)

			if (Physics->eta_b[iCell]<Numerics->etaMin) {

				Physics->eta_b[iCell] = Numerics->etaMin;
			}
			else if (Physics->eta_b[iCell]>Numerics->etaMax) {
				Physics->eta_b[iCell] = Numerics->etaMax;
			}


#endif



		}
	}


	Physics_copyValuesToSides(Physics->eta, Grid, BCStokes);
#if (DARCY)
	Physics_copyValuesToSides(Physics->eta_b, Grid, BCStokes);
#endif













	if (DEBUG) {
		printf("===== Compute Eta =====\n");
		printf("=== Check eta0  ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->eta0[C]);
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
				printf("%.2e  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
		}

#if (DARCY)
		printf("=== Check eta_b ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->eta_b[C]);
				C++;
			}
			printf("\n");
		}
#endif


		printf("=== Check n ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->n[C]);
				C++;
			}
			printf("\n");
		}
	}


	//free(EIIGrid);

}



void Physics_changePhaseOfFaults(Physics* Physics, Grid* Grid, MatProps* MatProps, Particles* Particles)
{
	/*
	compute* EII = (compute*) malloc(Grid->nECTot*sizeof(compute));
	Physics_computeStrainRateInvariant(Physics, Grid, EII);

	SingleParticle* thisParticle = NULL;

	int ix, iy, iNode;
	compute EII_node;

#pragma omp parallel for private(iy, ix, iNode, thisParticle, EII_node) schedule(static,32)
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
	 */

}



void Physics_updateDt(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics)
{
	/*
	if (fabs(Physics->maxV)<1E-6)
		Physics->maxV = 1E-6;
	 */

	Physics->dtAdv 	= Numerics->CFL_fac*Numerics->dLmin/(Physics->maxV); // note: the min(dx,dy) is the char length, so = 1
	Physics->dtT 	= 10.0*Numerics->CFL_fac*fmin(Grid->dx, Grid->dy)/(3*min(MatProps->k,MatProps->nPhase));

#if (DARCY)
	Physics->dtDarcy 	= 10.00*Numerics->CFL_fac*fmin(Grid->dx, Grid->dy)/(3*Physics->minPerm);
#endif


	printf("dtAdv = %.2e, dtT = %.2e, Numerics->dLmin = %.2e, Numerics->CFL = %.2e, (Physics->maxV) = %.2e\n", Physics->dtAdv, Physics->dtT, Numerics->dLmin, Numerics->CFL_fac, (Physics->maxV));


	//printf("maxV = %.3em, Physics->dt = %.3e, Physics->dt(SCALED)= %.3e yr, dtmin = %.2e, dtmax = %.2e, dtMax = %.2e\n",fabs(Physics->maxV), Physics->dt, Physics->dt*Char.time/3600/24/365, dtmin, dtmax, dtMax);

	/*
	// Doing this forbids to effectively deactivate the elasticity
	if (Physics->dtAdv>10*Physics->dtMaxwellMin && Physics->dtAdv<10*Physics->dtMaxwellMax) {
		// Physics dt is significantly larger than the minimum maxwell time and significantly lower than the maximum maxwell time
		// i.e. = OK
		Physics->dt = Physics->dtAdv;
	} else {

		Physics->dt = (Physics->dtMaxwellMin+Physics->dtMaxwellMax)/2;
	}
	 */

	Physics->dt = Physics->dtAdv;



	//Physics->dtAdv 	= fmin(Physics->dt,Physics->dtAdv);
	//Physics->dtT 	= fmin(Physics->dtT,Physics->dtAdv);

	Physics->dt  =  fmin(Physics->dtT,Physics->dtAdv);
#if (DARCY)
	Physics->dt  =  fmin(Physics->dt,Physics->dtDarcy);
#endif

	Physics->dtAdv = Physics->dt;
	Physics->dtT = Physics->dt;

#if (DARCY)
	Physics->dtDarcy = Physics->dt;
#endif






	if (Physics->dt<Numerics->dtMin) {
		Physics->dt = Numerics->dtMin;
	} else if (Physics->dt>Numerics->dtMax) {
		Physics->dt = Numerics->dtMax;
	}
	if (Physics->dtAdv<Numerics->dtMin) {
		Physics->dtAdv = Numerics->dtMin;
	} else if (Physics->dtAdv>Numerics->dtMax) {
		Physics->dtAdv = Numerics->dtMax;
	}
	if (Physics->dtT<Numerics->dtMin) {
		Physics->dtT = Numerics->dtMin;
	} else if (Physics->dtT>Numerics->dtMax) {
		Physics->dtT = Numerics->dtMax;
	}
#if (DARCY)
	if (Physics->dtDarcy<Numerics->dtMin) {
		Physics->dtDarcy = Numerics->dtMin;
	} else if (Physics->dtDarcy>Numerics->dtMax) {
		Physics->dtDarcy = Numerics->dtMax;
	}
#endif


	printf("Physics->dt = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, Physics->dtAdv = %.2e, Physics->dtT = %.2e, Physics->dtDarcy = %.2e\n", Physics->dt, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Physics->dtAdv, Physics->dtT, Physics->dtDarcy);



}


#if (DARCY)
void Physics_computePerm(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes)
{
	Physics->minPerm = 1E100;
	int iy, ix;
	int iCell;
	compute phi;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			phi = Physics->phi[iCell];
			/*
			if (phi>Numerics->phiMax) {
				phi = Numerics->phiMax;
			} else if (phi<Numerics->phiMin) {
				phi = Numerics->phiMin;
			}
			*/
			Physics->perm[iCell] = Physics->perm0[iCell] ;// *  phi*phi*phi  *  (1.0-phi)*(1.0-phi);

			if (Physics->perm[iCell]<Physics->minPerm) {
				Physics->minPerm = Physics->perm[iCell];
			}

		}
	}

	if (DEBUG) {
		printf("=== Check perm  ===\n");
		int C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->perm[C]);
				C++;
			}
			printf("\n");
		}
	}
}


void Physics_computePhi(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes)
{


	int iy, ix;
	int iCell;
	compute dt = Physics->dt;
	int nxVx = Grid->nxVx;
	int nxVy = Grid->nxVy;
	compute dx, dy;
	compute divV;
	compute sum = 0.0;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			dx = Grid->DXS[ix-1];
			dy = Grid->DYS[iy-1];
			divV  = (  Physics->Vx[ix+iy*nxVx] - Physics->Vx[ix-1+ iy   *nxVx]  )/dx;
			divV += (  Physics->Vy[ix+iy*nxVy] - Physics->Vy[ix  +(iy-1)*nxVy]  )/dy;
			//                    Physics->phi[iCell] = 1 - exp(-divV*dt)*(1-Physics->phi0[iCell]);
			// Dphi/Dt = (1-phi)*divV




			Physics->phi[iCell] = Physics->phi0[iCell];// + dt*0.5*(          (1.0-Physics->phi0[iCell])*Physics->divV0[iCell] + (1.0-Physics->phi[iCell])*divV         );


			if (Physics->phi[iCell] > Numerics->phiMax) {
				Physics->phi[iCell] = Numerics->phiMax;
			} else if (Physics->phi[iCell] < Numerics->phiMin) {
				Physics->phi[iCell] = Numerics->phiMin;
			}




			Physics->Dphi[iCell] = Physics->phi[iCell] - Physics->phi0[iCell];


			/*
			if (iCell == 150) {
			printf("    divV = %.2e, phi0 = %.2e, phi = %.2e, dt = %.2e\n", divV, Physics->phi0[iCell], Physics->phi[iCell], dt);
			}
			*/
			sum += Physics->phi[iCell];
		}
		//printf("divV = %.2e\n",divV);
	}
	printf("                    sum phi = %.2e\n", sum);


	Physics_copyValuesToSides(Physics->phi, Grid, BCStokes);
	Physics_copyValuesToSides(Physics->Dphi, Grid, BCStokes);

if (DEBUG) {
	printf("=== Check phi  ===\n");
	int C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->phi[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== Check Dphi  ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->Dphi[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== Check phi0  ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->phi0[C]);
			C++;
		}
		printf("\n");
	}
}

}







void Physics_initPhi(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics)
{

	Physics->PfGrad_Air_X = 0.0;
	Physics->PfGrad_Air_Y = 0*1E-2;

	Numerics->phiMin = 0.000;
	Numerics->phiMax = 1.0-Numerics->phiMin;


	printf("in InitPhi\n");
	int type = 0; // 0, porosity wave; 1, with ocean

	if (type==0) {
		compute xc = Grid->xmin + (Grid->xmax - Grid->xmin)/2.0;
		compute yc = Grid->ymin + (Grid->ymax - Grid->ymin)/3.0;
		compute phiBackground = 0.01;
		compute A = 0.0*phiBackground;
		compute x = Grid->xmin;
		compute y = Grid->ymin;
		compute w = (Grid->xmax - Grid->xmin)/8.0;
		compute XFac = 0.0;
		compute YFac = 1.0;
		int iCell;
		int iy, ix;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				iCell = ix+iy*Grid->nxEC;
				Physics->phi [iCell] = phiBackground + A*exp(   - XFac* (x-xc)*(x-xc)/(2*w*w) - YFac* (y-yc)*(y-yc)/(2*w*w)      );
				if (y==yc) {
					//printf("Physics->Dphi [iCell] = %.2e, x = %.2e, y = %.2e, xc, = %.2e, yc = %.2e, w = %.2e\n",Physics->Dphi [iCell], x, y, xc, yc, w);
				}
				Physics->Dphi  [iCell]  = Physics->phi[iCell];
				//Physics->phi0  [iCell] = Physics->phi[iCell];
				if (ix<Grid->nxEC-1) {
					x += Grid->DXEC[ix];
				} else {
					x = Grid->xmin;
				}

			}
			if (iy<Grid->nyEC-1) {
				y+=Grid->DYEC[iy];
			}
		}


	} else if (type == 1) {


		int iCell;
		int iy, ix;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				iCell = ix+iy*Grid->nxEC;

				if (MatProps->isWater[Physics->phase[iCell]]) {
					Physics->Dphi[iCell] = 1.0;
				} else {
					Physics->Dphi[iCell] = 0.0;
				}
				Physics->phi  [iCell] = Physics->Dphi[iCell];

			}
		}



	} else {
		printf("error in Physics_initPhi: unknwon type\n");
		exit(0);
	}
	printf("Out of InitPhi|n");



	if (DEBUG) {
	printf("\n=== Init phi  ===\n");
	printf("=== Check phi  ===\n");
	int C = 0;
	int iy, ix;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->phi[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== Check Dphi  ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->Dphi[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== Check phi0  ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.5e  ", Physics->phi0[C]);
			C++;
		}
		printf("\n");
	}
}

}





#endif


void Physics_copyValuesToSides(compute* ECValues, Grid* Grid, BC* BC)
{

	// Replace boundary values by their neighbours
	int INeigh, iy, ix, I;

	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (BC->SetupType==SimpleShearPeriodic) {
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

		ECValues[I] = ECValues[INeigh];

	}

	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (BC->SetupType==SimpleShearPeriodic) {
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
		ECValues[I] = ECValues[INeigh];

	}


	if (BC->SetupType==SimpleShearPeriodic) {
		int Iidentical; // index of the identical node
		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;
			Iidentical =   Grid->nxEC-2 + (iy)*Grid->nxEC  ; //
			ECValues[I] = ECValues[Iidentical];

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			Iidentical =   1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[Iidentical];


		}
	}
	else {
		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;
			INeigh =   ix+1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[INeigh];

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[INeigh];


		}
	}
	printf("end neighbour stuff");
}

void Physics_copyValuesToSidesi(int* ECValues, Grid* Grid, BC* BC)
{

	// Replace boundary values by their neighbours
	int INeigh, iy, ix, I;

	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (BC->SetupType==SimpleShearPeriodic) {
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

		ECValues[I] = ECValues[INeigh];

	}

	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (BC->SetupType==SimpleShearPeriodic) {
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
		ECValues[I] = ECValues[INeigh];

	}


	if (BC->SetupType==SimpleShearPeriodic) {
		int Iidentical; // index of the identical node
		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;
			Iidentical =   Grid->nxEC-2 + (iy)*Grid->nxEC  ; //
			ECValues[I] = ECValues[Iidentical];

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			Iidentical =   1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[Iidentical];


		}
	}
	else {
		// left boundary
		ix = 0;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;
			INeigh =   ix+1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[INeigh];

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			ECValues[I] = ECValues[INeigh];


		}
	}
	printf("end neighbour stuff");
}









void Physics_computeRho(Physics* Physics, Grid* Grid)
{

	int iCell;

	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->rho[iCell] = Physics->rho0[iCell];

#if (DARCY)
		//Physics->rho[iCell] = Physics->rho0[iCell];
		Physics->rho[iCell] = (1.0 - Physics->phi[iCell])*Physics->rho0[iCell] + Physics->phi[iCell]*Physics->rho_f;
		//Physics->rho[iCell] = Physics->rho0[iCell] * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell]);
#endif

	}

if (DEBUG) {
	printf("=== Check phi  ===\n");
	int C = 0;
	int iy, ix;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.8e  ", Physics->phi[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== Check rho  ===\n");
	C = 0;
	//int iy, ix;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.8e  ", Physics->rho[C]);
			C++;
		}
		printf("\n");
	}
}

}




void Physics_get_ECVal_FromSolution (compute* Val, int ISub, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
// Where Val is the value to extract from the solution, and DVal the increment since the last time step, IStep is the index of the subsystem of equations
	int I, IBC, INeigh, iy, ix;
	int INumMap0 = Numbering->subEqSystem0Dir[ISub];
	//printf("eq0 = %i, ISub = %i\n", INumMap0, ISub);
	int iCell;
#pragma omp parallel for private(iy, ix, I, iCell, IBC, INeigh) schedule(static,32)
	for (iy = 0; iy<Grid->nyEC; iy++) {
		for (ix = 0; ix<Grid->nxEC; ix++) {
			iCell = ix + iy*Grid->nxEC;
			I = Numbering->map[iCell + INumMap0];
			//printf("I = %i \n", I);
			if (I>=0) {
				Val[iCell]  = EqSystem->x[I];
			}
			else {

				IBC = abs(I)-1; // BC nodes are numbered -1 to -n

				// Get neighbours index
				if (iy==0) { // lower boundary
					if (BC->SetupType==SimpleShearPeriodic){
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxEC  + INumMap0 ];
					} else {
						if (ix==0) {
							INeigh = Numbering->map[  ix+1 + (iy+1)*Grid->nxEC + INumMap0 ];
						} else if (ix==Grid->nxEC-1) {
							INeigh = Numbering->map[  ix-1 + (iy+1)*Grid->nxEC + INumMap0 ];
						} else {
							INeigh = Numbering->map[  ix + (iy+1)*Grid->nxEC + INumMap0 ];
						}
					}
				} else if (iy==Grid->nyEC-1)  { //  upper boundary
					if (BC->SetupType==SimpleShearPeriodic){
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxEC + INumMap0 ];
					} else {
						if (ix==0) {
							INeigh = Numbering->map[  ix+1 + (iy-1)*Grid->nxEC  + INumMap0];
						} else if (ix==Grid->nxEC-1) {
							INeigh = Numbering->map[  ix-1 + (iy-1)*Grid->nxEC + INumMap0 ];
						} else {
							INeigh = Numbering->map[  ix + (iy-1)*Grid->nxEC + INumMap0 ];
						}
					}
				} else if (ix==0) { // left boundary
					INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxEC + INumMap0 ];
				} else if (ix==Grid->nxEC-1) { // right boundary
					INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxEC + INumMap0 ];
				}


				if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Val[iCell] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
				//	printf("IBC %i is Dir Ghost\n",IBC);
				}
				else if (BC->type[IBC]==NeumannGhost) { // Neumann
					if (ix==0)  {// left or bottom boundary
						Val[iCell] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->DXEC[0];
					} else if (ix==Grid->nxEC-1) {
						Val[iCell] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->DXEC[Grid->nxEC-2];
					}
					if (iy==0) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->DYEC[0];
					} else if (iy==Grid->nyEC-1) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->DYEC[Grid->nyEC-2];
					}
				}
				else {
					printf("error in Physics_get_ECVal_FromSolution: unknown boundary type\n");
					exit(0);
				}
				//printf("C=%i, IBC=%i, Type=%i, value=%.2e, valueNeigh=%.2e, FinalValue=%.2e\n",C, IBC,BC->type[IBC], BC->value[IBC], EqThermal->x[INeigh], Physics->T[C]);


				//Physics->T[C] = BC->value[abs(I)];
			}
		}
	}

}



void Physics_getPhase (Physics* Physics, Grid* Grid, Particles* Particles, MatProps* MatProps, BC* BCStokes)
{
	int ix, iy, iCell, iNode;
	coord depth, y;

	SingleParticle* thisParticle;
	compute locX, locY;
	int IxNode[] = {-1,  0, -1, 0};
	int IyNode[] = {-1, -1,  0, 0};
	int iPhase;
	compute contribPhase[NB_PHASE_MAX];
	compute maxContrib;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {

			iCell = ix+iy*Grid->nxEC;

			// Reinitialize contribs
			// ===================
			for (iPhase=0;iPhase<MatProps->nPhase;++iPhase) {
				contribPhase[iPhase] = 0;
			}


			// Count contribs
			// ===================
			for (iNode = 0; iNode < 4; ++iNode) {
				thisParticle = Particles->linkHead[ix+IxNode[iNode] + (iy+IyNode[iNode])*Grid->nxS];
				while (thisParticle != NULL) {
					++contribPhase[thisParticle->phase];
					thisParticle = thisParticle->next;
				}


			}

			// Find the most prominent phase
			// ===================
			maxContrib = 0;
			for (iPhase=0;iPhase<MatProps->nPhase;++iPhase) {
				if (contribPhase[iPhase] > maxContrib) {
					Physics->phase[iCell] = iPhase;
				}
			}


		}
	}


	Physics_copyValuesToSidesi(Physics->phase,Grid,BCStokes);

	if (DEBUG) {
		int C;
		printf("=== Check Phase ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%i  ", Physics->phase[C]);
				C++;
			}
			printf("\n");
		}

	}

}



void Physics_computePlitho(Physics* Physics, Grid* Grid)
{
	int iy, ix, iCell, iCellS, iCellN, iCellW, iCellE;
	compute rho_g_h;
	int ixStart, ixEnd, ixInc;
	int iyStart, iyEnd, iyInc;

	int C;



printf("enter Plitho\n");

	// Contribution of gy
	if (Physics->g[1]>0){
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				iCell = ix + iy*Grid->nxEC;
				iCellS = ix + (iy-1)*Grid->nxEC;
				if (iy==0) {
					rho_g_h = Physics->rho[iCell] * Physics->g[1] * (-0.5*Grid->DYEC[iy] );
				} else {
					rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellS]) * Physics->g[1] * Grid->DYEC[iy-1] ;
				}
				Physics->Plitho[iCell] = rho_g_h;
			}
		}

	} else {

		for (ix = 0; ix < Grid->nxEC; ++ix) {
			for (iy = Grid->nyEC-1; iy >= 0; --iy) {

				iCell = ix + iy*Grid->nxEC;
				iCellN = ix + (iy+1)*Grid->nxEC;
				iCellS = ix + (iy-1)*Grid->nxEC;
				if (iy==Grid->nyEC-1) {
					rho_g_h = Physics->rho[iCell] * -Physics->g[1] * (-0.5*Grid->DYEC[iy-1] );
				} else {
					rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellN]) * -Physics->g[1] * Grid->DYEC[iy] ;
				}
				//printf("ix = %i, iy = %i, rhogh = %.2e, Physics->rho[iCell] = %.2e\n", ix, iy, rho_g_h,Physics->rho[iCell]);
				Physics->Plitho[iCell] = rho_g_h;
			}
		}
	}


	if (abs(Physics->g[0])>1E-8) {
		// Contribution of gx
		if (Physics->g[0]>0){
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					iCellW = ix-1 + (iy)*Grid->nxEC;
					if (ix==0) {
						rho_g_h = Physics->rho[iCell] * Physics->g[0] * (-0.5*Grid->DXEC[ix] );
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellW]) * Physics->g[0] * Grid->DXEC[ix-1] ;
					}
					Physics->Plitho[iCell] += rho_g_h;
				}
			}
		} else {

			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = Grid->nxEC-1; ix >= 0; --ix) {
					iCell = ix + iy*Grid->nxEC;
					iCellE = ix+1 + (iy)*Grid->nxEC;
					iCellW = ix-1 + (iy)*Grid->nxEC;
					if (ix==Grid->nxEC-1) {
						rho_g_h = Physics->rho[iCell] * -Physics->g[0] * (-0.5*Grid->DXEC[ix-1] );
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellE]) * -Physics->g[0] * Grid->DXEC[ix] ;
					}
					Physics->Plitho[iCell] += rho_g_h;
				}
			}
		}
	}



	if (DEBUG) {
		printf("out Plitho\n");

		printf("=== compute P litho ===\n");
		//int C;
		// Check P
		// =========================
		printf("=== Plitho here ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->Plitho[C]);
				C++;
			}
			printf("\n");
		}
	}


}







