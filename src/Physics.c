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
	Physics->khi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->etaVisc		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->eta0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->n 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->rho_g 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->rho0_g 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );



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

	Physics->DeltaP0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->DDeltaP 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->Pf 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->phi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->phi0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->Dphi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // fluid phase fraction
	Physics->perm0_eta_f 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // permeability
	Physics->perm_eta_f 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // permeability
	Physics->eta_b 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // bulk viscosity
	Physics->khi_b 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // bulk viscosity
	//Physics->B				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // elastic bulk modulus

#endif





	//Physics->DeltaP0 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	//Physics->DDeltaP 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->G 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->khiShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->etaShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );


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

		Physics->khi[i] = 0;

#if (HEAT)
		Physics->T[i]  = 0;
		Physics->DT[i] = 0;
#endif

		Physics->P[i] = 0.0;

		//Physics->eta[i] = 0;
		//Physics->rho[i] = 0;
#if (DARCY)
		Physics->divV0[i] = 0;

		Physics->Pf  [i] = 0;
		Physics->Pc  [i] = 0;
		Physics->DeltaP0 [i] = 0;
		Physics->DDeltaP [i] = 0;
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



	Physics->dtMaxwellMin = 1E+100;
	Physics->dtMaxwellMax = 1E-100;





}


void Physics_freeMemory(Physics* Physics)
{
	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P );

	free( Physics->eta );

	free( Physics->eta0 );

	free(Physics->etaVisc);

	free(Physics->etaShear);
	free( Physics->khi );
	free( Physics->khiShear );

	free( Physics->n );
	free( Physics->rho_g );

	free( Physics->rho0_g );



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

	free(Physics->DeltaP0);
	free(Physics->DDeltaP);
	free(Physics->Pf);
	free(Physics->phi);
	free(Physics->Dphi);
	free(Physics->phi0);
	free(Physics->perm0_eta_f);
	free(Physics->perm_eta_f);
	free(Physics->eta_b);
	//free(Physics->B);
	free(Physics->khi_b);
#endif






}



void Physics_initPToLithostatic(Physics* Physics, Grid* Grid)
{
	int ix, iy, iCell;
	//printf("=== P ===\n");
	//compute rho_g_h, rhof_g_h;
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
		Physics->DeltaP0[iCell] = Physics->P[iCell] - Physics->Pf[iCell];
		Physics->DDeltaP[iCell] = Physics->P[iCell] - Physics->Pf[iCell];
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
			Physics->DeltaP0[iCell] = 0.0;
			Physics->DDeltaP[iCell] = 0.0;
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

	//coord dx = Grid->dx;
	//coord dy = Grid->dy;


	compute* sumOfWeights 	= (compute*) calloc(nNeighbours * Grid->nECTot , sizeof(compute));


	compute* eta0 			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* n    			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* rho0_g  		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* G    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* cohesion 		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* frictionAngle 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* sigma_xx_0   	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));


#if (HEAT)
	compute* T    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* k    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#endif

#if (DARCY)
	compute* DeltaP0   		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* phi0   		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* perm0_eta_f  		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	//compute* eta_b  		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	//compute* B  			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
#endif


	compute* sigma_xy_0   	= (compute*) malloc(nNeighbours * Grid->nSTot * sizeof(compute));


	// Reinitialize Physics array

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		eta0[i] = 0;
		n[i] = 0;
		rho0_g[i] = 0;


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
		DeltaP0[i] 		= 0;
		phi0[i] 		= 0;
		perm0_eta_f[i] 	= 0;
		//eta_b[i] 	= 0;
		//B[i] 		= 0;
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
	//int iArr;
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




					rho0_g			[iCell*4+i] += MatProps->rho0_g[phase]*weight;//* (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;

#if (HEAT)
					rho0_g			[iCell*4+i] += MatProps->rho0_g[phase] * weight * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell]);
					k				[iCell*4+i] += MatProps->k   [phase] * weight;
					T 				[iCell*4+i] += thisParticle->T * weight;
#else

#endif


#if (DARCY)
					DeltaP0				[iCell*4+i] += thisParticle->DeltaP0 * weight;
					phi0			[iCell*4+i] += thisParticle->phi * weight;

					perm0_eta_f			[iCell*4+i] += MatProps->perm0_eta_f   [phase] * weight;
					//eta_b			[iCell*4+i] += MatProps->eta_b   [phase] * weight;
					//B				[iCell*4+i] += MatProps->B   	 [phase] * weight;

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
	nPointers += 3;
	nPointersArithm += 3;
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

	ArrayOfPointers			[ 4] = rho0_g;
	ArrayOfPointersPhysics	[ 4] = Physics->rho0_g;
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
	ArrayOfPointers			[ i] = DeltaP0;
	ArrayOfPointersPhysics	[ i] = Physics->DeltaP0;
	i++;
	ArrayOfPointers			[ i] = phi0;
	ArrayOfPointersPhysics	[ i] = Physics->phi0;
	i++;
	ArrayOfPointers			[ i] = perm0_eta_f;
	ArrayOfPointersPhysics	[ i] = Physics->perm0_eta_f;
	/*
	i++;
	ArrayOfPointers			[ i] = eta_b;
	ArrayOfPointersPhysics	[ i] = Physics->eta_b;

	i++;
	ArrayOfPointers			[ i] = B;
	ArrayOfPointersPhysics	[ i] = Physics->B;
	*/
#endif
	// =======================


	nPointersArithm = i+1;

	// In this section goes the arrays for which harmonic averaging is used
	// =======================
	i++;
	ArrayOfPointers			[ i] = G;
	ArrayOfPointersPhysics	[ i] = Physics->G;

	// =======================

	nPointers = i+1;
	i++;
	ArrayOfPointers			[ i] = sumOfWeights;



	int iPtr;
	if(Grid->isPeriodic) {
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
		if (Grid->isPeriodic) {
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
		if (Grid->isPeriodic) {
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


	if (!Grid->isPeriodic) {
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
			if (!Grid->isPeriodic) {
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
	//int Counter;
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

	if(Grid->isPeriodic) {
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
				printf("%.2e  ", Physics->rho0_g[C]);
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
		printf("=== Check perm0_eta_f 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->perm0_eta_f[C]);
				C++;
			}
			printf("\n");
		}


#endif

	}



	free(sumOfWeights);
	free(eta0);
	free(n);
	free(rho0_g);


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
	free(DeltaP0);
	free(phi0);

	free(perm0_eta_f);
	//free(eta_b);
	//free(B);

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


				DT_sub_OnThisPart = ( TFromNodes - thisParticle->T ) * ( 1 - exp(-d * Physics->dt/dtDiff) );


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

	//compute dx = Grid->dx;
	//compute dy = Grid->dy;





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


				thisParticle->DeltaP0 += ( .25*(1.0-locX)*(1.0-locY)*Physics->DDeltaP[ix  +(iy  )*Grid->nxEC]
																			   + .25*(1.0-locX)*(1.0+locY)*Physics->DDeltaP[ix  +(iy+1)*Grid->nxEC]
																														 + .25*(1.0+locX)*(1.0+locY)*Physics->DDeltaP[ix+1+(iy+1)*Grid->nxEC]
																																								   + .25*(1.0+locX)*(1.0-locY)*Physics->DDeltaP[ix+1+(iy  )*Grid->nxEC] );
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

	compute khi, eta_vp;

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

				khi  				  = ( .25*(1.0-locX)*(1.0-locY)*Physics->khi[ix  +(iy  )*Grid->nxEC]
										+ .25*(1.0-locX)*(1.0+locY)*Physics->khi[ix  +(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0+locY)*Physics->khi[ix+1+(iy+1)*Grid->nxEC]
										+ .25*(1.0+locX)*(1.0-locY)*Physics->khi[ix+1+(iy  )*Grid->nxEC] );

				eta_vp = 1.0 / (1.0/eta + 1.0/khi);

				//printf("eta_vp = %.2e\n",eta_vp);
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

				dtMaxwell = eta_vp/G;

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
	compute scale;
#pragma omp parallel for private(iy, ix, I, InoDir, IBC, INeigh, scale) schedule(static,32) // maxVx would conflict
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			I = ix + iy*Grid->nxVx;
			InoDir = Numbering->map[I];



			if (InoDir>=0) { // Not a Dirichlet node
				scale = 1.0;//EqSystem->S[InoDir];
				Physics->Vx[I] = EqSystem->x[InoDir]*scale;
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

					scale = 1.0;//EqSystem->S[INeigh];

					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vx[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh]*scale;
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (iy==0)  // lower boundary
							Physics->Vx[I] = EqSystem->x[INeigh]*scale - BC->value[IBC]*Grid->dy;
						if (iy==Grid->nyVx-1)  // lower boundary
							Physics->Vx[I] = EqSystem->x[INeigh]*scale + BC->value[IBC]*Grid->dy;
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
#pragma omp parallel for private(iy, ix, I, IMap, InoDir, IBC, INeigh, scale) schedule(static,32) // maxVx would conflict
	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			IMap = ix + iy*Grid->nxVy + Grid->nVxTot;
			I = ix + iy*Grid->nxVy;

			InoDir = Numbering->map[IMap];



			if (InoDir>=0) { // Not a Dirichlet node
				scale = 1.0;//EqSystem->S[InoDir];
				Physics->Vy[I] = EqSystem->x[InoDir]*scale;
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

					scale = 1.0;//EqSystem->S[INeigh];

					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vy[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh]*scale;
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (ix==0)  // lower boundary
							Physics->Vy[I] = EqSystem->x[INeigh]*scale - BC->value[IBC]*Grid->dx;
						if (ix==Grid->nxVy-1)  // lower boundary
							Physics->Vy[I] = EqSystem->x[INeigh]*scale + BC->value[IBC]*Grid->dx;
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
	int iy, ix, iCell;
	//int iy, ix, I, InoDir, IBC, iCell;
	//compute * thisP;
	//int eq0;








#if (!DARCY)


	// /!\ For visuit's better if all sides are Neumann
	Physics_get_ECVal_FromSolution (Physics->P, 2, Grid, BCStokes, NumStokes, EqStokes);

	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	compute RefPressure = Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure;
	}

#else

	int i;


	// save the value from the previous time step
	/*
	if (Numerics->itNonLin == -1) {
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DeltaP0[i] = Physics->Pc[i];
		}
	}
	*/
	//printf("Pf\n");
	Physics_get_ECVal_FromSolution (Physics->Pf, 2, Grid, BCStokes, NumStokes, EqStokes);
	//printf("Pc\n");
	Physics_get_ECVal_FromSolution (Physics->Pc, 3, Grid, BCStokes, NumStokes, EqStokes);



	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)

	// Ref = average top row
	/*
	compute RefPressure = 0.0;//Physics->Pf[1 + (Grid->nyEC-2)*Grid->nxEC];
	for (ix = 0; ix < Grid->nxEC; ++ix) {
		iCell = ix + (Grid->nyEC-2)*Grid->nxEC;
		RefPressure += Physics->Pf[iCell];
	}
	RefPressure /= Grid->nxEC;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->Pf [iCell] 	= Physics->Pf [iCell] - RefPressure;
	}
	*/



	/*
	// Pressure at the surface of the solid is 0
	ix = 1;
		iy = Grid->nyEC-1;
		compute RefPressurePc;
		compute RefPressurePf;
		do  {
			iy--;
			iCell = ix + iy*Grid->nxEC;
		} while (Physics->phase[iCell]==Physics->phaseAir || Physics->phase[iCell]==Physics->phaseWater );
		//iy--;
		//iCell = ix + iy*Grid->nxEC;
		RefPressurePc = Physics->Pc[iCell];
		RefPressurePf = Physics->Pf[iCell];

			for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
				Physics->Pc [iCell] 	= Physics->Pc [iCell] - RefPressurePc;
				Physics->Pf [iCell] 	= Physics->Pf [iCell] - RefPressurePf;
			}

	*/




	/*
	// Pressure at the surface of the solid is 0
	for (ix = 0; ix < Grid->nxEC; ++ix) {
		iy = Grid->nyEC-1;
		compute RefPressurePc;
		compute RefPressurePf;
		do  {
			iy--;
			iCell = ix + iy*Grid->nxEC;
		} while (Physics->phase[iCell]==Physics->phaseAir || Physics->phase[iCell]==Physics->phaseWater );
		//iy--;
		//iCell = ix + iy*Grid->nxEC;
		RefPressurePc = Physics->Pc[iCell];
		RefPressurePf = Physics->Pf[iCell];
		int iySurface = iy;
		//for (iy = iyStart; iy >= 0; --iy) {
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			iCell = ix + iy*Grid->nxEC;
			if (iy <= iySurface) {


			//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
				Physics->Pc [iCell] 	= Physics->Pc [iCell] - RefPressurePc;
				Physics->Pf [iCell] 	= Physics->Pf [iCell] - RefPressurePf;
			} else {
				Physics->Pc [iCell] 	= 0.0;
				Physics->Pf [iCell] 	= 0.0;
			}
		}
	}
	*/






	/*
	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	RefPressure = Physics->Pf[1 + (Grid->nyEC-2)*Grid->nxEC];
	for (ix = 0; ix < Grid->nxEC; ++ix) {
		iCell = ix + (Grid->nyEC-2)*Grid->nxEC;
		RefPressure += Physics->Pf[iCell];
	}
	RefPressure /= Grid->nxEC;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->Pf [iCell] 	= Physics->Pf [iCell] - RefPressure;
	}

	 */




	// Fill P, the total pressure
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = Physics->Pc[iCell] + Physics->Pf[iCell];
	}



	/*
	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	int iCellTop;
	compute RefPressure = Physics->P[1 + (Grid->nyEC-1)*Grid->nxEC];
	for (ix = 0; i < Grid->nxEC; ++ix) {
		iCell 		= (Grid->nyEC-1)*Grid->nxEC + ix;
		//iCellTop 	= (Grid->nyEC-1)*Grid->nxEC + ix;
		RefPressure += (Physics->Pc[iCell] + Physics->Pc[iCellTop]);
	}
	RefPressure /= Grid->nxEC;


	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure;
		Physics->Pc [iCell] = Physics->Pc [iCell] - RefPressure;
		//Physics->Pf [iCell] = Physics->Pf [iCell] - RefPressure;
	}
	*/




	// get the increment from the previous time step DT
	//if (Numerics->itNonLin == -1) {
	compute phi;
		for (i = 0; i < Grid->nECTot; ++i) {
			phi = Physics->phi[i];
			Physics->DDeltaP[i] = Physics->Pc[i]/(1.0-phi) - Physics->DeltaP0[i];
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
				printf("%.3e  ", Physics->DeltaP0[C]);
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
				printf("%.3e  ", Physics->DDeltaP[C]);
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
	compute G, eta, khi;
	//compute phi;
	// compute stress
//#pragma omp parallel for private(iy, ix, iCell, Eps_xx, Z) schedule(static,32)
	compute dt = Physics->dt;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;
			Eps_xx 	= (Physics->Vx[ix + iy*Grid->nxVx] - Physics->Vx[ix-1 + iy*Grid->nxVx])/Grid->dx;//Grid->DXS[ix-1];


			Z = 1.0/(1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*dt));

			//Physics->Dsigma_xx_0[iCell] = ( 2.0*Physics->eta[iCell] * Eps_xx  -  Physics->sigma_xx_0[iCell] ) * Z;
			Physics->Dsigma_xx_0[iCell] = Z*(2.0*Eps_xx + Physics->sigma_xx_0[iCell]/(Physics->G[iCell]*dt)) - Physics->sigma_xx_0[iCell];

		}
	}
	//printf("Physics->sigma_xx_0[2+10*Grid->nxC] = %.2e, Physics->Dsigma_xx_0[2+10*Grid->nxC] = %.2e  ", Physics->sigma_xx_0[2+10*Grid->nxC], Physics->Dsigma_xx_0[2+10*Grid->nxC]);



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





//#pragma omp parallel for private(iy, ix, iNode, dVxdy, dVydx, Eps_xy, GShear, etaShear, Z) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]
					- Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]
					- Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Eps_xy = 0.5*(dVxdy+dVydx);




			G 	 	= shearValue(Physics->G, ix, iy, Grid->nxEC);
			eta 	= Physics->etaShear[iNode];
			khi 	= Physics->khiShear[iNode];

			Z = 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));

			Physics->Dsigma_xy_0[iNode] = Z * (2.0*Eps_xy + Physics->sigma_xy_0[iNode]/(G*dt)) - Physics->sigma_xy_0[iNode];
		}
	}


	if (DEBUG) {
		printf("=== Check Dsigma_xx ===\n");
		int C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->Dsigma_xx_0[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check Dsigma_xy ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.2e  ", Physics->Dsigma_xy_0[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check khi ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->khi[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check eta ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
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

void Physics_computeStrainRateInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* EII)
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


void Physics_computeStrainRateInvariantForOneNode(Physics* Physics, BC* BCStokes, Grid* Grid, int ix, int iy, compute* EII)
{
	compute dVxdy, dVydx, dVxdx, dVydy;




	if (Grid->isPeriodic) {
		if (ix == 0 || ix == Grid->nxS-1) {
			dVxdx = ( Physics->Vx[(1)+(iy+1)*Grid->nxVx] - Physics->Vx[(Grid->nxVx-1 -1)+(iy+1)*Grid->nxVx] +
					Physics->Vx[(1)+(iy  )*Grid->nxVx] - Physics->Vx[(Grid->nxVx-1 -1)+(iy  )*Grid->nxVx] )/4./Grid->dx;
		}
	}

	else {
		if (ix == 0) {
			dVxdx = ( Physics->Vx[(ix+1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] +
					  Physics->Vx[(ix+1)+(iy  )*Grid->nxVx] - Physics->Vx[(ix  )+(iy  )*Grid->nxVx] )/2./Grid->dx;
		} else if (ix == Grid->nxS-1) {
			dVxdx = ( Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] +
					  Physics->Vx[(ix  )+(iy  )*Grid->nxVx] - Physics->Vx[(ix-1)+(iy  )*Grid->nxVx] )/2./Grid->dx;
		} else {
			dVxdx = ( Physics->Vx[(ix+1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] +
					  Physics->Vx[(ix+1)+(iy  )*Grid->nxVx] - Physics->Vx[(ix-1)+(iy  )*Grid->nxVx] )/4./Grid->dx;
		}
	}

	if (iy == 0) {
		dVydy = ( Physics->Vy[(ix+1)+(iy+1)*Grid->nxVy] - Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] +
				  Physics->Vy[(ix  )+(iy+1)*Grid->nxVy] - Physics->Vy[(ix  )+(iy  )*Grid->nxVy] )/2./Grid->dy;
	} else if (iy == Grid->nyS-1) {
		dVydy = ( Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] +
				  Physics->Vy[(ix  )+(iy  )*Grid->nxVy] - Physics->Vy[(ix  )+(iy-1)*Grid->nxVy] )/2./Grid->dy;
	} else {
		dVydy = ( Physics->Vy[(ix+1)+(iy+1)*Grid->nxVy] - Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] +
			      Physics->Vy[(ix  )+(iy+1)*Grid->nxVy] - Physics->Vy[(ix  )+(iy-1)*Grid->nxVy] )/4./Grid->dy;
	}


	// the top and bottom row should never be needed

	dVxdy = (Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx]
		   - Physics->Vx[(ix  ) + (iy  )*Grid->nxVx])/Grid->dy;

	dVydx = (Physics->Vy[(ix+1) + (iy  )*Grid->nxVy]
		   - Physics->Vy[(ix  ) + (iy  )*Grid->nxVy])/Grid->dx;


	*EII = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  0.5*dVxdx*dVxdx  +  0.5*dVydy*dVydy);

}

void Physics_computeStressInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* SII) {


	compute dVxdy, dVydx, dVxdx, EII;
	compute sigma_xy0, sigma_xx0, sigmaII0;

	compute khi, eta, G, dt, phi, Z;
	compute Eff_strainRate;

	compute Eps_xx, Eps_xy;

	int iCell = ix + iy*Grid->nxEC;

	// Strain rates

	dVxdy = ( Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy-1)*Grid->nxVx] +
			Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy-1)*Grid->nxVx] )/4./Grid->dy;


	dVydx = ( Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] - Physics->Vy[(ix-1)+(iy-1)*Grid->nxVy] +
			Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix-1)+(iy  )*Grid->nxVy] )/4./Grid->dx;


	dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;


	Eps_xx = dVxdx;
	Eps_xy = 0.5*(dVxdy+dVydx);

	EII = sqrt(  Eps_xx*Eps_xx + Eps_xy*Eps_xy );




	// Old stress
	//
	sigma_xy0 = centerValue(Physics->sigma_xy_0,ix,iy,Grid->nxS);
	sigma_xx0 = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];

	sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + (sigma_xy0)*(sigma_xy0));


	khi 		= Physics->khi[iCell];
	eta 		= Physics->eta[iCell];
	G 		    = Physics->G[iCell];
	dt 			= Physics->dt;
	phi 		= 0.0;
#if (DARCY)
	phi = Physics->phi[iCell];
#endif


	Z 	= 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
	Eff_strainRate = sqrt(EII*EII + 1.0*Eps_xx*sigma_xx0/(G*dt) + 1.0*Eps_xy*sigma_xy0/(G*dt) + 1.0/4.0*(1.0/(G*dt))*(1.0/(G*dt))*sigmaII0*sigmaII0   );
	*SII = (1.0-phi)*2.0*Z*Eff_strainRate;


}





void Physics_initEta(Physics* Physics, Grid* Grid, BC* BCStokes) {

	int iy, ix, iCell;

		// =======================================================
		// Initial viscosity
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			for (ix = 1; ix<Grid->nxEC-1; ix++) {
				iCell = ix + iy*Grid->nxEC;
				Physics->etaVisc[iCell] = Physics->eta0[iCell];
				Physics->eta[iCell] = Physics->etaVisc[iCell];
				Physics->khi[iCell] = 1E30;
#if (DARCY)
				Physics->khi_b[iCell] = 1E30;
				Physics->eta_b[iCell] 	=  	Physics->eta0[iCell]/Physics->phi[iCell];
#endif
			}
		}
		Physics_copyValuesToSides(Physics->eta, Grid, BCStokes);
		Physics_copyValuesToSides(Physics->khi, Grid, BCStokes);
#if (DARCY)
		Physics_copyValuesToSides(Physics->eta_b, Grid, BCStokes);
		Physics_copyValuesToSides(Physics->khi_b, Grid, BCStokes);
#endif


		int iNode;
		for (iy = 0; iy<Grid->nyS; iy++) {
			for (ix = 0; ix<Grid->nxS; ix++) {
				iNode = ix + iy*Grid->nxS;
				Physics->etaShear[iNode] = shearValue(Physics->eta,  ix   , iy, Grid->nxEC);
				Physics->khiShear[iNode] = shearValue(Physics->khi,  ix   , iy, Grid->nxEC);
			}
		}

}




void Physics_computeEta(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes,MatProps* MatProps)
{
	int iCell, iy, ix;

	//compute* EIIGrid = (compute*) malloc(Grid->nECTot*sizeof(compute));
	//Physics_computeStrainRateInvariant(Physics, Grid, EIIGrid);

	//int C = 0;





	//printf("timeStep = %i, itNonLin = %i\n", Numerics->timeStep, Numerics->itNonLin);

	Physics->dtMaxwellMin = 1E100;
	Physics->dtMaxwellMax = 0;

	compute corr;
	//compute tolerance = 1e-8;
	compute etaVisc0;


	// local copies of values
	compute eta0, etaVisc, eta, cohesion, frictionAngle;
	compute sigma_xx, sigma_xy;
	compute sigma_y, sigmaII;
	compute EII_visc, EII;
	compute phi = 0.0;
	compute n;

	compute phiViscFac = 1.0; // Viscosity correction factor based on porosity

	compute dt = Physics->dt;

//#if (DARCY)
	compute eta_b;
	compute Pe;

	//compute phiMin = Numerics->phiMin;
	compute phiCrit = Numerics->phiCrit;

//#endif
	//compute etaMin = Numerics->etaMin;
	//compute etaMax = Numerics->etaMax;


	compute epsRef = Physics->epsRef;
	//compute dVxdx, dVydx, dVxdy, E_xx, E_xy;
	compute dVydx, dVxdy;

	compute sigmaII0;
	compute Z;
	compute sigma_xx0, sigma_xy0;

	compute sigmaT;//, PeSwitch;
	compute R = 2.0; // radius of the griffith curve

	compute G;

	compute khi;

	compute khi_b, Zb, Py;
	compute B, divV, DeltaP0;

	compute Eff_strainRate;
	//compute Pmin;

	compute tol;

	compute DeltaP;

	//compute sigma_y;
//#pragma omp parallel for private(iy,ix, iCell, sigma_xy, sigma_xx, EII, sigmaII, eta0, etaVisc, n, cohesion, frictionAngle, phi, eta_b, phiViscFac, Pe, sigma_y, etaVisc0, corr, eta) schedule(static,32)
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;




			// Assign local copies
			eta0  			= Physics->eta0 		[iCell];
			etaVisc 		= Physics->eta			[iCell];
			//etaVisc 		= Physics->etaVisc		[iCell];
			n 				= Physics->n			[iCell];
			cohesion 		= Physics->cohesion		[iCell];
			frictionAngle 	= Physics->frictionAngle[iCell];
			G 				= Physics->G[iCell];


#if (DARCY)
			phi = Physics->phi[iCell];
			khi_b = 1E30;
			eta_b = eta0/phi;
#endif



			//  Compute new stresses:
			// xy from node to cell center
			sigma_xy0 = centerValue(Physics->sigma_xy_0,ix,iy,Grid->nxS);
			sigma_xx0 = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];

			sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + (sigma_xy0)*(sigma_xy0));







			// Compute the viscous viscosity
			// ====================================
			// (this section is kinda legacy code, EtaVisc is always equal to eta, now , so I probably don't need to iterate)
			sigma_xy = sigma_xy0;
			sigma_xy += centerValue(Physics->Dsigma_xy_0,ix,iy,Grid->nxS);
			sigma_xx = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];
			sigmaII = sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);

			etaVisc0 = etaVisc;
			corr = 2*etaVisc0; // dummy initial value, just needs to be higher than etaVisc0
			//while (fabs(corr/etaVisc0)>tolerance) {
				EII_visc = sigmaII/(2*etaVisc);
				corr = phiViscFac  *  eta0 * pow(EII_visc/epsRef     ,    1.0/n - 1.0)     -    etaVisc ;
				etaVisc += 1.0*corr;
			//}
			//etaVisc = (etaVisc + etaVisc0)/2;



			eta = etaVisc;




#if (DARCY)
			// Limit the effective pressure
			if (phi>=phiCrit) {
				B = G/sqrt(Physics->phi[iCell]);

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];


				Zb 	= 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(B*dt));






				DeltaP = Zb * ( - divV + DeltaP0/(B*dt) ); // Pc
				Pe =  (1.0-phi) * DeltaP;
				//Pe = Physics->Pc[iCell];
			} else {
				Pe 		= Physics->P [iCell];
			}


#else
			Pe 		= Physics->P [iCell];
#endif






			sigma_y = cohesion * cos(frictionAngle)   +   Pe * sin(frictionAngle);






#if (DARCY)
			sigmaT = cohesion/R;
			tol = 1e-8;
			//Pmin = ((sigmaT-tol) - cohesion * cos(frictionAngle)) / sin(frictionAngle);
			//printf("iCell = %i, Pe = %.2e, Pmin  = %.2e, -sigmaT = %.2e\n", iCell, Pe, Pmin, -sigmaT);


			if (Pe < 0.0) {
				sigma_y = +Pe+(sigmaT-tol); // Pe will be shifted to 0 (arbitrary)
				//printf("A iCell = %i, Pe = %.2e, sigma_y = %.2e\n", iCell, Pe, sigma_y);
			}
			if (Pe<-sigmaT) {
				Pe = -sigmaT;
				sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle

				//printf("B iCell = %i, Pe limited = %.2e, sigma_y = %.2e\n", iCell, Pe, sigma_y);
			}
#else
			// Since there is no griffiths handling for negative pressure for the non darcy case yet
			// here I assume a flat Mohr Coulomb when Pe <0
			if (Pe<0) {
				sigma_y = cohesion * cos(frictionAngle);
			}

			/*
			if (sigma_y<1e-8) {
				Pmin = -cohesion*cos(frictionAngle)/sin(frictionAngle)*0.95;
				if (Pe<Pmin){
					Pe = Pmin;
				}
				sigma_y = 1e-8;
			}
			*/
#endif







			// Update sigmaII according to the current visco-plastic eta
			// ====================================
			//compute khi_old = Physics->khi[iCell];
			khi = 1E30; // first assume that Eps_pl = 0, (therefore the plastic "viscosity" khi is inifinite)
			Z 	= 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));


			dVxdy = ( Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy-1)*Grid->nxVx] +
					Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy-1)*Grid->nxVx] )/4./Grid->dy;


			dVydx = ( Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy] - Physics->Vy[(ix-1)+(iy-1)*Grid->nxVy] +
					Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix-1)+(iy  )*Grid->nxVy] )/4./Grid->dx;

			compute Eps_xy = 0.5*(dVxdy + dVydx);
			//compute Eps_xx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]								 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			compute Eps_xx 	= (Physics->Vx[ix + iy*Grid->nxVx] - Physics->Vx[ix-1 + iy*Grid->nxVx])/Grid->dx;//Grid->DXS[ix-1];


			EII = sqrt(Eps_xx*Eps_xx + Eps_xy*Eps_xy);

			//printf("EII = %.2e, EIbb = %.2e\n", EII, sqrt(Eps_xx*Eps_xx + Eps_xy*Eps_xy));

			Eff_strainRate = sqrt(EII*EII + 1.0*Eps_xx*sigma_xx0/(G*dt) + 1.0*Eps_xy*sigma_xy0/(G*dt) + 1.0/4.0*(1.0/(G*dt))*(1.0/(G*dt))*sigmaII0*sigmaII0   );
			sigmaII = (1.0-phi)*2.0*Z*Eff_strainRate;

			sigma_xx = (Z*(2.0*Eps_xx + sigma_xx0/(G*dt)));
			sigma_xy = (Z*(2.0*Eps_xy + sigma_xy0/(G*dt)));

			if (sigmaII > sigma_y) {
				khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );


				Z 	= 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
				sigmaII = (1.0-phi)*2.0*Z*Eff_strainRate;

			}

			//khi = khi_old + 0.5*(khi-khi_old);

			if (khi<0.0) {
				printf("khi = %.2e, eta = %.2e, G = %.2e, dt = %.2e, Eff_Strainrate = %.2e, 1-phi = %.2e, sigma_y = %.2e, Pe = %.2e, Pmin = %.2e\n", khi, eta, G, dt, Eff_strainRate, 1.0-phi, sigma_y, Pe, -cohesion*cos(frictionAngle)/sin(frictionAngle));
				printf("WTF!\n");
				exit(0);
			}

			// Copy updated values back
			Physics->eta[iCell] = eta;
			Physics->etaVisc[iCell] = etaVisc;// obsolete
			Physics->khi[iCell] = khi;





			// Griffiths


#if (DARCY)

			//compute khi_b_old = Physics->khi_b[iCell];
			khi_b = 1E30;

			// Limit the effective pressure
			Py = sigmaII - sigmaT;
			B = G/sqrt(Physics->phi[iCell]);




			Zb 	= 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(B*dt));



			if (phi>=phiCrit) {
				//compute PeOld = Pe;

				compute DeltaP = Zb * ( - divV + DeltaP0/(B*dt) ); // Pc
				Pe =  (1.0-phi) * DeltaP;
				//printf("iCell = %i, Pe = %.2e\n", iCell, Pe);
				//printf("Pe = %.2e, PeOld = %.2e\n", Pe, PeOld);
				// if sign is opposite, then Pe = 0
				// otherwise khi_b has to be negative in order to make the equation switch side
				// next iteration Pe and Py might have the same sign, and it will be ok







				if (Pe < Py) {
					if (Pe/Py<0) {
						printf("icell = %i, Pe = %.2e, Py = %.2e, sigmaII = %.2e\n", iCell, Pe, Py, sigmaII);
						printf("Pe and Py have opposite sense.\n");
						exit(0);//
						//Py = -1e-5;
					}

					//compute Pe_old = Pe;

					khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(B*dt))   - 1.0/(B*dt) - 1.0/eta_b    );
					Zb 	= 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(B*dt));
					Pe = (1.0-phi) * Zb * ( - divV + DeltaP0/(B*dt) ); // Pc

					//printf("Pe/Py-1.0 = %.2e\n", Pe/Py-1.0);



				}
				Physics->Pc[iCell] = Pe;
			}

			//khi_b = khi_b_old + 0.5*(khi_b - khi_b_old);

			Physics->khi_b[iCell] = khi_b;
			Physics->eta_b[iCell] = eta_b;
#endif








			//printf("%.2e  ", khi);
		}
		//printf("\n");
	}



	Physics_copyValuesToSides(Physics->etaVisc, Grid, BCStokes);
	Physics_copyValuesToSides(Physics->eta, Grid, BCStokes);
	Physics_copyValuesToSides(Physics->khi, Grid, BCStokes);
	//Physics_copyValuesToSides(sigma_xy_EC, Grid, BCStokes);
#if (DARCY)
	Physics_copyValuesToSides(Physics->eta_b, Grid, BCStokes);
	Physics_copyValuesToSides(Physics->khi_b, Grid, BCStokes);
#endif





	// ================================================================================
	// 									Shear nodes viscosity
	int iNode;
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = shearValue(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = shearValue(Physics->khi,  ix   , iy, Grid->nxEC);
		}

	}

	// 									Shear nodes viscosity
	// ================================================================================










}









void Physics_computeEta_applyPlasticity(compute* eta, compute* Pe, compute* phi, compute* cohesion, compute* frictionAngle, compute* EII, compute* sigmaII_phiFac, compute* sigma_xx, compute* sigma_xy)
{


}














void Physics_updateDt(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics)
{

	Physics->dtDarcy = 0.0;
	compute dtOld = Physics->dt;
	/*
	if (fabs(Physics->maxV)<1E-6)
		Physics->maxV = 1E-6;
	 */

	// Although a really large value can be useful to go to the steady state of phi directly
	if (Numerics->timeStep<=0 && Numerics->itNonLin == 0) {
		Physics->maxV = fabs(Physics->epsRef*(Grid->xmax-Grid->xmin)/2.0);
		if (fabs(Physics->maxV)<1E-8)
			Physics->maxV = 1E-8;
	}


	Physics->dtAdv 	= Numerics->CFL_fac*Numerics->dLmin/(Physics->maxV); // note: the min(dx,dy) is the char length, so = 1
	Physics->dtT 	= 10.0*Numerics->CFL_fac*fmin(Grid->dx, Grid->dy)/(3*min(MatProps->k,MatProps->nPhase));
int iCell, iy, ix;
#if (DARCY)
/*
	if (Numerics->timeStep==1 & Numerics->itNonLin == 0) {

	}
	*/
	Physics->dtDarcy 	= 1e10;//10.00*Numerics->CFL_fac*fmin(Grid->dx, Grid->dy)/(3*Physics->minPerm);


	compute CompactionLength;
	compute DarcyVelX, DarcyVelY;
	compute phi;
	compute perm_eta_f;
	compute dPfdx, dPfdy;
	compute CompactionTime;
	//compute CFLtime;
	compute VelSolidX, VelSolidY;
	compute VelFluidX, VelFluidY;
	compute VelFluid;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			phi = Physics->phi[iCell];
			perm_eta_f = Physics->perm_eta_f[iCell];
			dPfdx = (Physics->Pf[ix+1 + iy*Grid->nxEC] - Physics->Pf[ix-1 + iy*Grid->nxEC])/2.0/Grid->dx;
			dPfdy = (Physics->Pf[ix + (iy+1)*Grid->nxEC] - Physics->Pf[ix + (iy-1)*Grid->nxEC])/2.0/Grid->dy;
			CompactionLength = sqrt(4.0/3.0*perm_eta_f * (Physics->eta[iCell]/phi));
			DarcyVelX = perm_eta_f * (-dPfdx + Physics->rho_f_g*Physics->gFac[0]);
			DarcyVelY = perm_eta_f * (-dPfdy + Physics->rho_f_g*Physics->gFac[1]);
			VelSolidX = (Physics->Vx[ix-1+ iy*Grid->nxVx] +  Physics->Vx[ix + iy*Grid->nxVx])/2.0;
			VelSolidY = (Physics->Vx[ix+ (iy-1)*Grid->nxVx] +  Physics->Vx[ix + iy*Grid->nxVx])/2.0;

			VelFluidX = DarcyVelX/phi  + VelSolidX;
			VelFluidY = DarcyVelY/phi  + VelSolidY;

			VelFluid = sqrt(VelFluidX*VelFluidX + VelFluidY*VelFluidY);

			CompactionTime = CompactionLength/VelFluid;
			//CFLtime =Numerics->dLmin/VelFluid;



			//printf("CompactionLength = %.2e, DarcyVel = %.2e, Vx = %.2e, VelFluid = %.2e\n",CompactionLength, (sqrt(DarcyVelX*DarcyVelX + DarcyVelY*DarcyVelY)), Physics->Vx[10], VelFluid);

			if (CompactionTime*Numerics->CFL_fac<Physics->dtDarcy ) {
				Physics->dtDarcy = CompactionTime*Numerics->CFL_fac;
				//printf("Compaction time = %.2e, grid time = %.2e\n", CompactionTime, Grid->dx/(sqrt(DarcyVelX*DarcyVelX + DarcyVelY*DarcyVelY)));
			}

			/*
			if (CFLtime/2.0<Physics->dtDarcy ) {
				Physics->dtDarcy = CFLtime/2.0;
				//printf("Physics->maxV = %.2e, DarcyVel = %.2e\n",Physics->maxV,(sqrt(DarcyVelX*DarcyVelX + DarcyVelY*DarcyVelY)));
				//printf("ix = %i, iy = %i, Compaction time = %.2e, perm0 = %.2e, perm = %.2e, phi = %.2e, phase = %i,CompactionLength = %.2e, CFLtime = %.2e\n", ix, iy, CompactionTime, Physics->perm0[iCell], perm, Physics->phi[iCell], Physics->phase[iCell], CompactionLength, CFLtime);
			}
			*/





		}
	}



#endif
	Physics->dtMaxwellMin = 1E+100;
	Physics->dtMaxwellMax = 1E-100;
	compute EII, dtElastic;
	compute eta;
	compute sigma_xy0, sigma_xx0, sigmaII0;
	//compute DSigmaMax = 0.1; // DeltaSigma/SigmaII0 (maximum change relative to the previous value)
	compute DSigmaMax = 1.0; // DeltaSigma (maximum change relative to the characteristic stress)
	compute this_dt;
	compute Dsigma_xx,Dsigma_xy, DsigmaII;
	compute dtMaxwell;
	dtElastic = 1E10;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;

			Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix,iy, &EII);
			eta = Physics->eta[iCell];
			sigma_xy0  = Physics->sigma_xy_0[ix-1 + (iy-1)*Grid->nxS];// + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
			sigma_xy0 += Physics->sigma_xy_0[ix   + (iy-1)*Grid->nxS];// + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
			sigma_xy0 += Physics->sigma_xy_0[ix-1 + (iy  )*Grid->nxS];// + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
			sigma_xy0 += Physics->sigma_xy_0[ix   + (iy  )*Grid->nxS];// + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
			sigma_xy0 /= 4.0;

			sigma_xx0 = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + (sigma_xy0)*(sigma_xy0));

			Dsigma_xy  = Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];// + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
			Dsigma_xy += Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];// + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
			Dsigma_xy += Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];// + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
			Dsigma_xy += Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];// + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
			Dsigma_xy /= 4.0;

			Dsigma_xx = Physics->Dsigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			DsigmaII = sqrt((Dsigma_xx)*(Dsigma_xx)    + (Dsigma_xy)*(Dsigma_xy));


			if (sigmaII0 == 0) {
				sigmaII0 =0.0;
			}
			//this_dt = fabs(eta / ( ( (2*eta*EII/sigmaII0 - 1)/DSigmaMax -1 )*Physics->G[iCell] )); // relative dsigma
			if (DsigmaII>DSigmaMax) {
				this_dt = fabs(eta / ( ( (2.0*eta*EII - sigmaII0)/DSigmaMax -1.0 )*Physics->G[iCell] )    ); // Absolute Dsigma

				if (this_dt<dtElastic) {
					dtElastic = this_dt;
					//compute estDsigma = (2.0*eta*EII - sigmaII0) * (Physics->G[iCell]*this_dt/(eta + Physics->G[iCell]*this_dt));
					//printf("DsigmaII = %.2e, Dsigma_xx = %.2e, dtElastic = %.2e, sigmaII0 = %.2e, sigma_xx0 = %.2e, estDsigma = %.2e\n",DsigmaII, Physics->Dsigma_xx_0[iCell], dtElastic, sigmaII0, Physics->sigma_xx_0[iCell],estDsigma);
				}
			}

			if (Physics->phase[iCell]!=Physics->phaseAir && Physics->phase[iCell]!=Physics->phaseWater) {
				dtMaxwell = eta/Physics->G[iCell];
				if (dtMaxwell > Physics->dtMaxwellMax) {
					Physics->dtMaxwellMax = dtMaxwell;
				}
				if (dtMaxwell < Physics->dtMaxwellMin) {
					Physics->dtMaxwellMin = dtMaxwell;
				}
			}
		}
	}


	//printf("dtAdv = %.2e, dtT = %.2e, Numerics->dLmin = %.2e, Numerics->CFL = %.2e, (Physics->maxV) = %.2e, dtElastic = %.2e\n", Physics->dtAdv, Physics->dtT, Numerics->dLmin, Numerics->CFL_fac, (Physics->maxV), dtElastic);


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

	if (MatProps->G[Physics->phaseRef] < 1E10) { // to enable switching off the elasticity
		Physics->dt  =  fmin(Physics->dt,dtElastic);
	}

	//Physics->dtAdv 	= fmin(Physics->dt,Physics->dtAdv);
	//Physics->dtT 	= fmin(Physics->dtT,Physics->dtAdv);
#if (HEAT)
	Physics->dt  =  fmin(Physics->dtT,Physics->dtAdv);
#endif

#if (DARCY)
	Physics->dt  =  fmin(Physics->dt,Physics->dtDarcy);
#endif
	//printf("A0- Physics->dt = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, Physics->dtAdv = %.2e, Physics->dtT = %.2e, Physics->dtDarcy = %.2e\n", Physics->dt, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Physics->dtAdv, Physics->dtT, Physics->dtDarcy);

	compute corr;
	if (Numerics->timeStep<=0){// && Numerics->itNonLin==0) {
		Physics->dt = 1E-10;
	} else if (Numerics->timeStep==1){// && Numerics->itNonLin==0) {
		//Physics->dt = Physics->dt;
	} else {
		corr = (Physics->dt-dtOld);
		if (corr>dtOld) {
			corr = dtOld;
		} else if (corr< -dtOld) {
			corr = -dtOld;
		}
		Physics->dt = dtOld + 0.1*corr;
	}


	//printf("A - Physics->dt = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, Physics->dtAdv = %.2e, Physics->dtT = %.2e, Physics->dtDarcy = %.2e\n", Physics->dt, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Physics->dtAdv, Physics->dtT, Physics->dtDarcy);



	//Physics->dtAdv = Physics->dt;
	//Physics->dtT = Physics->dt;

#if (DARCY)
	Physics->dtDarcy = Physics->dt;
#endif

	//Numerics->dtMax  = 1e-2;
	//Numerics->dtMin = Physics->dtMaxwellMin + 0.2*(Physics->dtMaxwellMax - Physics->dtMaxwellMin);
	//Numerics->dtMax = Physics->dtMaxwellMax - 0.2*(Physics->dtMaxwellMax - Physics->dtMaxwellMin);

	if (MatProps->G[Physics->phaseRef] < 1E10) { // to enable switching off the elasticity
		Numerics->dtMin = pow(10,log10(Physics->dtMaxwellMin) + 0.1*(log10(Physics->dtMaxwellMax) - log10(Physics->dtMaxwellMin) ));
		Numerics->dtMax = pow(10,log10(Physics->dtMaxwellMax) - 0.0*(log10(Physics->dtMaxwellMax) - log10(Physics->dtMaxwellMin) ));
	}


	if (Physics->dt<Numerics->dtMin) {
		Physics->dt = Numerics->dtMin;
	} else if (Physics->dt>Numerics->dtMax) {
		Physics->dt = Numerics->dtMax;
	}
	/*
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
*/

	printf("B - Physics->dt = %.2e, dtAdv = %.2e, dtT = %.2e, PdtDarcy = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, dtMin = %.2e, dtMax = %.2e\n", Physics->dt, Physics->dtAdv, Physics->dtT, Physics->dtDarcy, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Numerics->dtMin, Numerics->dtMax);





}


#if (DARCY)
void Physics_computePerm(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes)
{
	Physics->minPerm = 1E100;
	int iy, ix;
	int iCell;
	compute phi;
	compute phiRef = 0.0001;
	compute PermEffRef = Physics->perm0_eta_f[0]  *  phiRef*phiRef*phiRef  / ( (1.0-phiRef)*(1.0-phiRef));
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
			if (Physics->phase[iCell] != Physics->phaseAir && Physics->phase[iCell] != Physics->phaseWater) {
				Physics->perm_eta_f[iCell] = Physics->perm0_eta_f[iCell]  *  phi*phi*phi  * ( (1.0-phi)*(1.0-phi));
			} else {
				Physics->perm_eta_f[iCell]=1e6*PermEffRef;
			}

			/*
			if (Physics->perm[iCell]<Physics->minPerm) {
				Physics->minPerm = Physics->perm[iCell];
			}
			*/


			/*
			//printf("Physics->perm[iCell] = %.2e, PermRef = %.2e, Physics->eta_f = %.2e\n",Physics->perm[iCell], PermEffRef, Physics->eta_f);
			if (Physics->perm[iCell]>1e6*PermEffRef) {
				Physics->perm[iCell] = 1e6*PermEffRef;
			} else if (Physics->perm[iCell]<1e-3*PermEffRef) {
				Physics->perm[iCell] = 1e-3*PermEffRef;
			}
			*/


		}
	}

	Physics_copyValuesToSides(Physics->perm_eta_f, Grid, BCStokes);

	if (DEBUG) {
		printf("=== Check perm  ===\n");
		int C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->perm_eta_f[C]);
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
	compute maxDiv = 0;
	compute maxPhi = 0;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			dx = Grid->DXS[ix-1];
			dy = Grid->DYS[iy-1];
			divV  = (  Physics->Vx[ix+iy*nxVx] - Physics->Vx[ix-1+ iy   *nxVx]  )/dx;
			divV += (  Physics->Vy[ix+iy*nxVy] - Physics->Vy[ix  +(iy-1)*nxVy]  )/dy;


			//printf("divV[%i, %i] = %.2e, dy = %.2e, dx = %.2e\n",iy, ix, divV, dy, dx);
			//                    Physics->phi[iCell] = 1 - exp(-divV*dt)*(1-Physics->phi0[iCell]);
			// Dphi/Dt = (1-phi)*divV

			//if (Physics->phase[iCell] != Physics->phaseAir && Physics->phase[iCell] != Physics->phaseWater) {
			Physics->phi[iCell] = Physics->phi0[iCell] + dt*0.5*(    (1.0-Physics->phi0[iCell])*Physics->divV0[iCell] + (1.0-Physics->phi[iCell])*divV   );
			//Physics->phi[iCell] = Physics->phi0[iCell] + dt*(  (1.0-Physics->phi[iCell])*divV   );


			//Physics->phi[iCell] = Physics->phi0[iCell] + dt*(    (1.0-Physics->phi[iCell])*divV   );
			//}
			/*else {
				Physics->phi[iCell] = Numerics->phiMax;
			}
			*/

			if (Physics->phi[iCell] > Numerics->phiMax) {
				Physics->phi[iCell] = Numerics->phiMax;
			} else if (Physics->phi[iCell] < Numerics->phiMin) {
				Physics->phi[iCell] = Numerics->phiMin;
			}




			Physics->Dphi[iCell] = Physics->phi[iCell] - Physics->phi0[iCell];

			if (fabs(divV)>maxDiv) {
				maxDiv = fabs(divV);
			}

			if (fabs(Physics->phi[iCell])>maxPhi) {
				maxPhi = fabs(Physics->phi[iCell]);
			}



			/*
			if (iCell == 150) {
			printf("    divV = %.2e, phi0 = %.2e, phi = %.2e, dt = %.2e\n", divV, Physics->phi0[iCell], Physics->phi[iCell], dt);
			}
			*/
			sum += Physics->phi[iCell];
		}
		//printf("divV = %.2e\n",divV);
	}
	//printf("                    sum phi = %.2e\n", sum);


	Physics_copyValuesToSides(Physics->phi, Grid, BCStokes);
	Physics_copyValuesToSides(Physics->Dphi, Grid, BCStokes);


	//printf("maxDiv = %.2e, maxPhi = %.2e\n", maxDiv, maxPhi);



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

	Numerics->phiMin = 1e-6;
	Numerics->phiCrit = 0.0001; // i.e. value above which Pe = Pc
	Numerics->phiMax = 0.8;

	printf("in InitPhi\n");
	int type = 0; // 0, porosity wave; 1, with ocean

	if (type==0) {
		compute xc = Grid->xmin + (Grid->xmax - Grid->xmin)/2.0;
		compute yc = Grid->ymin + (Grid->ymax - Grid->ymin)/6.0;
		printf("xc = %.2e, yc = %.2e", xc, yc);

		//compute xc = Grid->xmax - (Grid->xmax - Grid->xmin)/25.0;
		//compute yc = Grid->ymin + (Grid->ymax - Grid->ymin)/12.0;
		compute phiBackground = 0.1;//Numerics->phiMin;
		compute A = 00.0*phiBackground;
		compute x = Grid->xmin-Grid->DXEC[0]/2.0;
		compute y = Grid->ymin-Grid->DYEC[0]/2.0;
		compute w = 2.0;//(Grid->xmax - Grid->xmin)/15.0;
		compute XFac = 1.0;
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





				if (Physics->phase[iCell] == Physics->phaseAir || Physics->phase[iCell] == Physics->phaseWater) {
					Physics->phi [iCell] = Numerics->phiMax;
				}

				if (Physics->phase[iCell] == 2) {
					Physics->phi [iCell] = Numerics->phiMin;
				}






				Physics->Dphi  [iCell]  = Physics->phi[iCell];

				//Physics->phi0  [iCell] = Physics->phi[iCell];
				if (ix<Grid->nxEC-1) {
					x += Grid->DXEC[ix];
				} else {
					x = Grid->xmin-Grid->DXEC[0]/2.0;
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
		if (Grid->isPeriodic) {
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
		if (Grid->isPeriodic) {
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


	if (Grid->isPeriodic) {
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
	//printf("end neighbour stuff");
}

void Physics_copyValuesToSidesi(int* ECValues, Grid* Grid, BC* BC)
{

	// Replace boundary values by their neighbours
	int INeigh, iy, ix, I;

	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (Grid->isPeriodic) {
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
		if (Grid->isPeriodic) {
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


	if (Grid->isPeriodic) {
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
//	printf("end neighbour stuff");
}









void Physics_computeRho(Physics* Physics, Grid* Grid)
{

	int iCell;

	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->rho_g[iCell] = Physics->rho0_g[iCell];

#if (DARCY)
		//Physics->rho[iCell] = Physics->rho0[iCell];
		Physics->rho_g[iCell] = (1.0 - Physics->phi[iCell])*Physics->rho0_g[iCell] + Physics->phi[iCell]*Physics->rho_f_g;
		//Physics->rho[iCell] = Physics->rho0[iCell] * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell]);
#endif

	}

if (DEBUG) {
#if (DARCY)
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
	printf("=== Check rho_g  ===\n");
	C = 0;
	//int iy, ix;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.8e  ", Physics->rho_g[C]);
			C++;
		}
		printf("\n");
	}

#endif
}

}




void Physics_get_ECVal_FromSolution (compute* Val, int ISub, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
// Where Val is the value to extract from the solution, and DVal the increment since the last time step, IStep is the index of the subsystem of equations
	int I, IBC, INeigh, iy, ix;
	int INumMap0 = Numbering->subEqSystem0Dir[ISub];
	//printf("eq0 = %i, ISub = %i\n", INumMap0, ISub);
	int iCell;
/*
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
				else if (BC->type[IBC]==Dirichlet) {
					Val[iCell] = BC->value[IBC];
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
	*/

	compute scale;

#pragma omp parallel for private(iy, ix, I, iCell, IBC, INeigh, scale) schedule(static,32)
	for (iy = 0; iy<Grid->nyEC; iy++) {
		for (ix = 0; ix<Grid->nxEC; ix++) {
			iCell = ix + iy*Grid->nxEC;
			I = Numbering->map[iCell + INumMap0];


			//printf("I = %i \n", I);
			if (I>=0) {
				scale = 1.0;//EqSystem->S[I];
				Val[iCell]  = EqSystem->x[I]*scale;
			}
			else {

				IBC = abs(I)-1; // BC nodes are numbered -1 to -n

				// Get neighbours index
				if (iy==0) { // lower boundary
					if (Grid->isPeriodic){
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
					if (Grid->isPeriodic){
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


				scale = 1.0;//EqSystem->S[INeigh];

				if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Val[iCell] = 2.0*BC->value[IBC] - EqSystem->x[INeigh]*scale;
				//	printf("IBC %i is Dir Ghost\n",IBC);
				}
				else if (BC->type[IBC]==NeumannGhost) { // Neumann
					if (ix==0)  {// left or bottom boundary
						Val[iCell] = EqSystem->x[INeigh]*scale - BC->value[IBC]*Grid->DXEC[0];
					} else if (ix==Grid->nxEC-1) {
						Val[iCell] = EqSystem->x[INeigh]*scale + BC->value[IBC]*Grid->DXEC[Grid->nxEC-2];
					}
					if (iy==0) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh]*scale - BC->value[IBC]*Grid->DYEC[0];
					} else if (iy==Grid->nyEC-1) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh]*scale + BC->value[IBC]*Grid->DYEC[Grid->nyEC-2];
					}
				}
				else if (BC->type[IBC]==Dirichlet) {
					Val[iCell] = BC->value[IBC];
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
	//coord depth, y;

	SingleParticle* thisParticle;
	//compute locX, locY;
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
					maxContrib = contribPhase[iPhase];
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
	//int ixStart, ixEnd, ixInc;
	//int iyStart, iyEnd, iyInc;

	int C;


//printf("enter Plitho\n");

	// Contribution of gy
	if (Physics->g[1]>0){
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				iCell = ix + iy*Grid->nxEC;
				iCellS = ix + (iy-1)*Grid->nxEC;
				if (iy==0) {
					rho_g_h = Physics->rho_g[iCell] * Physics->gFac[1] * (-0.5*Grid->DYEC[iy] );
				} else {
					rho_g_h += 0.5*(Physics->rho_g[iCell]+Physics->rho_g[iCellS]) * Physics->gFac[1] * Grid->DYEC[iy-1] ;
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
					rho_g_h = Physics->rho_g[iCell] * -Physics->gFac[1] * (-0.5*Grid->DYEC[iy-1] );
				} else {
					rho_g_h += 0.5*(Physics->rho_g[iCell]+Physics->rho_g[iCellN]) * -Physics->gFac[1] * Grid->DYEC[iy] ;
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
						rho_g_h = Physics->rho_g[iCell] * Physics->gFac[0] * (-0.5*Grid->DXEC[ix] );
					} else {
						rho_g_h += 0.5*(Physics->rho_g[iCell]+Physics->rho_g[iCellW]) * Physics->gFac[0] * Grid->DXEC[ix-1] ;
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
						rho_g_h = Physics->rho_g[iCell] * -Physics->gFac[0] * (-0.5*Grid->DXEC[ix-1] );
					} else {
						rho_g_h += 0.5*(Physics->rho_g[iCell]+Physics->rho_g[iCellE]) * -Physics->gFac[0] * Grid->DXEC[ix] ;
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







