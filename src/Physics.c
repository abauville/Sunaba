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
	Physics->dt = 1.0e-100;


	Physics->phaseListHead 	= (SinglePhase**) malloc( Grid->nECTot 		* sizeof(  SinglePhase*  ) ); // array of pointers to particles
	for (i=0;i<Grid->nECTot;i++) {
		Physics->phaseListHead[i] = (SinglePhase*) malloc( 1 		* sizeof(  SinglePhase  ) );
		Physics->phaseListHead[i]->phase = -1;
		Physics->phaseListHead[i]->weight = 0;
		Physics->phaseListHead[i]->next = NULL;
	}
	Physics->sumOfWeightsCells = (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->sumOfWeightsNodes = (compute*) 	malloc( Grid->nSTot  * sizeof(compute) );

	Physics->Vx 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
	Physics->P 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->Z 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->eta 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->khi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	//Physics->etaVisc		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	//Physics->eta0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	//Physics->n 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->rho_g 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	//Physics->rho0_g 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );



#if (HEAT)

	Physics->k 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->DT 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
#endif

	//Physics->Plitho 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

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


	//Physics->cohesion 		= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	//Physics->frictionAngle 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );


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
		Physics->T[i]  = 1.0;
		Physics->DT[i] = 0.0;
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


void Physics_freeMemory(Physics* Physics, Grid* Grid)
{

	// Free phaseList
	int iCell;
	SinglePhase* temp;
	for (iCell=0;iCell<Grid->nECTot;iCell++) {
		while (Physics->phaseListHead[iCell] != NULL)
		{
			temp = Physics->phaseListHead[iCell];
			Physics->phaseListHead[iCell] = Physics->phaseListHead[iCell]->next;
			free(temp);
		}
	}
	free( Physics->phaseListHead );

	free(Physics->Z);

	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P );

	free( Physics->eta );

	//free( Physics->eta0 );

	//free(Physics->etaVisc);

	free(Physics->etaShear);
	free( Physics->khi );
	free( Physics->khiShear );

	//free( Physics->n );
	free( Physics->rho_g );

	//free( Physics->rho0_g );



#if (HEAT)
	free( Physics->k );
	free(Physics->T );
	free(Physics->T0);
	free(Physics->DT );
#endif

	free(Physics->phase);

	free(Physics->G );

	free(Physics->sigma_xx_0 );
	free(Physics->sigma_xy_0 );
	free(Physics->Dsigma_xx_0 );
	free(Physics->Dsigma_xy_0 );

	//free(Physics->cohesion);
	//free(Physics->frictionAngle);

	//free(Physics->Plitho);
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


	free(Physics->sumOfWeightsCells);
	free(Physics->sumOfWeightsNodes);



}





void addSinglePhase(SinglePhase** pointerToHead, int phase)
{
	// Adds a Particle at the beginning of a linked list
	SinglePhase* thisPhase = (SinglePhase*) malloc(sizeof(SinglePhase));
	thisPhase->phase = phase;
	thisPhase->weight = 0.0;

	thisPhase->next = NULL;
	if (*pointerToHead != NULL) {
		thisPhase->next = *pointerToHead;
	}
	*pointerToHead = thisPhase;
}






void Physics_initPToLithostatic(Physics* Physics, Grid* Grid)
{

	int iy, ix, iCell, iCellS, iCellN, iCellW, iCellE;
	compute rho_g_h;
	//int ixStart, ixEnd, ixInc;
	//int iyStart, iyEnd, iyInc;



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
				Physics->P[iCell] = rho_g_h;
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

				Physics->P[iCell] = rho_g_h;
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
					Physics->P[iCell] += rho_g_h;
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
					Physics->P[iCell] += rho_g_h;
				}
			}
		}
	}






	//Physics_computePlitho(Physics, Grid);
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {




#if (DARCY)
		Physics->Pf[iCell] = Physics->P[iCell];
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
	int iCell, iNode;
	//int nNeighbours = 4;
	coord locX, locY;


	// Reinitialize Physics array
	bool* changedHead = (bool*) malloc(Grid->nECTot * sizeof(bool));

#pragma omp parallel for private(iCell) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->sigma_xx_0 [iCell] = 0;
		Physics->sumOfWeightsCells [iCell] = 0;
#if (HEAT)

		Physics->T[iCell] = 0;
#endif
#if (DARCY)
		Physics->DeltaP0[iCell] 		= 0;
		Physics->phi0[iCell] 		= 0;
#endif

		changedHead[iCell] = false;

	}


#pragma omp parallel for private(iNode) schedule(static,32)
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		Physics->sigma_xy_0 [iNode] = 0;
		Physics->sumOfWeightsNodes [iNode] = 0;
	}


	Physics_reinitPhaseList(Physics, Grid);




	printf("Init ok\n");


	compute weight;


	int phase;

	int nxEC = Grid->nxEC;
	compute xMod[4], yMod[4];
	int ix,  iy;

	int C;

	xMod[0] = -1; yMod[0] = -1;
	xMod[1] =  1; yMod[1] = -1;
	xMod[2] = -1; yMod[2] =  1;
	xMod[3] =  1; yMod[3] =  1;



	// Index of neighbouring cells, with respect to the node ix, iy
	int i;
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right


	SingleParticle* thisParticle = NULL;




	int iColor; // indexing of the color group for nodes. Nodes of the same color don't collide with each other. i.e. similar to matrix coloring
	int ixStart[4] = {0,0,1,1};
	int iyStart[4] = {0,1,0,1};
	SinglePhase* thisPhaseInfo;
	printf("I'm in\n");
	compute sum54 = 0.0;
	compute sum54_0 = 0.0;
	compute sum54_1 = 0.0;
	compute sum55 = 0.0;;
	compute sum55_0 = 0.0;
	compute sum55_1 = 0.0;
	for (iColor = 0; iColor < 4; ++iColor) {
		//#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight) schedule(static,32)
		for (iy = iyStart[iColor]; iy < Grid->nyS; iy+=2) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ix+=2) { // I don't get why though
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

					for (i=0; i<4; i++) {


						iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
						weight = fabs((locX + xMod[i]*1.0)   *   (locY + yMod[i]*1.0));


						/*
						if (iCell == 4+5*Grid->nxEC) {
							sum54 += weight;
							if (phase==0) {
								sum54_0 += weight;
							} else {
								sum54_1 += weight;
							}
						}
						if (iCell == 5+5*Grid->nxEC) {
							sum55 += weight;
							if (phase==0) {
								sum55_0 += weight;
							} else {
								sum55_1 += weight;
							}
						}
						 */


						// Get the phase and weight of phase contribution for each cell
						thisPhaseInfo = Physics->phaseListHead[iCell];
						bool oldchanged = changedHead[iCell];
						bool added = false;
						while (thisPhaseInfo->phase != phase) {
							if (thisPhaseInfo->next == NULL) {
								//thisPhaseInfo->phase = phase;


								if (!changedHead[iCell]) {
									//printf("koko\n");
									thisPhaseInfo->phase = phase;
									changedHead[iCell] = true;
									break;
								} else {
									//printf("asoko\n");
									//thisPhaseInfo->phase = phase;
									//printf("koko\n");
									addSinglePhase(&Physics->phaseListHead[iCell],phase);
									thisPhaseInfo = Physics->phaseListHead[iCell];
									added = true;
									break;
									//printf("soko\n");

								}





							} else {
								thisPhaseInfo = thisPhaseInfo->next;
							}
						}
						thisPhaseInfo->weight += weight;

						/*
						//if (iCell == 4+5*Grid->nxEC || iCell == 5+5*Grid->nxEC) {
						if (iCell == 5+5*Grid->nxEC && phase == 0) {
							//printf("iCell = %i, iNode = %i, locX = %.2e, locY = %.2e, phase = %i, weight = %.2e, sum55 = %.2e, sum55_0 = %.2e, thisPhaseInfoWeight = %.2e\n",iCell, iNode, locX,locY,phase,weight,sum55,sum55_0, thisPhaseInfo->weight);
							printf("iCell = %i, iNode = %i, locX = %.2e, locY = %.2e, phase = %i, weight = %.2e, sum55 = %.2e, sum55_0 = %.2e, thisPhaseInfoWeight = %.2e, oldchanged = %i, changedHead = %i, added = %i\n",iCell, iNode, locX,locY,phase,weight,sum55,sum55_0, thisPhaseInfo->weight, oldchanged , changedHead[iCell], added);
						}
						 */


						// For properties that are stored on the markers, sum contributions
						Physics->sigma_xx_0		[iCell] += thisParticle->sigma_xx_0 * weight;
#if (HEAT)
						Physics->T				[iCell] += thisParticle->T * weight;
#endif
#if (DARCY)
						Physics->DeltaP0		[iCell] += thisParticle->DeltaP0 * weight;
						Physics->phi0			[iCell] += thisParticle->phi * weight;
#endif
						Physics->sumOfWeightsCells	[iCell] += weight;

					}
					thisParticle = thisParticle->next;
				}
			}
		}
	}

	printf("I'm out\n");

	// Copy contribution from one side to the other in case of periodic BC
	if(Grid->isPeriodic) {
		int iCellS, iCellD, j;
#pragma omp parallel for private(iy, j, iCellS, iCellD) schedule(static,32)
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (j = 0; j<2; ++j) {
				iCellS = j + iy*Grid->nxEC;
				iCellD = Grid->nxEC-2+j + iy*Grid->nxEC;

				Physics->sigma_xx_0		[iCellD] += Physics->sigma_xx_0		[iCellS];
				Physics->sigma_xx_0		[iCellS]  = Physics->sigma_xx_0		[iCellD];

#if (HEAT)
				Physics->T				[iCellD] += Physics->T				[iCellS];
				Physics->T				[iCellS]  = Physics->T				[iCellD];
#endif
#if (DARCY)
				Physics->DeltaP0		[iCellD] += Physics->DeltaP0		[iCellS];
				Physics->DeltaP0		[iCellS]  = Physics->DeltaP0		[iCellD];
				Physics->phi0			[iCellD] += Physics->phi0			[iCellS];
				Physics->phi0			[iCellS]  = Physics->phi0			[iCellD];
#endif
				Physics->sumOfWeightsCells	[iCellD] += Physics->sumOfWeightsCells	[iCellS];
				Physics->sumOfWeightsCells	[iCellS]  = Physics->sumOfWeightsCells	[iCellD];


			}
		}
	}


	free(changedHead);

	// Dividing by the sum of weights
#pragma omp parallel for private(iCell) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->sigma_xx_0	[iCell] /= Physics->sumOfWeightsCells	[iCell];
#if (HEAT)
		Physics->T			[iCell] /= Physics->sumOfWeightsCells	[iCell];
#endif
#if (DARCY)
		Physics->DeltaP0	[iCell] /= Physics->sumOfWeightsCells	[iCell];
		Physics->phi0		[iCell] /= Physics->sumOfWeightsCells	[iCell];
#endif
	}





	// Filling side values
	Physics_copyValuesToSides(Physics->sigma_xx_0, Grid);
#if (HEAT)
	Physics_getValuesToSidesFromBC(Physics->T, Grid, BCThermal, NumThermal);
#endif
#if (DARCY)
	Physics_copyValuesToSides(Physics->DeltaP0, Grid);
	Physics_copyValuesToSides(Physics->phi0, Grid);
#endif
	//Physics_copyValuesToSidesi(Physics->sumOfWeights, Grid* Grid);




#if (HEAT)
	// Should probably be moved to a specific function
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->k[iCell] = 0.0;
		thisPhaseInfo = Physics->phaseListHead[iCell];
		while (thisPhaseInfo != NULL) {
			Physics->k[iCell] += MatProps->k[thisPhaseInfo->phase] * thisPhaseInfo->weight;
			thisPhaseInfo = thisPhaseInfo->next;
		}
		Physics->k[iCell] /= Physics->sumOfWeightsCells[iCell];
	}
#endif










	// ==================================
	// Interpolate to nodes
	// ==================================
	int signX, signY, iNodeNeigh;
	xMod[0] =  1; yMod[0] =  1;
	xMod[1] =  0; yMod[1] =  1;
	xMod[2] =  1; yMod[2] =  0;
	xMod[3] =  0; yMod[3] =  0;

	for (iColor = 0; iColor < 4; ++iColor) {
		//#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, weight, signX, signY, iNodeNeigh) schedule(static,32)
		for (iy = iyStart[iColor]; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ++ix) { // I don't get why though
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];

				// Loop through the particles in the shifted cell
				// ======================================
				while (thisParticle!=NULL) {
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

						weight = (locX + xMod[i]*1.0)   *   (locY + yMod[i]*1.0);

						Physics->sigma_xy_0 		[iNodeNeigh] += thisParticle->sigma_xy_0 * weight;
						Physics->sumOfWeightsNodes	[iNodeNeigh] += weight; // using the same arrays

					}
					thisParticle = thisParticle->next;
				}
			}
		}
	}


	// Adding contribution to the other side for periodic BC
	if(Grid->isPeriodic) {
		int iCellS, iCellD;
#pragma omp parallel for private(iy, iCellS, iCellD,i) schedule(static,32)
		for (iy = 0; iy < Grid->nyS; ++iy) {
			iCellS = 0 + iy*Grid->nxS; // Source
			iCellD = Grid->nxS-1 + iy*Grid->nxS; // destination

			Physics->sigma_xy_0 		[iCellD] += Physics->sigma_xy_0  [iCellS];
			Physics->sigma_xy_0 		[iCellS]  = Physics->sigma_xy_0  [iCellD];
			Physics->sumOfWeightsNodes	[iCellD] += Physics->sumOfWeightsNodes[iCellS];
			Physics->sumOfWeightsNodes	[iCellS]  = Physics->sumOfWeightsNodes[iCellD];
		}
	}




	// Dividing by the sum of weights
#pragma omp parallel for private(iNode) schedule(static,32)
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		Physics->sigma_xy_0 [iNode] /= Physics->sumOfWeightsNodes[iNode]; // harmonic average
	}









	if (DEBUG) {


#if (HEAT)

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
		printf("=== Check phi0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", Physics->phi0[C]);
				C++;
			}
			printf("\n");
		}


#endif

	}




}


















#if (HEAT)
void Physics_interpTempFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps, BC* BCThermal)
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


				rhoParticle = MatProps->rho0[phase];// * (1+MatProps->beta[phase]*PFromNodes) * (1-MatProps->alpha[phase]*thisParticle->T);
				if (rhoParticle<0) {
					printf("error: Negative density on Particles in Physisc_interTempFromCellsParticle\n");
					exit(0);
				}

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

				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->T = BCThermal->TT;
				}

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
					if (iy==0) { // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVx  ];
					} else if (iy==Grid->nyVx-1) { // lower boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVx  ];
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet.");
						exit(0);
					}

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
					if (ix==0) {  // lower boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVy + Grid->nVxTot ];
					} else if (ix==Grid->nxVy-1) { // lower boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVy + Grid->nVxTot  ];
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet.");
						exit(0);
					}
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





	/*
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
	 */






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
	compute RefPressure = Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];//Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];
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


	//compute RefPressure = Physics->Pf[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];
	/*
	for (ix = 0; ix < Grid->nxEC; ++ix) {
		iCell = ix + (Grid->nyEC-2)*Grid->nxEC;
		RefPressure += Physics->Pf[iCell];
	}
	RefPressure /= Grid->nxEC;
	 */


	compute RefPressure = 0.0;//Physics->Pf[1 + (Grid->nyEC-2)*Grid->nxEC];
	for (iy = 0; iy < Grid->nyEC-1; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			Physics->Pf [iCell] 	= Physics->Pf [iCell] - RefPressure;
		}
	}

	RefPressure = 0.0;//Physics->Pc[1 + (Grid->nyEC-2)*Grid->nxEC];
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->Pc [iCell] 	= Physics->Pc [iCell] - RefPressure;
	}





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

	//}






#endif


	/*

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

	 */







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


	/*
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

	 */

}
#endif





void Physics_computeStressChanges(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes)
{

	// see Taras' book p. 186
	int ix, iy, iCell, iNode;
	compute Z;
	compute Eps_xx, Eps_xy;
	compute dVxdy, dVydx, dVxdx, dVydy;
	compute G;//, eta, khi;
	//compute phi;
	// compute stress
	//#pragma omp parallel for private(iy, ix, iCell, Eps_xx, Z) schedule(static,32)
	compute dt = Physics->dt;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
								 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
								 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;


			Eps_xx = 0.5*(dVxdx-dVydy);


			Physics->Z[iCell] = 1.0/(1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*dt));

			//Physics->Dsigma_xx_0[iCell] = ( 2.0*Physics->eta[iCell] * Eps_xx  -  Physics->sigma_xx_0[iCell] ) * Z;
			Physics->Dsigma_xx_0[iCell] = Physics->Z[iCell]*(2.0*Eps_xx + Physics->sigma_xx_0[iCell]/(Physics->G[iCell]*dt)) - Physics->sigma_xx_0[iCell];

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
			/*
			eta 	= Physics->etaShear[iNode];
			khi 	= Physics->khiShear[iNode];

			Z = 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
			 */
			Z 	 	= shearValue(Physics->Z, ix, iy, Grid->nxEC);

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


#if (DARCY)
	compute phi;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		phi = Physics->phi[iCell];
		Physics->DDeltaP[iCell] = Physics->Pc[iCell]/(1.0-phi) - Physics->DeltaP0[iCell];
	}
#endif




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
	int iNode, Ix, Iy;
	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};


	compute ShearComp_sqr;


	//#pragma omp parallel for private(ix,iy, IE, dVxdx, dVydy, dVxdy, dVydx) schedule(static,32)
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			IE = ix+iy*Grid->nxEC;
			//I = (ix-1)+(iy-1)*Grid->nxC;

			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
								 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
								 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			// Method A: using the averageing of derivatives on the four nodes
			// Compute Eps_xy at the four nodes of the cell
			// 1. Sum contributions
			dVxdy = 0;
			dVydx = 0;
			ShearComp_sqr = 0.0;
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];

				dVxdy = ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
									  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


				dVydx = ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
									  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;
				//printf("koko\n");
				ShearComp_sqr += (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx)) ;

			}
			StrainRateInvariant[IE] = sqrt(  (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy))  +  0.5*ShearComp_sqr );


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

	compute ShearComp_sqr;
	int iNode, Ix, Iy;
	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};				dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
																 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
	dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

	// Method A: using the averageing of derivatives on the four nodes
	// Compute Eps_xy at the four nodes of the cell
	// 1. Sum contributions
	dVxdy = 0;
	dVydx = 0;
	ShearComp_sqr = 0.0;
	for (iNode = 0; iNode < 4; ++iNode) {
		Ix = (ix-1)+IxMod[iNode];
		Iy = (iy-1)+IyMod[iNode];

		dVxdy = ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
							  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


		dVydx = ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
							  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;
		//printf("koko\n");
		ShearComp_sqr += (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx)) ;

	}
	*EII = sqrt(  (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy))  +  0.5*ShearComp_sqr );


}


void Physics_computeStrainRateInvariantForOneNode(Physics* Physics, BC* BCStokes, Grid* Grid, int ix, int iy, compute* EII)
{
	// Be careful, Anton's trick not in!!

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


	*EII = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy)) );

}

void Physics_computeStressInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* SII) {


	compute EII;
	compute sq_sigma_xy0,sigma_xy0, sigma_xx0, sigmaII0;

	compute khi, eta, G, dt, phi, Z;
	compute Eff_strainRate;

	compute Eps_xx, Eps_xy;

	int iCell = ix + iy*Grid->nxEC;

	Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix, iy, &EII);




	// Old stress
	//
	sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
	sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
	sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
	sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
	sigma_xx0     = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];

	sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.5*sq_sigma_xy0);

	sigma_xy0 = centerValue(Physics->sigma_xy_0, ix, iy, Grid->nxS);

	khi 		= Physics->khi[iCell];
	eta 		= Physics->eta[iCell];
	G 		    = Physics->G[iCell];
	dt 			= Physics->dt;
	phi 		= 0.0;
#if (DARCY)
	phi = Physics->phi[iCell];
#endif


	Z 	= 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
	//Eff_strainRate = sqrt(EII*EII + 1.0*Eps_xx*sigma_xx0/(G*dt) + 1.0*Eps_xy*sigma_xy0/(G*dt) + 1.0/4.0*(1.0/(G*dt))*(1.0/(G*dt))*sigmaII0*sigmaII0   );
	Eff_strainRate = EII + (1.0/(G*dt))*sigmaII0;
	*SII = (1.0-phi)*2.0*Z*Eff_strainRate;


}





void Physics_initEta(Physics* Physics, Grid* Grid, MatProps* MatProps) {

	int iy, ix, iCell;
	SinglePhase* thisPhaseInfo;
	compute P, T;
	int phase;
	compute EII, weight;
	compute B, E, V, Binc, n, taup, q, s, gamma;;
	compute invEtaDiff, invEtaDisl, invEtaPei;
	compute R = Physics->R;

	compute sumOfWeights;

	// =======================================================
	// Initial viscosity
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			//Physics->etaVisc[iCell] = Physics->eta0[iCell];
			//Physics->eta[iCell] = Physics->eta0[iCell];

			Physics->eta[iCell] = 0.0;
			Physics->G[iCell] = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];

			EII = fabs(Physics->epsRef);
#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];

			invEtaDiff = 0.0;
			invEtaDisl = 0.0;
			invEtaPei = 0.0;
			while (thisPhaseInfo != NULL) {
				weight = thisPhaseInfo->weight;
				phase = thisPhaseInfo->phase;



				if (MatProps->vDiff[phase].isActive) {
					B 			 = MatProps->vDiff[phase].B;
					E 			 = MatProps->vDiff[phase].E;
					V 			 = MatProps->vDiff[phase].V;
					Binc 		 = B*exp( - (E+V*P)/(R*T)   );
					invEtaDiff 	+= weight / (2.0*(Binc));
					//printf("B = %.2e, E = %.2e, V = %.2e, Binc = %.2e, P = %.2e, T = %.2e, - (E+V*P)/(R*T) = %.2e\n", B, E, V, Binc, P, T,- (E+V*P)/(R*T));
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					Binc 		 = B*exp( - (E+V*P)/(R*T)   );
					invEtaDisl 	 += weight / (2.0*pow(Binc,1.0/n)*pow(EII,-1.0/n+1.0));

				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					Binc 		 = B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) );
					invEtaPei 	+= weight / (2.0*pow(Binc ,1.0/s)*pow(EII,-1.0/s+1.0) );

				}


				Physics->G[iCell]	+= thisPhaseInfo->weight/MatProps->G[phase];
				thisPhaseInfo = thisPhaseInfo->next;

			}

			// Compute "isolated viscosities (i.e. viscosity as if all strain rate was caused by a single mechanism see Anton's talk)"
			if (MatProps->vDiff[phase].isActive) {
				invEtaDiff 		 = sumOfWeights	/ invEtaDiff;
			} else {
				invEtaDiff 		 = 0.0;
			}
			if (MatProps->vDisl[phase].isActive) {
				invEtaDisl 		 = sumOfWeights	/ invEtaDisl;
			} else {
				invEtaDisl 		 = 0.0;
			}
			if (MatProps->vPei[phase].isActive) {
				invEtaPei 		 = sumOfWeights	/ invEtaPei ;
			} else {
				invEtaPei	 	 = 0.0;
			}

			Physics->eta[iCell] = 1.0 / (invEtaDiff + invEtaDisl + invEtaPei);
			//printf("Physics->eta[iCell] = %.2e, invEtaDiff = %.2e, invEtaDisl = %.2e, invEtaPei = %.2e\n",Physics->eta[iCell],invEtaDiff, invEtaDisl, invEtaPei);
			Physics->G[iCell]  = Physics->sumOfWeightsCells[iCell]/Physics->G[iCell];
			Physics->khi[iCell] = 1E30;

			Physics->Z[iCell] = 1.0/( 1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );
#if (DARCY)
			Physics->eta_b[iCell] = 1E30;
			Physics->khi_b[iCell] = 1E30;
#endif

		}
	}
	Physics_copyValuesToSides(Physics->eta, Grid);
	Physics_copyValuesToSides(Physics->khi, Grid);
	Physics_copyValuesToSides(Physics->G, Grid);
	Physics_copyValuesToSides(Physics->Z, Grid);
#if (DARCY)
	Physics_copyValuesToSides(Physics->khi_b, Grid);
	Physics_copyValuesToSides(Physics->eta_b, Grid);
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


	// local copies of values
	compute eta, cohesion, frictionAngle;
	compute sigma_y, sigmaII;
	compute EII;
	compute phi = 0.0;


	compute dt = Physics->dt;

	compute Pe;
#if (DARCY)
	compute eta_b;

	compute Rad = 2.0; // radius of the griffith curve
	//compute phiMin = Numerics->phiMin;
	compute phiCrit = Numerics->phiCrit;
	compute sigmaT;//, PeSwitch;
	compute khi_b, Zb, Py;
	compute Bulk, divV, DeltaP0;
	compute DeltaP;
	//compute tol;
#endif
	//compute etaMin = Numerics->etaMin;
	//compute etaMax = Numerics->etaMax;


	//compute epsRef = Physics->epsRef;
	//compute dVxdx, dVydx, dVxdy, E_xx, E_xy;

	compute sigmaII0;
	compute Z;
	compute sigma_xx0, sq_sigma_xy0;




	compute G;

	compute khi;



	compute Eff_strainRate;
	//compute Pmin;

	//compute tol;

	//compute etaOld;


	SinglePhase* thisPhaseInfo;
	int phase;
	compute weight, sumOfWeights;

	compute B, E, V, n, taup, gamma, q, s;


	compute P, T;


	compute invEtaDiff, invEtaDisl, invEtaPei;


	compute R = Physics->R;
	compute ZUpper, ZLower;


	compute alpha = 0.25;


	compute BDiff[NB_PHASE_MAX], BDisl[NB_PHASE_MAX], BPei[NB_PHASE_MAX];

	compute maxInvVisc;

	compute tol = 1e-7;
	compute PrevZcorr, Zcorr;
	compute eta_thisPhase;

	//compute sigma_y;
#if (!DARCY)
#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi) schedule(static,16) collapse(2)
#else
//#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi, sigmaT, Bulk, khi_b, eta_b, divV, DeltaP0, Zb, DeltaP, Py) schedule(static,16) collapse(2)
#endif
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;

			//  Compute sigmaII0
			sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
			sigma_xx0  = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.5*(sq_sigma_xy0));

			Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix, iy, &EII);
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];
#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
#if (DARCY)
			phi = Physics->phi[iCell];
			//phi = 0.5;
#else
			phi = 0.0;
#endif
			//printf("MatProps->vDisl[0] = %.2e, MatProps->vDisl[1] = %.2e\n", MatProps->vDisl[0].B, MatProps->vDisl[1].B);

			alpha = 1.0;

			// Precompute B and viscosities using EII
			eta = 0.0;
			G = 0.0;
			maxInvVisc = 0.0;
			cohesion = 0.0;
			frictionAngle = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				invEtaDiff = 0.0;
				invEtaDisl = 0.0;
				invEtaPei  = 0.0;
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				G 				+= weight/MatProps->G[phase];
				cohesion 		+= MatProps->cohesion[phase] * weight;
				frictionAngle 	+= MatProps->frictionAngle[phase] * weight;
				if (MatProps->vDiff[phase].isActive) {
					B 			 = MatProps->vDiff[phase].B;
					E 			 = MatProps->vDiff[phase].E;
					V 			 = MatProps->vDiff[phase].V;
					BDiff[phase] = B*exp( - (E+V*P)/(R*T)   );
					invEtaDiff   = (2.0*(BDiff[phase]));
					maxInvVisc = fmax(invEtaDiff,maxInvVisc);
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					BDisl[phase] = B*exp( - (E+V*P)/(R*T)   );
					invEtaDisl 	 = (2.0*pow(BDisl[phase],1.0/n)*pow(EII,-1.0/n+1.0));
					maxInvVisc = fmax(invEtaDisl,maxInvVisc);
				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					BPei[phase]	 = B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) );
					invEtaPei 	 = (2.0*pow(BPei[phase] ,1.0/s)*pow(EII,-1.0/s+1.0) );
					maxInvVisc = fmax(invEtaPei,maxInvVisc);
				}
				thisPhaseInfo 	= thisPhaseInfo->next;
				eta_thisPhase = (1.0-phi)*(1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
				if (MatProps->isAir[phase] || MatProps->isWater[phase]) {
					eta_thisPhase = Numerics->StickyAirStress/(2*EII);
					eta_thisPhase = fmin(eta_thisPhase, 1e-2); // eta in the Air should not be larger than the characteristic viscosity
				}
				eta += weight * eta_thisPhase;





			}
			G 				 = sumOfWeights	/ G;
			eta 			/= sumOfWeights;
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;


			// limit eta;
			if (eta>Numerics->etaMax) {
				eta = Numerics->etaMax;
			}
			if (eta<Numerics->etaMin) {
				eta = Numerics->etaMin;
			}

			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}
			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);

			Z = 0.5*(ZUpper+ZLower);
			Zcorr = Z;

			Eff_strainRate = EII + (1.0/(G*dt))*sigmaII0;
			sigmaII = (1.0-phi)*2.0*Z*Eff_strainRate;

			// compute viscosities using sigmaII
			while (fabs(Zcorr/Z)>tol) {
				eta = 0.0;
				thisPhaseInfo = Physics->phaseListHead[iCell];


				while (thisPhaseInfo != NULL) {
					invEtaDiff = 0.0;
					invEtaDisl = 0.0;
					invEtaPei  = 0.0;
					phase = thisPhaseInfo->phase;

					if (MatProps->vDiff[phase].isActive) {
						invEtaDiff 	= (2.0*(BDiff[phase]));
					}
					if (MatProps->vDisl[phase].isActive) {
						n 			 = MatProps->vDisl[phase].n;
						invEtaDisl 	 = (2.0*BDisl[phase]*pow(sigmaII,-1.0+n));
					}
					if (MatProps->vPei[phase].isActive) {
						E 			 = MatProps->vPei[phase].E;
						V 			 = MatProps->vPei[phase].V;
						gamma 		 = MatProps->vPei[phase].gamma;
						q 			 = MatProps->vPei[phase].q;
						s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
						invEtaPei  	= ( 2.0*BPei[phase]*pow(sigmaII,-1.0+s) );
					}

					eta_thisPhase = (1.0-phi)*(1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
					if (MatProps->isAir[phase] || MatProps->isWater[phase]) {
						eta_thisPhase = Numerics->StickyAirStress/(2*EII);
						eta_thisPhase = fmin(eta_thisPhase, 1e-2); // eta in the Air should not be larger than the characteristic viscosity
						//printf("isAir, etathisPhase = %.2e\n", eta_thisPhase);
					}
					/*
					if (iy>12) {
						printf("iCell = %i, iy = %i, phase = %i, eta_thisPhase = %.2e\n", iCell, iy, phase, eta_thisPhase);
					}
					*/
					eta += weight * eta_thisPhase;
					//eta += thisPhaseInfo->weight * (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));

					thisPhaseInfo 	= thisPhaseInfo->next;
				}

				eta 			/= sumOfWeights;
				if (eta>Numerics->etaMax) {
					eta = Numerics->etaMax;
				}
				if (eta<Numerics->etaMin) {
					eta = Numerics->etaMin;
				}

				PrevZcorr = Zcorr;
				Zcorr = (1.0/(1.0/(G*dt) + 1.0/eta) - Z);
				if (Zcorr/PrevZcorr<-0.9) {
					alpha = alpha/2.0;
				}
				Z += alpha*Zcorr;

				sigmaII = (1.0-phi)*2.0*Z*Eff_strainRate;
			}

			//printf("eta = %.2e, etaMax = %.2e \n",eta, Numerics->etaMax);


			// Compute the effective Pressure Pe
#if (DARCY)
			// Limit the effective pressure
			if (phi>=phiCrit) {
				Bulk = G/sqrt(phi);
				khi_b = 1E30;
				eta_b = eta/phi;
				//if (eta_b)

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];

				Zb 	= 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));

				DeltaP = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc
				Pe =  (1.0-phi) * DeltaP;
				//Pe = Physics->Pc[iCell];
			} else {
				Pe 		= Physics->P [iCell];
				khi_b = 1E30;
				eta_b = eta/phi;
			}

			/*
			if (eta_b<1e-2) {
				eta_b = 1e-2;
			}
			if (eta_b>1e5) {
				eta_b = 1e5;
			}
			*/

			/*
			if (Physics->phase[iCell] == Physics->phaseAir || Physics->phase[iCell] == Physics->phaseWater ) {
				eta_b = 1e30;
			}
			*/


#else
			Pe 		= Physics->P [iCell];
#endif

			// compute the yield stress
			sigma_y = cohesion * cos(frictionAngle)   +   Pe * sin(frictionAngle);


#if (DARCY)
			sigmaT = (cohesion*cos(frictionAngle))/Rad;
			//tol = 1e-8;

			/*
			if (Pe < 0.0) {
				//printf("Pe<0.0\n");
				sigma_y = +Pe+(sigmaT-tol); // Pe will be shifted to 0 (arbitrary)
				//printf("A iCell = %i, Pe = %.2e, sigma_y = %.2e\n", iCell, Pe, sigma_y);
			}
			*/
			/*
			if (Pe<0) {
				sigma_y = cohesion * cos(frictionAngle);
			}
			*/

			/*
			if (Pe<-sigmaT) {
				//printf("Pe<-sigmaT\n");
				Pe = -sigmaT;
				sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle
			}
			*/

			if (Pe<0) {
				sigma_y = cohesion * cos(frictionAngle);
			}


#else
			// Since there is no griffiths handling for negative pressure for the non darcy case yet
			// here I assume a flat Mohr Coulomb when Pe <0
			if (Pe<0) {
				sigma_y = cohesion * cos(frictionAngle);
			}
#endif



			// Compute khi
			// ====================================
			khi = 1e30;
			if (sigmaII > sigma_y) {
				khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );

				if (khi<0.0) {
					// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
					// if it happens then I consider the case where there are == , which means khi -> inf
					printf("khi = %.2e, eta = %.2e, G = %.2e, dt = %.2e, Eff_Strainrate = %.2e, 1-phi = %.2e, sigma_y = %.2e, Pe = %.2e, Pmin = %.2e\n", khi, eta, G, dt, Eff_strainRate, 1.0-phi, sigma_y, Pe, -cohesion*cos(frictionAngle)/sin(frictionAngle));
					printf("WTF!\n");
					khi = 1e30;
					//exit(0);
				}


				Z 	= 1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
				sigmaII = (1.0-phi)*2.0*Z*Eff_strainRate;
			}





			// Copy updated values back
			//printf("iCell = %i, eta = %.2e, Z = %.2e, khi = %.2e, G = %.2e\n",iCell, eta, Z, khi, G);
			Physics->eta[iCell] = eta;
			Physics->khi[iCell] = khi;
			Physics->G	[iCell] = G;
			Physics->Z	[iCell] = Z;


#if (DARCY)
			Py = sigmaII - sigmaT;
			//printf("iCell = %i, khi_b before = %.2e\n", iCell, khi_b);
			khi_b = 1e30;
			if (phi>=phiCrit) {
				if (Pe < Py) {
					/*
					if (Pe/Py<0) {
						printf("icell = %i, Pe = %.2e, Py = %.2e, sigmaII = %.2e, Pc = %.2e\n", iCell, Pe, Py, sigmaII, Physics->Pc[iCell]);
						printf("WARNING: Pe and Py have opposite sense. Increasing (or decreasing?) the time step might solve the problem. iCell = %i, Pe = %.2e, Py = %.2e\n", iCell, Pe, Py);

						//exit(0);//
						//Py = 0.0;
					}
					*/
					khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
					//printf("khi_b = %.2e, phi = %.2e, Py = %.2e, (- divV + DeltaP0/(B*dt))  = %.2e, - 1.0/(Bulk*dt) = %.2e, - 1.0/eta_b = %.2e\n",khi_b, phi, Py, (- divV + DeltaP0/(B*dt)), - 1.0/(Bulk*dt), - 1.0/eta_b  );
					//printf("khi_b = %.2e, phi = %.2e, Py = %.2e\n",khi_b, phi, Py);
					if (khi_b<0.0) {
						/*
						// quite rare case where (1.0-phi)/Py * (- divV + DeltaP0/(B*dt)) <  - 1.0/(B*dt) - 1.0/eta_b
						// if it happens then I consider the case where there are == , which means khi -> inf
						printf("iCell = %i, phase = %i, khi_b = %.2e, eta_b = %.2e, Bulk = %.2e, G = %.2e, phi = %.2e\n", iCell, Physics->phase[iCell], khi_b, eta_b, Bulk, Physics->G[iCell], phi);
						printf("WTF!\n");

						//eta_b *= 0.5;
						//dt = 0.1*Physics->dt;
						//phi = 0.5*phi;
						khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
						printf("correction: khi_b = %.2e, eta_b = %.2e, (1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))= %.2e,  - 1.0/(Bulk*dt) - 1.0/eta_b = %.2e\n", khi_b, eta_b, (1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt)),  - 1.0/(Bulk*dt) - 1.0/eta_b);
						printf("Py = %.2e, sigmaT = %.2e, cohesion = %.2e, cohesion*cos(a) = %.2e, DeltaP0 = %.2e\n",Py, sigmaT, cohesion, cohesion * cos(frictionAngle), DeltaP0);
						//khi_b = 1e30;
						//exit(0);
						 */

					}



					Zb 	= 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
					Pe = (1.0-phi) * Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc
					Physics->Pc[iCell] 	  = Pe;
					/*
					if (khi_b<0.0) {
						printf("Pe = %.2e\n",Pe);
					}
					*/

				}

			}

			//printf("iCell = %i, khi_b after  = %.2e\n", khi_b);

			//printf("Pc = %.2e, Pe = %.2e\n", Physics->Pc[iCell], Pe);
			Physics->eta_b[iCell] = eta_b;
			Physics->khi_b[iCell] = khi_b;
#endif








			//printf("%.2e  ", khi);
		}
		//printf("\n");
	}


	//Physics_copyValuesToSides(Physics->etaVisc, Grid);
	Physics_copyValuesToSides(Physics->eta, Grid);
	Physics_copyValuesToSides(Physics->khi, Grid);
	Physics_copyValuesToSides(Physics->G, Grid);
	Physics_copyValuesToSides(Physics->Z, Grid);
#if (DARCY)
	Physics_copyValuesToSides(Physics->khi_b, Grid);
	Physics_copyValuesToSides(Physics->eta_b, Grid);
#endif





	// ================================================================================
	// 									Shear nodes viscosity
	int iNode;
#pragma omp parallel for private(iy,ix, iNode) schedule(static,32)
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


















void Physics_updateDt(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics)
{

	Physics->dtDarcy = 1e100;
	printf("In: Physics->dt = %.2e\n", Physics->dt);
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


	Physics->dtAdv 	= Numerics->CFL_fac_Stokes*Numerics->dLmin/(Physics->maxV); // note: the min(dx,dy) is the char length, so = 1
	compute Kappa = 0.0;
	compute minKappa = 1e10;
	int i;
	for (i = 0; i < MatProps->nPhase; ++i) {
		if (MatProps->isAir[i] == false && MatProps->isWater[i] == false) {
			Kappa = MatProps->k[i] / (MatProps->rho0[i]*Physics->Cp);
			if (Kappa<minKappa) {
				minKappa = Kappa;
			}
		}
	}

	Physics->dtT 	= Numerics->CFL_fac_Thermal*fmin(Grid->dx, Grid->dy)*fmin(Grid->dx, Grid->dy)/(minKappa);
	//printf("WTF   ===  minKappa = %.2e, k = %.2e, rho = %.2e, Cp = %.2e, Numerics->CFL_fac_Thermal = %.2e\n",minKappa, MatProps->k[0], MatProps->rho0[0], Physics->Cp, Numerics->CFL_fac_Thermal);
	//printf("perm_eta_f = %.2e, phi = %.2e Physics->Pf[0] = %.2e\n",Physics->perm_eta_f[0],Physics->phi[0],Physics->Pf[0]);
	int iCell, iy, ix;
#if (DARCY)
	/*
	if (Numerics->timeStep==1 & Numerics->itNonLin == 0) {

	}
	 */
	Physics->dtDarcy 	= 1e100;//10.00*Numerics->CFL_fac*fmin(Grid->dx, Grid->dy)/(3*Physics->minPerm);


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

	compute saveL, saveT, saveV;

	compute minCompactionLength = 1e100;

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

			VelFluidX = DarcyVelX/phi ;// + VelSolidX;
			VelFluidY = DarcyVelY/phi ;// + VelSolidY;

			VelFluid = sqrt(VelFluidX*VelFluidX + VelFluidY*VelFluidY);

			CompactionTime = CompactionLength/VelFluid;
			//CFLtime =Numerics->dLmin/VelFluid;


			if (CompactionLength<minCompactionLength) {
				minCompactionLength = CompactionLength;
			}



			if (CompactionTime*Numerics->CFL_fac_Darcy<Physics->dtDarcy ) {

				Physics->dtDarcy = CompactionTime*Numerics->CFL_fac_Darcy;
				//saveV = VelFluid;
				//saveT = CompactionTime;
				//saveL = CompactionLength;
				//printf("CompactionLength = %.2e, DarcyVel = %.2e, Vx = %.2e, VelFluid = %.2e\n",CompactionLength, (sqrt(DarcyVelX*DarcyVelX + DarcyVelY*DarcyVelY)), Physics->Vx[10], VelFluid);
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
		//printf("phase = %i, CompactionTime = %.2e, CompactionLength =%.2e, VelFluid = %.2e, VelSolid = %.2e, DarcyVelY/phi = %.2e, DarcyVelY = %.2e\n", Physics->phase[iCell], CompactionTime, CompactionLength, VelFluid, sqrt(VelSolidX*VelSolidX + VelSolidY*VelSolidY), DarcyVelY/phi, DarcyVelY);
	}

	//printf("C.L = %.2e, C.time = %.2e, FluidVel = %.2e\n",saveL, saveT, saveV);
	printf("dx = %.2e, dy = %.2e, minCompactionLength = %.2e\n",Grid->dx, Grid->dy, minCompactionLength);


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
	compute khi;
	dtElastic = 1E10;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;

			Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix,iy, &EII);
			eta = Physics->eta[iCell];
			khi = Physics->khi[iCell];
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
				this_dt = fabs(  eta / ( ( (2.0*eta*EII - sigmaII0)/DSigmaMax -1.0 )*Physics->G[iCell] )    ); // Absolute Dsigma

				if (this_dt<dtElastic) {
					dtElastic = this_dt;
					//compute estDsigma = (2.0*eta*EII - sigmaII0) * (Physics->G[iCell]*this_dt/(eta + Physics->G[iCell]*this_dt));
					//printf("DsigmaII = %.2e, Dsigma_xx = %.2e, dtElastic = %.2e, sigmaII0 = %.2e, sigma_xx0 = %.2e, estDsigma = %.2e\n",DsigmaII, Physics->Dsigma_xx_0[iCell], dtElastic, sigmaII0, Physics->sigma_xx_0[iCell],estDsigma);
				}
			}

			if (Physics->phase[iCell]!=Physics->phaseAir && Physics->phase[iCell]!=Physics->phaseWater) {
				dtMaxwell = (1.0/(1.0/eta+1.0/khi))/Physics->G[iCell];
				if (dtMaxwell > Physics->dtMaxwellMax) {
					Physics->dtMaxwellMax = dtMaxwell;
					//printf("iy = %i, dtMaxwell = %.2e, eta = %.2e, khi = %.2e, G = %2.e \n",iy, dtMaxwell, eta, khi, Physics->G[iCell]);
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
	/*
	if (Numerics->timeStep<=0){// && Numerics->itNonLin==0) {
		Physics->dt = 1E-1;
	} else if (Numerics->timeStep==1){// && Numerics->itNonLin==0) {
		//Physics->dt = Physics->dt;
	} else {
	*/
	//if (Numerics->timeStep>0){

		corr = (Physics->dt-dtOld);
		if (Numerics->timeStep>0){
			if (corr>dtOld) {
				corr = dtOld;
			} else if (corr< -dtOld) {
				corr = -dtOld;
			}
		}
		Physics->dt = dtOld + 0.2*corr;

		// Relimit
		//Physics->dt  =  fmin(Physics->dt,1.0*Physics->dtAdv);
		//Physics->dt  =  fmin(Physics->dt,1.0*Physics->dtDarcy);

		//printf("dtOld = %.2e, dt = %.2e, corr = %.2e, 0.1*corr = %.2e, dtOld+0.1*corr = %.2e\n", dtOld, Physics->dt, corr, 0.1*corr, dtOld + 0.1*corr);
		//printf("A - Physics->dt = %.2e, dtAdv = %.2e, dtT = %.2e, PdtDarcy = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, dtMin = %.2e, dtMax = %.2e\n", Physics->dt, Physics->dtAdv, Physics->dtT, Physics->dtDarcy, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Numerics->dtMin, Numerics->dtMax);

	//}

	//printf("A - Physics->dt = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, Physics->dtAdv = %.2e, Physics->dtT = %.2e, Physics->dtDarcy = %.2e\n", Physics->dt, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Physics->dtAdv, Physics->dtT, Physics->dtDarcy);



	//Physics->dtAdv = Physics->dt;
	//Physics->dtT = Physics->dt;

#if (DARCY)
	//Physics->dtDarcy = Physics->dt;
#endif

	//Numerics->dtMax  = 1e-2;
	//Numerics->dtMin = Physics->dtMaxwellMin + 0.2*(Physics->dtMaxwellMax - Physics->dtMaxwellMin);
	//Numerics->dtMax = Physics->dtMaxwellMax - 0.2*(Physics->dtMaxwellMax - Physics->dtMaxwellMin);

	if (Numerics->use_dtMaxwellLimit) {
		if (MatProps->G[Physics->phaseRef] < 1E30) { // to enable switching off the elasticity
			compute dtMinOld, dtMaxOld;
			dtMinOld = Numerics->dtMin;
			dtMaxOld = Numerics->dtMax;
			Numerics->dtMin = pow(10,log10(Physics->dtMaxwellMin) + 0.1*(log10(Physics->dtMaxwellMax) - log10(Physics->dtMaxwellMin) ));
			Numerics->dtMax = pow(10,log10(Physics->dtMaxwellMax) - 0.0*(log10(Physics->dtMaxwellMax) - log10(Physics->dtMaxwellMin) ));
			if (Numerics->timeStep>0 || Numerics->itNonLin>1){
				corr = Numerics->dtMin - dtMinOld;
				Numerics->dtMin = dtMinOld + 0.1*corr;
				corr = Numerics->dtMax - dtMaxOld;
				Numerics->dtMax = dtMaxOld + 0.1*corr;
			}
		}
	}

	/*
	if (Physics->dtAdv<Numerics->dtMin) {
		Physics->dtAdv = Numerics->dtMin;
	} else if (Physics->dtAdv>Numerics->dtMax) {
		Physics->dtAdv = Numerics->dtMax;
	}
	*/

	if (Physics->dt<Numerics->dtMin) {
		Physics->dt = Numerics->dtMin;
	} else if (Physics->dt>Numerics->dtMax) {
		Physics->dt = Numerics->dtMax;
	}
	printf("B - Physics->dt = %.2e, dtAdv = %.2e, dtT = %.2e, PdtDarcy = %.2e, dtMaxwellMin = %.2e, dtMaxwellMax = %.2e, dtMin = %.2e, dtMax = %.2e\n", Physics->dt, Physics->dtAdv, Physics->dtT, Physics->dtDarcy, Physics->dtMaxwellMin ,Physics->dtMaxwellMax, Numerics->dtMin, Numerics->dtMax);
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







}


#if (DARCY)
void Physics_computePerm(Physics* Physics, Grid* Grid, Numerics* Numerics, MatProps* MatProps)
{
	Physics->minPerm = 1E100;
	int iy, ix;
	int iCell;
	compute phi;
	compute phiRef = 0.0001;
	compute PermEffRef = MatProps->perm0_eta_f[0]  *  phiRef*phiRef*phiRef  / ( (1.0-phiRef)*(1.0-phiRef));

	compute perm0;
	SinglePhase* thisPhaseInfo;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		phi = Physics->phi[iCell];

		//if (Physics->phase[iCell] != Physics->phaseAir && Physics->phase[iCell] != Physics->phaseWater) {

			perm0 = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				perm0 += MatProps->perm0_eta_f[thisPhaseInfo->phase] * thisPhaseInfo->weight;
				thisPhaseInfo = thisPhaseInfo->next;
			}
			perm0 /= Physics->sumOfWeightsCells[iCell];

			Physics->perm_eta_f[iCell] = perm0  *  phi*phi*phi  * ( (1.0-phi)*(1.0-phi));

		//} else {
		//	Physics->perm_eta_f[iCell]=1e6*PermEffRef;
		//}

	}



	Physics_copyValuesToSides(Physics->perm_eta_f, Grid);

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


void Physics_computePhi(Physics* Physics, Grid* Grid, Numerics* Numerics)
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
			//printf("iCell = %i, phi0 = %.2e, phi = %.2e, dt =%.2e, divV0=%.2e, divV = %.2e, dt*divV = %.2e\n", iCell, Physics->phi0[iCell], Physics->phi[iCell], dt, Physics->divV0[iCell], divV, dt*divV);

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




			//if (iCell == 150) {
			//printf("Physics->Vx[ix+iy*nxVx] = %.2e, divV = %.2e, phi0 = %.2e, phi = %.2e, dt = %.2e\n", Physics->Vx[ix+iy*nxVx], divV, Physics->phi0[iCell], Physics->phi[iCell], dt);
			//}

			sum += Physics->phi[iCell];
		}
		//printf("divV = %.2e\n",divV);
	}
	//printf("                    sum phi = %.2e\n", sum);


	Physics_copyValuesToSides(Physics->phi, Grid);
	Physics_copyValuesToSides(Physics->Dphi, Grid);


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













#endif


void Physics_getValuesToSidesFromBC(compute* ECValues, Grid* Grid, BC* BC, Numbering* Numbering) {
	// Replace boundary values by their neighbours
	int INeigh, iy, ix, I;
	int IBC;
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
		IBC = abs(Numbering->map[I])-1; // BC nodes are numbered -1 to -n
		ECValues[I] = Physics_computeSideValuesFromBC_ForOneCell(ECValues[INeigh], BC, IBC, ix, iy, Grid);
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
		IBC = abs(Numbering->map[I])-1; // BC nodes are numbered -1 to -n
		ECValues[I] = Physics_computeSideValuesFromBC_ForOneCell(ECValues[INeigh], BC, IBC, ix, iy, Grid);
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
			IBC = abs(Numbering->map[I])-1; // BC nodes are numbered -1 to -n
			ECValues[I] = Physics_computeSideValuesFromBC_ForOneCell(ECValues[INeigh], BC, IBC, ix, iy, Grid);

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			IBC = abs(Numbering->map[I])-1; // BC nodes are numbered -1 to -n
			ECValues[I] = Physics_computeSideValuesFromBC_ForOneCell(ECValues[INeigh], BC, IBC, ix, iy, Grid);

		}
	}
}

compute Physics_computeSideValuesFromBC_ForOneCell(compute neighValue, BC* BC, int IBC, int ix, int iy, Grid* Grid)
{
	compute sideValue = -1.0;
	// BCtype: BC->type[IBC]
	// sideValue: Val[iCell]
	// neighValue: EqSystem->x[INeigh]*scale
	// BCValue: BC->value[IBC]
	if (BC->type[IBC]==DirichletGhost) { // Dirichlet
		sideValue = 2.0*BC->value[IBC] - neighValue;
		//	printf("IBC %i is Dir Ghost\n",IBC);
	}
	else if (BC->type[IBC]==NeumannGhost) { // Neumann
		if (ix==0)  {// left or bottom boundary
			sideValue = neighValue - BC->value[IBC]*Grid->DXEC[0];
		} else if (ix==Grid->nxEC-1) {
			sideValue = neighValue + BC->value[IBC]*Grid->DXEC[Grid->nxEC-2];
		}
		if (iy==0) { // right or top boundary
			sideValue = neighValue - BC->value[IBC]*Grid->DYEC[0];
		} else if (iy==Grid->nyEC-1) { // right or top boundary
			sideValue = neighValue + BC->value[IBC]*Grid->DYEC[Grid->nyEC-2];
		}
	}
	else if (BC->type[IBC]==Dirichlet) {
		sideValue = BC->value[IBC];
	}
	else if (BC->type[IBC]==Infinity) {
		if (ix==0)  {// left or bottom boundary
			sideValue = neighValue * BC->DeltaL/(BC->DeltaL+Grid->DXEC[0]) + BC->value[IBC] * Grid->DXEC[0]/(BC->DeltaL+Grid->DXEC[0]);
		} else if (ix==Grid->nxEC-1) {
			sideValue = neighValue * BC->DeltaL/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]) + BC->value[IBC] * Grid->DXEC[Grid->nxEC-2]/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]);
		}
		if (iy==0) { // right or top boundary
			sideValue = neighValue * BC->DeltaL/(BC->DeltaL+Grid->DYEC[0]) + BC->value[IBC] * Grid->DYEC[0]/(BC->DeltaL+Grid->DYEC[0]);
		} else if (iy==Grid->nyEC-1) { // right or top boundary
			sideValue = neighValue * BC->DeltaL/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]) + BC->value[IBC] * Grid->DYEC[Grid->nyEC-2]/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]);
		}
	}
	else {
		sideValue = 0.0;
		printf("error in Physics_get_ECVal_FromSolution: unknown boundary type\n");
		exit(0);
	}
	return sideValue;


}


void Physics_copyValuesToSides(compute* ECValues, Grid* Grid)
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

void Physics_copyValuesToSidesi(int* ECValues, Grid* Grid)
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









void Physics_computeRho(Physics* Physics, Grid* Grid, MatProps* MatProps)
{

	int iCell;
	SinglePhase* thisPhaseInfo;
#pragma omp parallel for private(iCell, thisPhaseInfo) schedule(dynamic,16)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->rho_g[iCell] = 0.0;
		thisPhaseInfo = Physics->phaseListHead[iCell];
		while (thisPhaseInfo != NULL) {
			Physics->rho_g[iCell] += MatProps->rho0_g[thisPhaseInfo->phase] * thisPhaseInfo->weight;
			thisPhaseInfo = thisPhaseInfo->next;
		}
		Physics->rho_g[iCell] /= Physics->sumOfWeightsCells[iCell];

#if (DARCY)

		Physics->rho_g[iCell] = (1.0 - Physics->phi[iCell])*Physics->rho_g[iCell] + Physics->phi[iCell]*Physics->rho_f_g;

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

	Physics_copyValuesToSides(Physics->rho_g, Grid);


}



compute Physics_getFromMatProps_ForOneCell(Physics* Physics, compute* ListFromMatProps, MatProps* MatProps, int iCell) {
	SinglePhase* thisPhaseInfo;
	compute value = 0.0;

	thisPhaseInfo = Physics->phaseListHead[iCell];
	while (thisPhaseInfo != NULL) {
		value += ListFromMatProps[thisPhaseInfo->phase] * thisPhaseInfo->weight;
		thisPhaseInfo = thisPhaseInfo->next;
	}
	return value /= Physics->sumOfWeightsCells[iCell];
}





void Physics_get_ECVal_FromSolution (compute* Val, int ISub, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
	// Where Val is the value to extract from the solution, and DVal the increment since the last time step, IStep is the index of the subsystem of equations
	int I, IBC, INeigh, iy, ix;
	int INumMap0 = Numbering->subEqSystem0Dir[ISub];
	//printf("eq0 = %i, ISub = %i\n", INumMap0, ISub);
	int iCell;


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

				Val[iCell] = Physics_computeSideValuesFromBC_ForOneCell(EqSystem->x[INeigh]*scale, BC, IBC, ix, iy, Grid);
				/*
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
				else if (BC->type[IBC]==Infinity) {
					if (ix==0)  {// left or bottom boundary
						Val[iCell] = EqSystem->x[INeigh]*scale * BC->DeltaL/(BC->DeltaL+Grid->DXEC[0]) + BC->value[IBC] * Grid->DXEC[0]/(BC->DeltaL+Grid->DXEC[0]);
					} else if (ix==Grid->nxEC-1) {
						Val[iCell] = EqSystem->x[INeigh]*scale * BC->DeltaL/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]) + BC->value[IBC] * Grid->DXEC[Grid->nxEC-2]/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]);
					}
					if (iy==0) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh]*scale * BC->DeltaL/(BC->DeltaL+Grid->DYEC[0]) + BC->value[IBC] * Grid->DYEC[0]/(BC->DeltaL+Grid->DYEC[0]);
					} else if (iy==Grid->nyEC-1) { // right or top boundary
						Val[iCell] = EqSystem->x[INeigh]*scale * BC->DeltaL/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]) + BC->value[IBC] * Grid->DYEC[Grid->nyEC-2]/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]);
					}
				}
				else {
					printf("error in Physics_get_ECVal_FromSolution: unknown boundary type\n");
					exit(0);
				}
				 */
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

	int phaseAir = Physics->phaseAir;
	int phaseWater;
	if (Physics->phaseWater==-1) {
		phaseWater = Physics->phaseAir;
	} else {
		phaseWater = Physics->phaseWater;
	}
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

			if (contribPhase[phaseAir]>0) {
				Physics->phase[iCell] = phaseAir;
			} else if (contribPhase[phaseWater]>0) {
				Physics->phase[iCell] = phaseWater;
			} else {
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
	}


	Physics_copyValuesToSidesi(Physics->phase,Grid);

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





void Physics_reinitPhaseList(Physics* Physics, Grid* Grid) {
	int iCell;
	SinglePhase* temp;
#pragma omp parallel for private(iCell, temp) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		while (Physics->phaseListHead[iCell]->next!=NULL) {
			temp = Physics->phaseListHead[iCell];
			Physics->phaseListHead[iCell] = Physics->phaseListHead[iCell]->next;
			free(temp);
		}
		Physics->phaseListHead[iCell]->phase = -1;
		Physics->phaseListHead[iCell]->weight = 0.0;
	}



}


void Physics_check(Physics* Physics, Grid* Grid, Char* Char) {

	printf("=== Physics_check ===\n");
	int iCell, ix, iy;
	compute* Data;
	int iData;
	int nData = 10;
#if (HEAT)
	nData +=1;
#endif
#if (DARCY)
	nData +=6;
#endif

	compute s 	= Char->time;			// second
	compute m 	= Char->length; 		// meter
	compute kg 	= Char->mass; 			// kilogram
	compute K 	= Char->temperature; 	// Kelvin

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s); 		// Watt

	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	compute mol = 1.0;


	compute norm_g = sqrt(Physics->g[0]*Physics->g[0] + Physics->g[1]*Physics->g[1]);

	bool Dim = true;
	compute unit = 1.0;

	for (iData = 0; iData < nData; ++iData) {
		switch (iData) {
		case 0:
			printf("=====    G    =====\n");
			Data = Physics->G;
			if (Dim) unit = Pa;
			break;
		case 1:
			printf("=====   eta   =====\n");
			Data = Physics->eta;
			if (Dim) unit = Pas;
			break;
		case 2:
			printf("=====   khi   =====\n");
			Data = Physics->khi;
			if (Dim) unit = Pas;
			break;
		case 3:
			printf("=====    Z    =====\n");
			Data = Physics->Z;
			if (Dim) unit = Pas;
			break;
		case 4:
			printf("=====  rho_g  =====\n");
			Data = Physics->rho_g;
			if (Dim) unit = kg/m/m/m   *m/s/s ;
			break;
		case 5:
			printf("=====  rho  =====\n");
			Data = Physics->rho_g;
			if (Dim) unit =  1.0/norm_g * kg/m/m/m ;
			break;
		case 6:
			printf("=====  sigma_xx_0  =====\n");
			Data = Physics->sigma_xx_0;
			if (Dim) unit = Pa;
			break;
		case 7:
			printf("=====  Dsigma_xx_0  =====\n");
			Data = Physics->Dsigma_xx_0;
			if (Dim) unit = Pa;
			break;
		case 8:
			printf("=====  sumOfWeightsCells  =====\n");
			Data = Physics->sumOfWeightsCells;
			if (Dim) unit = 1.0;
			break;
		case 9:
			printf("=====  	 P    =====\n");
			Data = Physics->P;
			if (Dim) unit = Pa;
			break;
		case 10:
#if (HEAT)
			printf("=====    T    =====\n");
			Data = Physics->T;
			if (Dim) unit = K;
#endif
			break;
		case 11:
#if (DARCY)
			printf("=====   phi   =====\n");
			Data = Physics->phi;
			if (Dim) unit = 1.0;
#endif
			break;
		case 12:
#if (DARCY)
			printf("=====    Pc    =====\n");
			Data = Physics->Pc;
			if (Dim) unit = Pa;
#endif
			break;
		case 13:
#if (DARCY)
			printf("=====    Pf    =====\n");
			Data = Physics->Pf;
			if (Dim) unit = Pa;
#endif
			break;
		case 14:
#if (DARCY)
			printf("=====    khi_b    =====\n");
			Data = Physics->khi_b;
			if (Dim) unit = Pas;
#endif
			break;
		case 15:
#if (DARCY)
			printf("=====    eta_b    =====\n");
			Data = Physics->eta_b;
			if (Dim) unit = Pas;
#endif
			break;
		case 16:
#if (DARCY)
			printf("=====    perm    =====\n");
			Data = Physics->perm_eta_f;
			if (Dim) unit = Physics->eta_f * m*m    ;
#endif
			break;
		}
		printf("Char unit = %.2e\n",unit);
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				iCell = ix+iy*Grid->nxEC;
				printf("%.2e  ", Data[iCell]*unit);
			}
			printf("\n");
		}
	}




}

