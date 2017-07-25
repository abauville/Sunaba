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

#if (CRANK_NICHOLSON_VEL || INERTIA)
	Physics->Vx0 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy0 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
#endif
#if (CRANK_NICHOLSON_VEL && CRANK_NICHOLSON_P)
	Physics->P0				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
#endif

	Physics->P 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->Z 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->eta 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->khi 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

	Physics->rho 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

#if (STRAIN_SOFTENING)
	Physics->strain 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->Dstrain 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
#endif

#if (HEAT)

	Physics->k 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->T0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->DT 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
#endif


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
	Physics->khi_b 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // bulk plasticity
	Physics->Zb				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) ); // bulk effective viscosity

#endif



	Physics->G 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->khiShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->etaShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->ZShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->phase 			= (int*) 	malloc( Grid->nECTot * sizeof(int) );


	// Initialize stuff
	//int i;
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0.0;
#if (CRANK_NICHOLSON_VEL || INERTIA)
		Physics->Vx0[i] = 0.0;
#endif
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0.0;
#if (CRANK_NICHOLSON_VEL || INERTIA)
		Physics->Vy0[i] = 0.0;
#endif
	}

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nECTot; ++i) {

		Physics->khi[i] = 0;
#if (STRAIN_SOFTENING)
		Physics->strain[i] = 0;
		Physics->Dstrain[i] = 0;
#endif

#if (HEAT)
		Physics->T[i]  = 1.0;
		Physics->DT[i] = 0.0;
#endif

		Physics->P[i] = 0.0;
#if (CRANK_NICHOLSON_VEL && CRANK_NICHOLSON_P)
		Physics->P0[i] = 0.0;
#endif

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


#if (CRANK_NICHOLSON_VEL)
	free(Physics->Vx0);
	free(Physics->Vy0);
#if (CRANK_NICHOLSON_P)
	free(Physics->P0);
#endif
#endif

	free(Physics->P );

	free( Physics->eta );


	free(Physics->etaShear);
	free( Physics->khi );
	free( Physics->khiShear );
	free( Physics->ZShear );

	free( Physics->rho );

#if (STRAIN_SOFTENING)
	free(Physics->strain);
	free(Physics->Dstrain);
#endif


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
	free(Physics->Zb);
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
					rho_g_h = Physics->rho[iCell] * -Physics->g[1] * (-0.5*Grid->DYEC[iy-1] );
				} else {
					rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellN]) * -Physics->g[1] * Grid->DYEC[iy] ;
				}
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
						rho_g_h = Physics->rho[iCell] * Physics->g[0] * (-0.5*Grid->DXEC[ix] );
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellW]) * Physics->g[0] * Grid->DXEC[ix-1] ;
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
						rho_g_h = Physics->rho[iCell] * -Physics->g[0] * (-0.5*Grid->DXEC[ix-1] );
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellE]) * -Physics->g[0] * Grid->DXEC[ix] ;
					}
					Physics->P[iCell] += rho_g_h;
				}
			}
		}
	}



	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {

#if (DARCY)
		Physics->Pf[iCell] = Physics->P[iCell];
		Physics->Pc[iCell] = 0.0;
		Physics->DeltaP0[iCell] = 0.0;
		Physics->DDeltaP[iCell] = 0.0;
#endif
	}




}











void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumStokes, Numbering* NumThermal, BC* BCThermal)
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
		Physics->sigma_xx_0 [iCell] = 0.0;
		Physics->sumOfWeightsCells [iCell] = 0.0;

#if (HEAT)

		Physics->T[iCell] = 0.0;
#endif
#if (DARCY)
		Physics->DeltaP0[iCell] 		= 0.0;
		Physics->phi0[iCell] 		= 0.0;
#endif
#if (STRAIN_SOFTENING)
		Physics->strain[iCell] 		= 0.0;
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
	int ix,  iy, C;

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

	for (iColor = 0; iColor < 4; ++iColor) {
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight, thisPhaseInfo) schedule(static,16)
		for (iy = iyStart[iColor]; iy < Grid->nyS; iy+=2) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ix+=2) { // I don't get why though
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				// Loop through the particles in the shifted cell
				// ======================================
				while (thisParticle!=NULL) {

					locX = 2.0*(thisParticle->x-Grid->X[ix]);
					locY = 2.0*(thisParticle->y-Grid->Y[iy]);

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
						weight = fabs((locX + xMod[i])   *   (locY + yMod[i]));


						// Get the phase and weight of phase contribution for each cell
						thisPhaseInfo = Physics->phaseListHead[iCell];

						while (thisPhaseInfo->phase != phase) {
							if (thisPhaseInfo->next == NULL) {

								if (!changedHead[iCell]) {
									thisPhaseInfo->phase = phase;
									changedHead[iCell] = true;
									break;
								} else {

									addSinglePhase(&Physics->phaseListHead[iCell],phase);
									thisPhaseInfo = Physics->phaseListHead[iCell];
									break;
								}


							} else {
								thisPhaseInfo = thisPhaseInfo->next;
							}
						}
						thisPhaseInfo->weight += weight;


						// For properties that are stored on the markers, sum contributions
						Physics->sigma_xx_0		[iCell] += thisParticle->sigma_xx_0 * weight;
#if (HEAT)
						Physics->T				[iCell] += thisParticle->T * weight;
#endif
#if (DARCY)
						Physics->DeltaP0		[iCell] += thisParticle->DeltaP0 * weight;
						Physics->phi0			[iCell] += thisParticle->phi * weight;
#endif
#if (STRAIN_SOFTENING)
						Physics->strain			[iCell] += thisParticle->strain * weight;
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
#if (STRAIN_SOFTENING)
				Physics->strain			[iCellD] += Physics->strain			[iCellS];
				Physics->strain			[iCellS]  = Physics->strain			[iCellD];
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
#if (STRAIN_SOFTENING)
		Physics->strain		[iCell] /= Physics->sumOfWeightsCells	[iCell];
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
#if (STRAIN_SOFTENING)
	Physics_copyValuesToSides(Physics->strain, Grid);
#endif




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
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, signX, signY, phase, i, iNodeNeigh, weight) schedule(static,16)
		for (iy = iyStart[iColor]; iy < Grid->nyS; iy+=2) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ix+=2) { // I don't get why though
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





	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->DT[i] = (Physics->T[i] - Physics->T0[i]) * Physics->dtAdv/Physics->dtT;
	}

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

				dtDiff = (Physics->Cp*rhoParticle)/(  MatProps->k[phase]*( 2.0/(Grid->dx*Grid->dx) + 2.0/(Grid->dy*Grid->dy) )  );


				DT_sub_OnThisPart = ( TFromNodes - thisParticle->T ) * ( 1.0 - exp(-d * Physics->dtAdv/dtDiff) );

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
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		DT_sub_OnThisCell = ( DT_sub_OnTheCells[I+0] + DT_sub_OnTheCells[I+1] + DT_sub_OnTheCells[I+2] + DT_sub_OnTheCells[I+3]) / sum ; // harmonic average
		DT_rem_OnTheCells[iCell] = Physics->DT[iCell] - DT_sub_OnThisCell;

	}

	Physics_copyValuesToSides(DT_rem_OnTheCells, Grid);



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

#if (DARCY)
		int i;
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DDeltaP[i] *=  Physics->dtAdv/Physics->dt;
			Physics->Dphi[i] *=  Physics->dtAdv/Physics->dt;
			//printf("Physics->dtAdv/Physics->d = %.2e\n",Physics->dtAdv/Physics->dt);
		}
#endif

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

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

				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {

					thisParticle->DeltaP0 = 0.0;


				} else {

				thisParticle->DeltaP0 += ( .25*(1.0-locX)*(1.0-locY)*Physics->DDeltaP[ix  +(iy  )*Grid->nxEC]
																					  + .25*(1.0-locX)*(1.0+locY)*Physics->DDeltaP[ix  +(iy+1)*Grid->nxEC]
																																   + .25*(1.0+locX)*(1.0+locY)*Physics->DDeltaP[ix+1+(iy+1)*Grid->nxEC]
																																												+ .25*(1.0+locX)*(1.0-locY)*Physics->DDeltaP[ix+1+(iy  )*Grid->nxEC] );
				}
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



void Physics_interpStrainFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

#if (STRAIN_SOFTENING)
	int iCell;
	compute SII;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			Physics_computeStressInvariantForOneCell(Physics, Grid, ix, iy, &SII);
			Physics->Dstrain[iCell] = SII/(2.0*Physics->khi[iCell])*Physics->dtAdv; // Recovering the incremental plastic strain
		}
	}
	Physics_copyValuesToSides(Physics->Dstrain, Grid);
#endif

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

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

#if (STRAIN_SOFTENING)


				thisParticle->strain += ( .25*(1.0-locX)*(1.0-locY)*Physics->Dstrain[ix  +(iy  )*Grid->nxEC]
																					 + .25*(1.0-locX)*(1.0+locY)*Physics->Dstrain[ix  +(iy+1)*Grid->nxEC]
																																  + .25*(1.0+locX)*(1.0+locY)*Physics->Dstrain[ix+1+(iy+1)*Grid->nxEC]
																																											   + .25*(1.0+locX)*(1.0-locY)*Physics->Dstrain[ix+1+(iy  )*Grid->nxEC] );

#endif

				thisParticle = thisParticle->next;
			}
		}
	}


}




void Physics_interpStressesFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal, MatProps* MatProps, Numerics* Numerics)
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

	compute d_ve_ini = 0.00;
	compute dtm = Physics->dtAdv;
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
#pragma omp parallel for private(iy, ix, i, iNode, thisParticle, locX, locY, signX, signY, sigma_xx_0_fromNodes, sigma_xy_0_fromNodes, eta, khi, G, dtMaxwell, Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart, iNodeNeigh, weight, iCell) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================

			compute d_ve = d_ve_ini;


			while (thisParticle!=NULL) {

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


				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {

					// Compute Dsigma sub grid
					Dsigma_xx_sub_OnThisPart =   0.0;
					Dsigma_xy_sub_OnThisPart =   0.0;

					// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
					thisParticle->sigma_xx_0 += Dsigma_xx_sub_OnThisPart;
					thisParticle->sigma_xy_0 += Dsigma_xy_sub_OnThisPart;

				}
				else {

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
				eta_vp = fmax(eta_vp,Numerics->etaMin);
				//printf("eta_vp = %.2e\n",eta_vp);
				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from

				locX = fabs(locX)-1.0;
				locY = fabs(locY)-1.0;


				sigma_xy_0_fromNodes = ( .25*(1.0-locX)*(1.0-locY)*Physics->sigma_xy_0[ix      +(iy  )    *Grid->nxS]
																					   + .25*(1.0-locX)*(1.0+locY)*Physics->sigma_xy_0[ix      +(iy+signY)*Grid->nxS]
																																	   + .25*(1.0+locX)*(1.0+locY)*Physics->sigma_xy_0[ix+signX+(iy+signY)*Grid->nxS]
																																													   + .25*(1.0+locX)*(1.0-locY)*Physics->sigma_xy_0[ix+signX+(iy  )    *Grid->nxS] );





				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				G = MatProps->G[thisParticle->phase];

				dtMaxwell = eta_vp/G;
				dtMaxwell = fmin(dtm,dtMaxwell);

				// Compute Dsigma sub grid
				Dsigma_xx_sub_OnThisPart =   ( sigma_xx_0_fromNodes - thisParticle->sigma_xx_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
				Dsigma_xy_sub_OnThisPart = ( sigma_xy_0_fromNodes - thisParticle->sigma_xy_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
				//if (( 1.0 - exp(-d_ve * dtm/dtMaxwell))<0.8) {
				//printf("( 1.0 - exp(-d_ve * dtm/dtMaxwell) = %.2e\n", ( 1.0 - exp(-d_ve * dtm/dtMaxwell)));
				//}

				// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
				thisParticle->sigma_xx_0 += Dsigma_xx_sub_OnThisPart;
				thisParticle->sigma_xy_0 += Dsigma_xy_sub_OnThisPart;
				}


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

					weight = (locX + xModNode[i]*0.5)   *   (locY + yModNode[i]*0.5);

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

				thisParticle = thisParticle->next;

			}
		}
	}


	int I;
	compute sum;
	compute Dsigma_xy_sub_OnThisNode;
#pragma omp parallel for private(iNode, I, sum, Dsigma_xy_sub_OnThisNode) schedule(static,32)
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		I = 4*iNode;
		sum = sumOfWeights_OnTheNodes[I+0] + sumOfWeights_OnTheNodes[I+1] + sumOfWeights_OnTheNodes[I+2] + sumOfWeights_OnTheNodes[I+3];
		if (sum==0) {
			printf("error in Physics_interpStressesFromCellsToParticles: node #%i received no contribution from particles\n", iNode );
			exit(0);
		}

		Dsigma_xy_sub_OnThisNode = ( Dsigma_xy_sub_OnTheNodes[I+0] +  Dsigma_xy_sub_OnTheNodes[I+1] +  Dsigma_xy_sub_OnTheNodes[I+2] +  Dsigma_xy_sub_OnTheNodes[I+3]) / sum ; // harmonic average

		Dsigma_xy_rem_OnTheNodes[iNode] = Physics->Dsigma_xy_0[iNode] - Dsigma_xy_sub_OnThisNode;
	}



	compute Dsigma_xx_sub_OnThisCell;
#pragma omp parallel for private(iCell, I, sum, Dsigma_xx_sub_OnThisCell) schedule(static,32)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		I = 4*iCell;
		sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		Dsigma_xx_sub_OnThisCell = ( Dsigma_xx_sub_OnTheCells[I+0] + Dsigma_xx_sub_OnTheCells[I+1] + Dsigma_xx_sub_OnTheCells[I+2] + Dsigma_xx_sub_OnTheCells[I+3]) / sum ; // harmonic average
		Dsigma_xx_rem_OnTheCells[iCell] = Physics->Dsigma_xx_0[iCell] - Dsigma_xx_sub_OnThisCell;
	}

	// Copy values to sides
	Physics_copyValuesToSides(Dsigma_xx_rem_OnTheCells, Grid);







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

				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->sigma_xx_0 = 0.0;
				} else {
					thisParticle->sigma_xx_0  += ( .25*(1.0-locX)*(1.0-locY)*Dsigma_xx_rem_OnTheCells[ix  +(iy  )*Grid->nxEC]
																								  + .25*(1.0-locX)*(1.0+locY)*Dsigma_xx_rem_OnTheCells[ix  +(iy+1)*Grid->nxEC]
																																					   + .25*(1.0+locX)*(1.0+locY)*Dsigma_xx_rem_OnTheCells[ix+1+(iy+1)*Grid->nxEC]
																																																			+ .25*(1.0+locX)*(1.0-locY)*Dsigma_xx_rem_OnTheCells[ix+1+(iy  )*Grid->nxEC] );
				}


				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from
				if (locX<0.0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0.0) {
					signY = -1;
				} else {
					signY = 1;
				}


				locX = fabs(locX)-1.0;
				locY = fabs(locY)-1.0;

				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->sigma_xy_0 = 0.0;
				} else {
					if (Grid->isFixed && ix<1) { // should be better optimized, and also include the right boundary, for the case where injection happens
						thisParticle->sigma_xy_0 = 0.0;
					} else {
						thisParticle->sigma_xy_0  += ( .25*(1.0-locX)*(1.0-locY)*Dsigma_xy_rem_OnTheNodes[ix      +(iy  )    *Grid->nxS]
																								  + .25*(1.0-locX)*(1.0+locY)*Dsigma_xy_rem_OnTheNodes[ix      +(iy+signY)*Grid->nxS]
																																					   + .25*(1.0+locX)*(1.0+locY)*Dsigma_xy_rem_OnTheNodes[ix+signX+(iy+signY)*Grid->nxS]
																																																			+ .25*(1.0+locX)*(1.0-locY)*Dsigma_xy_rem_OnTheNodes[ix+signX+(iy  )    *Grid->nxS] );
					}

				}

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





void Physics_eulerianAdvectVel(Grid* Grid, Physics* Physics, BC* BCStokes, Numbering* NumStokes)
{
#if (INERTIA || CRANK_NICHOLSON_VEL)
	int ix, iy;
	compute dVxdx, dVxdy, dVydx, dVydy;
	compute dVxdx0, dVxdy0, dVydx0, dVydy0;
	compute* VxNew = (compute*) malloc(Grid->nVxTot * sizeof(compute));
	compute* VyNew = (compute*) malloc(Grid->nVyTot * sizeof(compute));
	compute Vx, Vy;
	compute dt = Physics->dt;
	for (iy = 1; iy < Grid->nyVx-1; ++iy) {
		for (ix = 1; ix < Grid->nxVx-1; ++ix) {
			dVxdx = (Physics->Vx[ix+1 +  iy   *Grid->nxVx] - Physics->Vx[ix-1 +  iy   *Grid->nxVx])/(2.0*Grid->dx);
			dVxdy = (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   + (iy-1)*Grid->nxVx])/(2.0*Grid->dy);
			dVxdx0 = (Physics->Vx0[ix+1 +  iy   *Grid->nxVx] - Physics->Vx0[ix-1 +  iy   *Grid->nxVx])/(2.0*Grid->dx);
			dVxdy0 = (Physics->Vx0[ix   + (iy+1)*Grid->nxVx] - Physics->Vx0[ix   + (iy-1)*Grid->nxVx])/(2.0*Grid->dy);
			Vy = 0.25* (Physics->Vy[ix   + (iy  )*Grid->nxVy] + Physics->Vy[ix+1 + (iy  )*Grid->nxVy] + Physics->Vy[ix   + (iy-1)*Grid->nxVy] + Physics->Vy[ix+1 + (iy-1)*Grid->nxVy]);
			//VxNew[ix+iy*Grid->nxVx] = Physics->Vx[ix   +  iy   *Grid->nxVx]*(1.0-dt*dVxdx) - dt*Vy*dVxdy;
			VxNew[ix+iy*Grid->nxVx] = Physics->Vx[ix   +  iy   *Grid->nxVx]*(1.0-dt*.5*(dVxdx+dVxdx0)) - dt*Vy*.5*(dVxdy+dVxdy0);
		}
	}

	for (iy = 1; iy < Grid->nyVy-1; ++iy) {
		for (ix = 1; ix < Grid->nxVy-1; ++ix) {
			dVydx = (Physics->Vy[ix+1 +  iy   *Grid->nxVy] - Physics->Vy[ix-1 +  iy   *Grid->nxVy])/(2.0*Grid->dx);
			dVydy = (Physics->Vy[ix   + (iy+1)*Grid->nxVy] - Physics->Vy[ix   + (iy-1)*Grid->nxVy])/(2.0*Grid->dy);
			dVydx0 = (Physics->Vy0[ix+1 +  iy   *Grid->nxVy] - Physics->Vy0[ix-1 +  iy   *Grid->nxVy])/(2.0*Grid->dx);
			dVydy0 = (Physics->Vy0[ix   + (iy+1)*Grid->nxVy] - Physics->Vy0[ix   + (iy-1)*Grid->nxVy])/(2.0*Grid->dy);
			Vx = 0.25* (Physics->Vx[ix   + (iy  )*Grid->nxVx] + Physics->Vx[ix-1 + (iy  )*Grid->nxVx] + Physics->Vx[ix   + (iy+1)*Grid->nxVx] + Physics->Vx[ix-1 + (iy+1)*Grid->nxVx]);
			//VyNew[ix+iy*Grid->nxVy] =  Physics->Vy[ix   +  iy   *Grid->nxVy]*(1.0-dt*dVydy) - Vx*dt*dVydx;
			VyNew[ix+iy*Grid->nxVy] =  Physics->Vy[ix   +  iy   *Grid->nxVy]*(1.0-dt*.5*(dVydy+dVydy0)) - Vx*dt*.5*(dVydx+dVydx0);
		}
	}




	int iVx, iVy, InoDir;
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			iVx = ix + iy*Grid->nxVx;
			InoDir = NumStokes->map[iVx];
			if (Grid->isPeriodic) {
				printf("error:  in Physics_interpFromParticlestoCell: the implementation of the interpolation of velocities from particles to cell is not finished for the case of periodic BC");
			}
			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vx [iVx] = VxNew[iVx];
				Physics->Vx0[iVx] = VxNew[iVx];
			} else {
				Physics->Vx0[iVx] = Physics->Vx[iVx];
			}
		}
	}

	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			iVy = ix + iy*Grid->nxVy;
			InoDir = NumStokes->map[iVy + Grid->nVxTot];
			if (Grid->isPeriodic) {
				printf("error:  in Physics_interpFromParticlestoCell: the implementation of the interpolation of velocities from particles to cell is not finished for the case of periodic BC");
			}
			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vy [iVy] = VyNew[iVy];
				Physics->Vy0[iVy] = VyNew[iVy];
			} else {
				Physics->Vy0[iVy] = Physics->Vy[iVy];
			}
		}
	}


	free(VxNew);
	free(VyNew);
	#endif
}



void Physics_get_VxVy_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem, Numerics* Numerics)
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

			scale = 1.0;//EqSystem->S[InoDir];

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
				else if (BC->type[IBC]==Neumann) { // Neumann on normal node
					// Get neighbours index

					if (ix==0) { // left boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVx  ];
						if (INeigh<0) {
							if (iy==0) {
								INeigh = Numbering->map[  ix+1 + (iy+1)*Grid->nxVx  ];
							} else if (iy==Grid->nyVx-1) {
								INeigh = Numbering->map[  ix+1 + (iy-1)*Grid->nxVx  ];
							}
						}
						Physics->Vx[I] = EqSystem->x[INeigh]*scale;// - BC->value[IBC] *Grid->dx/(2*Physics->Z[ix+1 + (iy)*Grid->nxEC ]);
					} else if (ix==Grid->nxVx-1) { // right boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVx  ];
						if (INeigh<0) {
							if (iy==0) {
								INeigh = Numbering->map[  ix-1 + (iy+1)*Grid->nxVx  ];
							} else if (iy==Grid->nyVx-1) {
								INeigh = Numbering->map[  ix-1 + (iy-1)*Grid->nxVx  ];
							}
						}
						Physics->Vx[I] = EqSystem->x[INeigh]*scale;// + BC->value[IBC] *Grid->dx/(2*Physics->Z[ix + (iy)*Grid->nxEC ]);
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet. (Neumann Vx)\n");
						exit(0);
					}
				}
				else { // on a ghost node

					// Get neighbours index
					if (iy==0) { // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVx  ];
					} else if (iy==Grid->nyVx-1) { // lower boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVx  ];
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet. (Ghost Vx)\n");
						exit(0);
					}

					scale = 1.0;//EqSystem->S[INeigh];

					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vx[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh]*scale;
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (iy==0)  // lower boundary
							Physics->Vx[I] = EqSystem->x[INeigh]*scale - BC->value[IBC]/Physics->ZShear[ix + 0*Grid->nxS]*Grid->dy;
						if (iy==Grid->nyVx-1)  // top boundary
							Physics->Vx[I] = EqSystem->x[INeigh]*scale + BC->value[IBC]/Physics->ZShear[ix + (Grid->nyS-1)*Grid->nxS]*Grid->dy;
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
				else if (BC->type[IBC]==Neumann) {
					// Get neighbours index
					if (iy==0) { // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVy + Grid->nVxTot ];
						if (INeigh<0) {
							if (ix==0) {
								INeigh = Numbering->map[  ix+1 + (iy+1)*Grid->nxVy  ];
							} else if (ix==Grid->nxVy-1) {
								INeigh = Numbering->map[  ix-1 + (iy+1)*Grid->nxVy  ];
							}
						}
						Physics->Vy[I] = EqSystem->x[INeigh]*scale;// - BC->value[IBC] *Grid->dx/(2*Physics->Z[ix + (iy+1)*Grid->nxEC ]);
					} else if (iy==Grid->nyVy-1) { // top boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVy + Grid->nVxTot ];
						if (INeigh<0) {
							if (ix==0) {
								INeigh = Numbering->map[  ix+1 + (iy-1)*Grid->nxVy  ];
							} else if (ix==Grid->nxVy-1) {
								INeigh = Numbering->map[  ix-1 + (iy-1)*Grid->nxVy  ];
							}
						}
						Physics->Vy[I] = EqSystem->x[INeigh]*scale;// + BC->value[IBC] *Grid->dx/(2*Physics->Z[ix + (iy  )*Grid->nxEC ]);
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet. (Neumann Vy)\n");
						exit(0);
					}

				}
				else { // on a ghost node

					// Get neighbours index
					if (ix==0) {  // left boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVy + Grid->nVxTot ];
					} else if (ix==Grid->nxVy-1) { // right boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVy + Grid->nVxTot  ];
					} else {
						INeigh = 0;
						printf("error internal BC are not properly taken into account yet. (Ghost Vy)\n");
						exit(0);
					}
					scale = 1.0;//EqSystem->S[INeigh];

					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vy[I] = 2.0*BC->value[IBC] - EqSystem->x[INeigh]*scale;
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (ix==0)  // left boundary
							Physics->Vy[I] = EqSystem->x[INeigh]*scale - BC->value[IBC]/Physics->ZShear[0 + iy*Grid->nxS]*Grid->dx;
						if (ix==Grid->nxVy-1)  // right boundary
							Physics->Vy[I] = EqSystem->x[INeigh]*scale + BC->value[IBC]/Physics->ZShear[Grid->nxS-1 + iy*Grid->nxS]*Grid->dx;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}
		}
	}



#if (CRANK_NICHOLSON_VEL)
	compute weight[2];
	if (Numerics->timeStep>0) {
		weight[0] =  0.5;
		weight[1] =  0.5;
	} else {
		weight[0] =  1.0;
		weight[1] =  0.0;
	}

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = weight[0]*Physics->Vx[i] + weight[1]*Physics->Vx0[i];
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = weight[0]*Physics->Vy[i] + weight[1]*Physics->Vy0[i];
	}

#endif



	compute maxVx, maxVy;
	compute Vx, Vy;
	maxVx = 0.0;
	maxVy = 0.0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			Vx = (Physics->Vx[ix-1+  iy   *Grid->nxVx]+Physics->Vx[ix+  iy   *Grid->nxVx])/2.0;
			Vy = (Physics->Vy[ix  + (iy-1)*Grid->nxVy]+Physics->Vx[ix+ (iy-1)*Grid->nxVy])/2.0;
			maxVx = fmax(maxVx, Vx);
			maxVy = fmax(maxVy, Vy);
		}
	}
	Physics->maxVx = maxVx;
	Physics->maxVy = maxVy;



}

#if (CRANK_NICHOLSON_VEL || INERTIA)
void Physics_updateOldVel_P				(Physics* Physics, Grid* Grid)
{

	// A better method would be to intervert the pointers;
	int i;

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx0[i] = Physics->Vx[i];

	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy0[i] = Physics->Vy[i];
	}

#if (CRANK_NICHOLSON_P)
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P0[i] = Physics->P[i];
	}
#endif





}
#endif


void Physics_get_P_FromSolution(Physics* Physics, Grid* Grid, BC* BCStokes, Numbering* NumStokes, EqSystem* EqStokes, Numerics* Numerics)
{
	int ix, iCell;

#if (!DARCY)


	// /!\ For visu it's better if all sides are Neumann
	Physics_get_ECVal_FromSolution (Physics->P, 2, Grid, BCStokes, NumStokes, EqStokes);

	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	compute RefPressure = Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];// - 1.0;//Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];//Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure;
	}

#if (CRANK_NICHOLSON_VEL)
	#if (CRANK_NICHOLSON_P)
	compute weight[2];
	int i;
	if (Numerics->timeStep>0) {
		weight[0] =  0.5;
		weight[1] =  0.5;
	} else {
		weight[0] =  1.0;
		weight[1] =  0.0;
	}

#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P[i] = weight[0]*Physics->P[i] + weight[1]*Physics->P0[i];
	}

	#endif
#endif


#else

	int i;
	Physics_get_ECVal_FromSolution (Physics->Pf, 2, Grid, BCStokes, NumStokes, EqStokes);
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



	// Fill P, the total pressure
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = Physics->Pc[iCell] + Physics->Pf[iCell];
	}



#endif



}


#if (HEAT)
void Physics_get_T_FromSolution(Physics* Physics, Grid* Grid, BC* BCThermal, Numbering* NumThermal, EqSystem* EqThermal, Numerics* Numerics)
{
	Physics_get_ECVal_FromSolution (Physics->T, 0, Grid, BCThermal, NumThermal, EqThermal);
}
#endif





void Physics_computeStressChanges(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes, Numerics* Numerics)
{

	// see Taras' book p. 186
	int ix, iy, iCell, iNode;
	compute Z;
	compute Eps_xx, Eps_xy;
	compute dVxdy, dVydx, dVxdx, dVydy;
	compute G;
	compute phi;


	//compute phi;
	// compute stress
	//#pragma omp parallel for private(iy, ix, iCell, Eps_xx, Z) schedule(static,32)
	compute dt = Physics->dt;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;

#if (DARCY)
			phi = Physics->phi[iCell];
#else
			phi = 0.0;
#endif

			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
								 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
								 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;


			Eps_xx = 0.5*(dVxdx-dVydy);


			compute Ds0_old = Physics->Dsigma_xx_0[iCell];


			Physics->Dsigma_xx_0[iCell] = Physics->Z[iCell]/(1.0-phi)*(2.0*Eps_xx + Physics->sigma_xx_0[iCell]/(Physics->G[iCell]*dt)) - Physics->sigma_xx_0[iCell];

			Physics->Dsigma_xx_0[iCell] *= Physics->dtAdv/Physics->dt; // To update by the right amount according to the time step


			if (Numerics->timeStep>0) {
				//Physics->Dsigma_xx_0[iCell] = 1.0/2.0*Physics->Dsigma_xx_0[iCell] + 1.0/2.0*Ds0_old; // Crank-Nicolson
				//Physics->Dsigma_xx_0[iCell] = .7*Physics->Dsigma_xx_0[iCell] + .3*Ds0_old; // empirical
				//Physics->Dsigma_xx_0[iCell] = 1.0/sqrt(2.0)*Physics->Dsigma_xx_0[iCell] + (1.0-1.0/sqrt(2.0))*Ds0_old; // empirical
			}

			//Physics->Dsigma_xx_0[iCell] = 0.0;
		}
	}
	//printf("Physics->sigma_xx_0[2+10*Grid->nxC] = %.2e, Physics->Dsigma_xx_0[2+10*Grid->nxC] = %.2e  ", Physics->sigma_xx_0[2+10*Grid->nxC], Physics->Dsigma_xx_0[2+10*Grid->nxC]);
	printf("dtAdv = %.2e, dt = %.2e\n", Physics->dtAdv, Physics->dt);


	// Replace boundary values by their neighbours
	// lower boundary
	Physics_copyValuesToSides(Physics->Dsigma_xx_0, Grid);




	//#pragma omp parallel for private(iy, ix, iNode, dVxdy, dVydx, Eps_xy, GShear, etaShear, Z) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;

#if (DARCY)
			phi = shearValue(Physics->phi,  ix   , iy, Grid->nxEC);
#else
			phi = 0.0;
#endif

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]
								  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]
								  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Eps_xy = 0.5*(dVxdy+dVydx);


			G 	 	= shearValue(Physics->G, ix, iy, Grid->nxEC);

			Z 	 	= Physics->ZShear[iNode];

			compute Ds0_old = Physics->Dsigma_xy_0[iNode];

			Physics->Dsigma_xy_0[iNode] = Z/(1.0-phi) * (2.0*Eps_xy + Physics->sigma_xy_0[iNode]/(G*dt)) - Physics->sigma_xy_0[iNode];


			Physics->Dsigma_xy_0[iNode] *= Physics->dtAdv/Physics->dt;


			if (Numerics->timeStep>0) {
				//Physics->Dsigma_xy_0[iNode] = 1.0/2.0*Physics->Dsigma_xy_0[iNode] + 1.0/2.0* Ds0_old; // empirical
				//Physics->Dsigma_xy_0[iNode] = .7*Physics->Dsigma_xy_0[iNode] + .3* Ds0_old; // empirical
				//Physics->Dsigma_xy_0[iNode] = 1.0/sqrt(2.0)*Physics->Dsigma_xy_0[iNode] + (1.0-1.0/sqrt(2.0))* Ds0_old;

			}

			// Ensure free slip
			if (ix==0 && BC->IsFreeSlipLeft) {
				Physics->Dsigma_xy_0[iNode] = 0.0;
			}
			if (ix==Grid->nxS && BC->IsFreeSlipRight) {
				Physics->Dsigma_xy_0[iNode] = 0.0;
			}
			if (iy == 0 && BC->IsFreeSlipBot) {
				Physics->Dsigma_xy_0[iNode] = 0.0;
			}
			if (iy==Grid->nyS && BC->IsFreeSlipTop) {
				Physics->Dsigma_xy_0[iNode] = 0.0;
			}
		}
	}



#if (DARCY)
	compute Bulk, Zb, divV, DeltaP0, DeltaP;
	//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;

			phi = Physics->phi[iCell];
			Bulk = Physics->G[iCell]/sqrt(phi);

			divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
			divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
			DeltaP0 = Physics->DeltaP0[iCell];


			Zb 	= Physics->Zb[iCell];

			DeltaP = Zb/(1.0-phi) * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

			Physics->DDeltaP[iCell] = DeltaP - Physics->DeltaP0[iCell];
			Physics->DDeltaP[iCell] *= Physics->dtAdv/Physics->dt;
		}
	}
	Physics_copyValuesToSides(Physics->DDeltaP, Grid);
#endif


}









void Physics_computeStrainRateInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* EII)
{
	compute dVxdy, dVydx, dVxdx, dVydy;

	compute ShearComp_sqr;
	int iNode, Ix, Iy;
	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};
	dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
			 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
	dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

	// Method A: using the averaging of derivatives on the four nodes
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
	*EII = sqrt(  (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy))  +  0.25*ShearComp_sqr );


}


void Physics_computeStrainRateInvariantForOneNode(Physics* Physics, BC* BCStokes, Grid* Grid, int ix, int iy, compute* EII)
{
	// Be careful, Anton's trick not in!!

	compute dVxdy, dVydx, dVxdx, dVydy;

	dVxdy = (Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx]
						 - Physics->Vx[(ix  ) + (iy  )*Grid->nxVx])/Grid->dy;

	dVydx = (Physics->Vy[(ix+1) + (iy  )*Grid->nxVy]
						 - Physics->Vy[(ix  ) + (iy  )*Grid->nxVy])/Grid->dx;


	compute dVxdxCell[4], dVydyCell[4]; // order: NE, NW, SW, SE

	// use Anton's trick for the inner nodes
	if (ix>0 && ix<Grid->nxS-1 && iy>0 && iy<Grid->nyS-1) {
		dVxdxCell[0] = Physics->Vx[(ix+1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy+1)*Grid->nxVx];
		dVxdxCell[1] = Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx];
		dVxdxCell[2] = Physics->Vx[(ix  )+(iy  )*Grid->nxVx] - Physics->Vx[(ix-1)+(iy  )*Grid->nxVx];
		dVxdxCell[3] = Physics->Vx[(ix+1)+(iy  )*Grid->nxVx] - Physics->Vx[(ix  )+(iy  )*Grid->nxVx];

		dVydyCell[0] = Physics->Vy[(ix+1)+(iy+1)*Grid->nxVy] - Physics->Vy[(ix+1)+(iy  )*Grid->nxVy];
		dVydyCell[1] = Physics->Vy[(ix  )+(iy+1)*Grid->nxVy] - Physics->Vy[(ix  )+(iy  )*Grid->nxVy];
		dVydyCell[2] = Physics->Vy[(ix  )+(iy  )*Grid->nxVy] - Physics->Vy[(ix  )+(iy-1)*Grid->nxVy];
		dVydyCell[3] = Physics->Vx[(ix+1)+(iy  )*Grid->nxVx] - Physics->Vx[(ix+1)+(iy-1)*Grid->nxVx];
		compute NormalComp_sqr = 0.0;
		int iCell;
		for (iCell = 0; iCell < 4; ++iCell) {

			dVxdx = dVxdxCell[iCell];
			dVydy = dVydyCell[iCell];
			NormalComp_sqr += (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy)) ;

		}
		*EII = sqrt( 0.25*NormalComp_sqr +  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))   );

	} else {
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
		*EII = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  (0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy)) );

	}


}

void Physics_computeStressInvariantForOneCell(Physics* Physics, Grid* Grid, int ix, int iy, compute* SII) {




	int iCell = ix + iy*Grid->nxEC;



	int Method = 0; // 0 compute from strain invariant, 1 compute from Dsigma



	//sigma_xy0 = centerValue(Physics->sigma_xy_0, ix, iy, Grid->nxS);
	if (Method == 0) {
		compute EII;
		compute sq_sigma_xy0,sigma_xy0, sigma_xx0, sigmaII0;

		compute khi, eta, G, dt, phi, Z;
		compute Eff_strainRate;


		Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix, iy, &EII);

		// Old stress
		sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
		sigma_xx0     = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];

		sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*sq_sigma_xy0);

		compute dVxdy, dVydx, dVxdx, dVydy, Eps_xy, Eps_xx;
		dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
					 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

		dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
					 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

		Eps_xx = 0.5*(dVxdx-dVydy);


		//EII = sqrt(Eps_xx*Eps_xx + Eps_xy*Eps_xy);
		// Anton's trick
		dVxdy = 0;
		dVydx = 0;
		compute Exy_x_Sxy0 = 0.0;
		int iNode, Ix, Iy;
		int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
		int IyMod[4] = {0,0,1,1};
		for (iNode = 0; iNode < 4; ++iNode) {
			Ix = (ix-1)+IxMod[iNode];
			Iy = (iy-1)+IyMod[iNode];

			dVxdy = ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
								  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


			dVydx = ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
								  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;

			Exy_x_Sxy0 += (0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0[Ix+Iy*Grid->nxS];
		}
		Exy_x_Sxy0 /= 4.0; // Eps_xy*sigma_xy0



		khi 		= Physics->khi[iCell];
		eta 		= Physics->eta[iCell];
		G 		    = Physics->G[iCell];
		dt 			= Physics->dt;
		phi 		= 0.0;
#if (DARCY)
		phi = Physics->phi[iCell];
#endif


		Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(2.0*G*dt));
		Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(2.0*G*dt) + Exy_x_Sxy0/(2.0*G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
		*SII = 2.0*Z*Eff_strainRate;
	} else if (Method == 1) {
		compute sq_sigma_xy,sigma_xy, sigma_xx, sigmaII;
		sq_sigma_xy  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];

		sq_sigma_xy += Physics->Dsigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->Dsigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->Dsigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->Dsigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->Dsigma_xy_0[ix  +(iy  )*Grid->nxS];

		sigma_xx     = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];

		sigmaII = sqrt((sigma_xx)*(sigma_xx)    + 0.25*sq_sigma_xy);
	}

}





void Physics_initEta(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics) {

	int iy, ix, iCell;
	SinglePhase* thisPhaseInfo;
	compute P, T;
	int phase;
	compute EII, weight;
	compute B, E, V, Binc, n, taup, q, s, gamma;
	compute invEtaDiff, invEtaDisl, invEtaPei;
	compute R = Physics->R;
	compute eta, G, cohesion, frictionAngle, eta_thisPhase;
	compute sumOfWeights;

	// =======================================================
	// Initial viscosity
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			//Physics->etaVisc[iCell] = Physics->eta0[iCell];
			//Physics->eta[iCell] = Physics->eta0[iCell];

			thisPhaseInfo = Physics->phaseListHead[iCell];

			EII = fabs(Physics->epsRef)/1.0;
#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];

			eta = 0.0;
			G = 0.0;
			cohesion = 0.0;
			frictionAngle = 0.0;
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
					invEtaDiff   = (2.0*(B*exp( - (E+V*P)/(R*T)   )));
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					invEtaDisl 	 = (2.0*pow(B*exp( - (E+V*P)/(R*T)   ),1.0/n)*pow(EII,-1.0/n+1.0));
				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					invEtaPei 	 = (2.0*pow(B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) ) ,1.0/s)*pow(EII,-1.0/s+1.0) );
				}
				thisPhaseInfo 	= thisPhaseInfo->next;
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));

				eta += weight * eta_thisPhase;





			}
			eta = eta / sumOfWeights;
			if (eta>Numerics->etaMax) {
				eta = Numerics->etaMax;
			}
			if (eta<Numerics->etaMin) {
				eta = Numerics->etaMin;
			}

			Physics->eta[iCell] = eta;

			Physics->G[iCell]  = Physics->sumOfWeightsCells[iCell]/G;
			Physics->khi[iCell] = 1E30;

			Physics->Z[iCell] = 1.0/( 1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );

			if (Physics->Z[iCell]<Numerics->etaMin) {
				Physics->Z[iCell] = Numerics->etaMin;
			}
#if (DARCY)
			Physics->eta_b[iCell] = Physics->eta[iCell]/(Physics->phi[iCell]);
			Physics->khi_b[iCell] = 1e30;
			Physics->Zb[iCell] 	  = 1.0/( 1.0/Physics->eta_b[iCell] + 1.0/(Physics->G[iCell]/(sqrt(Physics->phi[iCell]))*Physics->dt) );
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
	Physics_copyValuesToSides(Physics->Zb, Grid);
#endif


	int iNode;
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = shearValue(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = shearValue(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = shearValue(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}

}




void Physics_computeEta(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes,MatProps* MatProps)
{
	int iCell, iy, ix;

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

	compute sigmaII0;
	compute Z;
	compute sigma_xx0, sq_sigma_xy0;




	compute G;

	compute khi;



	compute Eff_strainRate;

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

	compute* ZprePlasticity = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* sigma_y_Stored = (compute*) malloc(Grid->nECTot * sizeof(compute));

#if (STRAIN_SOFTENING)
	compute strainReductionFac = 0.9; // 1.0 stays the same
#endif


	//compute sigma_y;
#if (!DARCY)
#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi) schedule(static,16) collapse(2)
#else
#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi, sigmaT, Bulk, khi_b, eta_b, divV, DeltaP0, Zb, DeltaP, Py) schedule(static,16) collapse(2)
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
			sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

			Physics_computeStrainRateInvariantForOneCell(Physics, Grid, ix, iy, &EII);
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];



			compute dVxdy, dVydx, dVxdx, dVydy, Eps_xx;
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Eps_xx = 0.5*(dVxdx-dVydy);


			// Anton's trick
			dVxdy = 0;
			dVydx = 0;
			compute Exy_x_Sxy0 = 0.0;
			int iNode, Ix, Iy;
			int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
			int IyMod[4] = {0,0,1,1};
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];

				dVxdy = ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
									  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


				dVydx = ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
									  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;

				Exy_x_Sxy0 += (0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0[Ix+Iy*Grid->nxS];
			}
			Exy_x_Sxy0 /= 4.0; // Eps_xy*sigma_xy0














#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
#if (DARCY)
			phi = Physics->phi[iCell];
#else
			phi = 0.0;
#endif
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
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));


				eta += weight * eta_thisPhase;





			}
			G 				 = sumOfWeights	/ G;
			eta 			/= sumOfWeights;
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;


#if (STRAIN_SOFTENING)
			compute strainLimit = 1.0;
			compute coeff = (1.0-Physics->strain[iCell]/strainLimit);
			if (coeff<strainReductionFac){
				coeff = strainReductionFac;
			}
			frictionAngle *= coeff;
#endif



			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}
			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);

			Z = 0.5*((1.0-phi)*ZUpper+(1.0-phi)*ZLower);
			Zcorr = Z;

			Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(2.0*G*dt) + Exy_x_Sxy0/(2.0*G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
			sigmaII = 2.0*Z*Eff_strainRate;

			// compute viscosities using sigmaII
			while (fabs(Zcorr/Z)>tol) {
				eta = 0.0;
				thisPhaseInfo = Physics->phaseListHead[iCell];

				while (thisPhaseInfo != NULL) {
					invEtaDiff = 0.0;
					invEtaDisl = 0.0;
					invEtaPei  = 0.0;
					phase = thisPhaseInfo->phase;
					weight = thisPhaseInfo->weight;
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

					eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
					eta += weight * eta_thisPhase;
					thisPhaseInfo 	= thisPhaseInfo->next;



				}

				eta 			/= sumOfWeights;

				PrevZcorr = Zcorr;
				Zcorr = (1.0-phi)*(1.0/(1.0/(G*dt) + 1.0/eta)) - Z;
				if (Zcorr/PrevZcorr<-0.9) {
					alpha = alpha/2.0;
				}
				Z += alpha*Zcorr;

				sigmaII = 2.0*Z*Eff_strainRate;
			}


			// Compute the effective Pressure Pe
#if (DARCY)
			// Limit the effective pressure
			if (phi>=phiCrit) {
				Bulk = G/sqrt(phi);
				khi_b = 1E30;
				eta_b = eta/phi;

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];

				Zb 	= (1.0-phi)*1.0/(1.0/eta_b + 1.0/(Bulk*dt));
				Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

			} else {
				Pe 		= Physics->P [iCell];
				khi_b = 1E30;
				eta_b = eta/phi;

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];

				Bulk = G/sqrt(phi);
				Zb 	= (1.0-phi)* 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
			}


#else
			compute Pf = 0.0;
			Pe 		= Physics->P [iCell] - Pf;

#endif

			Z 	= (1.0-phi)*1.0/(1.0/eta + 1.0/(G*dt));


			sigmaII = 2.0*Z*Eff_strainRate;

			ZprePlasticity[iCell] = Z/(1.0-phi);


			khi = 1e30;
#if (DARCY)
			khi_b = 1e30;
#endif
			alpha = 1.0;



			sigma_y = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);


#if (DARCY)
			sigmaT = (cohesion*cos(frictionAngle))/Rad;

			if (Pe<-sigmaT) {

				Py = -sigmaT;
				khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
				Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
				Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

				sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle
			}
#else
			// Since there is no griffiths handling for negative pressure for the non darcy case yet
			// here I assume a flat Mohr Coulomb when Pe <0
			if (Pe<0) {
				sigma_y = cohesion * cos(frictionAngle);
			}
#endif

			sigma_y_Stored[iCell] = sigma_y;

			compute Pe0;

			sigmaII0 = sigmaII;
			Pe0 = Pe;
#if (DARCY)
			Py = sigmaII - sigmaT;
			khi_b = 1e30;
#endif
			compute yieldTol = 1e-4;
			int iDum = 0;

			do {
				sigmaII = sigmaII0;

				iDum++;


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
					Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
					sigmaII = 2.0*Z*Eff_strainRate;

				} else {
					khi = 1e30;
				}



#if (DARCY)
				Py = sigmaII - sigmaT;

#if (1)
				if (phi>=phiCrit) {

					Pe0 = Pe;
					if (Pe < Py && Pe!=0) {

						compute khi_bOld = khi_b;
						khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b );
						Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
						Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc
						if (isnan(Pe)!=0) {
							printf("error in computeEta: Pe  = %.2e\n", Pe);
							exit(0);
						}

					} else {
						khi_b = 1e30;
					}

					// sigma_y chunk
					// =================
					sigma_y = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);
					sigmaT = (cohesion*cos(frictionAngle))/Rad;
					if (Pe<-sigmaT) {
						Py = -sigmaT;
						khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
						Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
						Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc


						sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle
					}

					//printf("koko\n");
					// =================
				} else {
					break;
				}
#endif


				if (fabs(1.0-Pe/Py)<yieldTol || Pe>Py-yieldTol) {
					if (sigmaII0<sigma_y) {
						break;
					}
					if (fabs(1.0-sigmaII/sigma_y)<yieldTol) {
						break;
					}
				}
#else
				break;
#endif

			} while (iDum<1) ;

			// Copy updated values back
			Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt)); // this might not be needed, but it's an extra security
			if (Z<Numerics->etaMin) {
				Z = Numerics->etaMin;
			}


			Physics->eta[iCell] = eta;
			Physics->khi[iCell] = khi;
			Physics->G	[iCell] = G;
			Physics->Z	[iCell] = Z;
#if (DARCY)
			Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt)); // this might not be needed, but it's an extra security

			if (fabs(Zb)<Numerics->etaMin) {
				Zb = Zb/fabs(Zb) * Numerics->etaMin;
			}
			if (fabs(Zb)>Numerics->etaMax) {
				Zb = Zb/fabs(Zb) * Numerics->etaMax;
			}


			Physics->eta_b[iCell] = eta_b;
			Physics->khi_b[iCell] = khi_b;
			Physics->Zb[iCell] = Zb;
#endif

		}
	}


	//Physics_copyValuesToSides(Physics->etaVisc, Grid);
	Physics_copyValuesToSides(Physics->eta, Grid);
	Physics_copyValuesToSides(Physics->khi, Grid);
	Physics_copyValuesToSides(Physics->G, Grid);
	Physics_copyValuesToSides(Physics->Z, Grid);

	Physics_copyValuesToSides(sigma_y_Stored, Grid);
	Physics_copyValuesToSides(ZprePlasticity, Grid);

#if (DARCY)
	Physics_copyValuesToSides(Physics->khi_b, Grid);
	Physics_copyValuesToSides(Physics->eta_b, Grid);
	Physics_copyValuesToSides(Physics->Zb, Grid);
#endif





	// ================================================================================
	// 									Shear nodes viscosity
	compute sq_sigma_xx0;
	compute sigma_xy0;
	int iNode;
	//#pragma omp parallel for private(iy,ix, iNode) schedule(static,32)
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = shearValue(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = shearValue(Physics->khi,  ix   , iy, Grid->nxEC);
#if (DARCY)
			phi = shearValue(Physics->phi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = shearValue(Physics->Z,  ix   , iy, Grid->nxEC);
#else
			/*
			phi = 0.0;

			eta = shearValue(Physics->eta,  ix   , iy, Grid->nxEC);
			G = shearValue(Physics->G,  ix   , iy, Grid->nxEC);

			sigma_y = shearValue(sigma_y_Stored,  ix   , iy, Grid->nxEC);
			sq_sigma_xx0  = Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC];
			sigma_xy0  	  = Physics->sigma_xy_0[iNode];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xy0)*(sigma_xy0)    + 0.25*(sq_sigma_xx0));


			Physics_computeStrainRateInvariantForOneNode(Physics,BCStokes,Grid,ix,iy,&EII);


			Eff_strainRate = EII + (1.0/(2.0*G*dt))*sigmaII0;

			Z = shearValue(ZprePlasticity,  ix   , iy, Grid->nxEC);
			Z = Z*(1.0-phi);

			//printf("Z = %.2e, Zoth = %.2e, eta = %.2e, etaGrid = %.2e, G = %.2e, GGrid = %.2e\n",Z, 1.0/(1.0/eta + 1/(G*dt)), eta, Physics->eta[ix + iy*Grid->nxEC], G, Physics->G[ix + iy*Grid->nxEC]);

			sigmaII = 2.0*Z*Eff_strainRate;

			//printf("Z = %.2e, sigmaII = %.2e, G[0] = %.2e, sigmay = %.2e\n",Z,sigmaII,Physics->G[0],sigma_y);

			if (sigmaII > sigma_y) {
				//printf("iCell = %i, C = %i\n", iCell, C);
				//khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );
				khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );


				//printf("khi = %.2e, khiCorr = %.2e\n",khi, khiCorr);
				if (khi<0.0) {
					// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
					// if it happens then I consider the case where there are == , which means khi -> inf
					printf("WTF!\n");
					khi = 1e30;
					exit(0);
				}


				Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
			}


			if (Z<Numerics->etaMin) {
				Z = Numerics->etaMin;
			}

			Physics->ZShear[iNode] = Z;
			if (ix == 0 || iy == 0 || ix == Grid->nxS-1 || iy == Grid->nyS-1) {
				Physics->ZShear[iNode] = shearValue(Physics->Z,  ix   , iy, Grid->nxEC);
			}
			*/

			Physics->ZShear[iNode] = shearValue(Physics->Z,  ix   , iy, Grid->nxEC);
#endif

		}

	}

	// 									Shear nodes viscosity
	// ================================================================================

	free(ZprePlasticity);
	free(sigma_y_Stored);

}






void Physics_updateDt(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics)
{
	compute dtAdvOld = Physics->dtAdv;
	compute dtOld = Physics->dt;
	Physics->dtDarcy = 1e100;
	Physics->dtT	 = 1e100;


	Numerics->dtAlphaCorrIni = 1.0;


	Physics->dtAdv 	= Numerics->CFL_fac_Stokes*Grid->dx/(Physics->maxVx); // note: the min(dx,dy) is the char length, so = 1
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Numerics->CFL_fac_Stokes*Grid->dy/(Physics->maxVy));

	if (Numerics->dtVep>0.0) {
		Physics->dt = Numerics->dtVep;
	} else {
		Physics->dt = Physics->dtAdv;
	}


	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt); // dtAdv<=dtVep

	Numerics->lsGoingDown = false; // true if going down, false if going up or staying the same
	Numerics->lsGoingUp = false;



	int iCell;
	compute dtMaxwell_VP_ov_E, dtMaxwell_VP_ov_EP, dtMaxwell_EP_ov_E;
	compute min_dtMaxwell_EP_ov_E = 1e100; // above this time step, effectively elasto-plastic (EP), below is elastic (E)
	compute min_dtMaxwell_VP_ov_E = 1e100;
	compute min_dtMaxwell_VP_ov_EP = 1e100;

	compute eta_vp, eta_ep;


	compute a, b, c;
	compute dtImp, dtExp;
	compute min_dtImp = 1e100;
	compute min_dtExp = 1e100;

	compute maxwellFac = .05;
	compute stress_ep;

	compute eta, G, khi;

#if (DARCY)
	compute eta_b, Bulk, khi_b, phi;
	compute dtExp_b;
	compute min_dtExp_b = 1e100;
	compute dtMaxwell_VP_ov_E_b, dtMaxwell_VP_ov_EP_b, dtMaxwell_EP_ov_E_b;
	compute min_dtMaxwell_EP_ov_E_b = 1e100; // above this time step, effectively elasto-plastic (EP), below is elastic (E)
	compute min_dtMaxwell_VP_ov_E_b = 1e100;
	compute min_dtMaxwell_VP_ov_EP_b = 1e100;
#endif

	compute A = Numerics->dtMaxwellFac_EP_ov_E;
	compute B = Numerics->dtMaxwellFac_VP_ov_EP;
	if (Numerics->timeStep>-1) {
		//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		int iy, ix;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			for (ix = 1; ix<Grid->nxEC-1; ix++) {
				iCell = ix + iy*Grid->nxEC;
				if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]]) {
					eta = Physics->eta[iCell];
					G = Physics->G[iCell];
					khi = Physics->khi[iCell];

					eta_vp = 1.0/(1.0/eta + 1.0/khi);
					stress_ep = 1.0/(1.0/(G)+dtOld/khi);
					eta_ep = 1.0/(1.0/(G*dtOld)+1.0/khi);
					dtMaxwell_EP_ov_E = eta_ep/G;
					//dtMaxwell_VP_ov_E = eta_vp/G;
					dtMaxwell_VP_ov_EP = (eta_vp/stress_ep);

					dtExp = A*dtMaxwell_EP_ov_E + B*dtMaxwell_VP_ov_EP;

					min_dtExp = fmin(min_dtExp, dtExp);

					min_dtMaxwell_EP_ov_E 	= fmin(min_dtMaxwell_EP_ov_E,dtMaxwell_EP_ov_E);
					//min_dtMaxwell_VP_ov_E 	= fmin(min_dtMaxwell_VP_ov_E,dtMaxwell_VP_ov_E);
					min_dtMaxwell_VP_ov_EP 	= fmin(min_dtMaxwell_VP_ov_EP,dtMaxwell_VP_ov_EP);

					/*
					// in implicit form g is the solution of the second order polynomial a*dt^2 + b*dt + c = 0, with
					dtImp = khi*(A*khi + A*eta + 2*B*eta - khi - eta + sqrt((khi + eta)*(A*A*khi + A*A*eta + 4*A*B*eta - 2*A*khi - 2*A*eta + khi + eta)))/(2*G*(-B*eta + khi + eta));
					if (dtImp>0.0) {
						min_dtImp = fmin(min_dtImp, dtImp);
					}
					*/

#if (0)//(DARCY)
					compute eta_b, Bulk, khi_b, phi;
					phi = Physics->phi[iCell];
					eta_b = Physics->eta_b[iCell];
					Bulk = Physics->G[iCell]/sqrt(phi);
					khi_b = Physics->khi_b[iCell];

					compute divV;
					divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
					divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;

					eta_vp = 1.0/(1.0/eta_b + fabs(1.0/khi_b));
					//stress_ep = 1.0/(1.0/(Bulk)+dtOld/khi_b);
					eta_ep =1.0/(1.0/(Bulk*dtOld)+fabs(1.0/khi_b));

					//stress_ep = fabs(Physics->Zb[iCell] * ( - divV + Physics->DeltaP0[iCell]/(Bulk*dtOld) ) ); // Pc
					stress_ep = 1.0/(1.0/(Bulk)+fabs(dtOld/khi_b));
					dtMaxwell_EP_ov_E_b = eta_ep/Bulk;
					//dtMaxwell_VP_ov_E = eta_vp/G;
					dtMaxwell_VP_ov_EP_b = (eta_vp/stress_ep);

					dtExp_b = A*dtMaxwell_EP_ov_E_b + B*dtMaxwell_VP_ov_EP_b;
					min_dtExp_b = fmin(min_dtExp_b, dtExp_b);
					min_dtMaxwell_EP_ov_E_b 	= fmin(min_dtMaxwell_EP_ov_E_b,dtMaxwell_EP_ov_E_b);

#endif


				}
			}
		}


	}
	printf("min_dtExp = %.2e\n", min_dtExp);
	//printf("lastRes = %.2e, absTol = %.2e\n",Numerics->lsLastRes,Numerics->absoluteTolerance);

	if (Numerics->lsLastRes>100.0*Numerics->absoluteTolerance) {
		if (Numerics->timeStep>0) {
			Physics->dt = dtOld;
		}
	} else {

		if (Numerics->use_dtMaxwellLimit) {
			compute coeffA, coeffB, coeffC;

#if (0)//(DARCY)
			Physics->dt  	= fmin(Physics->dt, min_dtExp_b);
#else
			Physics->dt  	= fmin(Physics->dt, min_dtExp);
#endif
			if (fabs((Physics->dt-dtOld)/dtOld)>.02) {
				if (Numerics->timeStep <= 0) {
					Numerics->dtCorr = Physics->dt;
					Numerics->dtAlphaCorr = Numerics->dtAlphaCorrIni;
				}
				Numerics->dtPrevCorr = Numerics->dtCorr;
				Numerics->dtCorr = Physics->dt-dtOld;


				if (Numerics->dtCorr/Numerics->dtPrevCorr<-0.9) {
					Numerics->dtAlphaCorr /= 2.0;
				}

				Physics->dt = dtOld + Numerics->dtAlphaCorr * Numerics->dtCorr;

				// limit the amount of time step decrease from one iteration to another as a factor of dtOld
				/*
				if (Physics->dt/dtOld<0.5) {
					Physics->dt = dtOld * 0.5;
				}
				*/

				if ((Physics->dt-dtOld)<0.0) { 	// going down
					Numerics->lsGoingDown = true;
					//printf("GoingDown!!\n");
				} else { 						// going up
					Numerics->lsGoingUp = true;
					//printf("GoingUp!!\n");
				}

			} else {
				Physics->dt = dtOld;
			}


		}
	}



	// limit the amount upgoing dt as some factor of dt from the last time step;
	//if (Numerics->timeStep>0){
		compute MaxGoingUpFac = 1.5;
		if (Numerics->lsGoingUp && Physics->dt/Numerics->dtPrevTimeStep>MaxGoingUpFac) {
			if (Numerics->timeStep>0)
			Physics->dt = MaxGoingUpFac * Numerics->dtPrevTimeStep;
		}
	//}

	compute dtAdv0 = Physics->dtAdv;
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt);
#if (1)
	if (Numerics->timeStep>0) {
		if (Numerics->use_dtMaxwellLimit) {
#if (0)//(DARCY)
			compute min_EP_ov_E = fmin(min_dtMaxwell_EP_ov_E,min_dtMaxwell_EP_ov_E_b);
			if ((Physics->dtAdv>1.01*min_EP_ov_E)) {
				//Physics->dt = Physics->dtAdv;
			}
#else
			// Security: cannot go lower than EP_ov_E

			if ((Physics->dtAdv>1.01*min_dtMaxwell_EP_ov_E)) {
				Physics->dt = Physics->dtAdv;
			}


#endif
		}
	}
#endif



	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt); // dtAdv<=dtVep

	// Limit according to dtMin, dtMac
	Physics->dtAdv = fmin(Numerics->dtMax,  Physics->dtAdv);
	Physics->dtAdv = fmax(Numerics->dtMin,  Physics->dtAdv);

	Physics->dt = Physics->dtAdv;

	Physics->dtT = Physics->dt;

	if (Numerics->use_dtMaxwellLimit) {
		printf("min_dtExp = %.2e, Numerics->dtAlphaCorr = %.2e, dtAdv0 = %.2e, min_dtMaxwell_EP_ov_E = %.2e, min_dtMaxwell_VP_ov_EP = %.2e, dt = %.2e\n", min_dtExp, Numerics->dtAlphaCorr, dtAdv0, min_dtMaxwell_EP_ov_E, min_dtMaxwell_VP_ov_EP, Physics->dt);
#if (DARCY)
		printf("min_dtExp = %.2e, min_dtExp_b = %.2e, Numerics->dtAlphaCorr = %.2e, dtAdv0 = %.2e, min_dtMaxwell_EP_ov_E = %.2e, dt = %.2e\n", min_dtExp, min_dtExp_b, Numerics->dtAlphaCorr, dtAdv0, min_dtMaxwell_EP_ov_E, Physics->dt);

#endif
	} else {
		printf("dtVep = %.2e, dtAdv = %.2e, dt = %.2e,\n",Numerics->dtVep, Physics->dtAdv, Physics->dt);
	}

	if (fabs(Physics->dt-dtOld)/Physics->dt > 0.01 && Numerics->lsGoingDown) {
		Numerics->oneMoreIt = true;
	} else {
		Numerics->oneMoreIt = false;
	}




}


#if (DARCY)
void Physics_computePerm(Physics* Physics, Grid* Grid, Numerics* Numerics, MatProps* MatProps)
{
	Physics->minPerm = 1E100;
	int iy, ix;
	int iCell;
	compute phi;
	compute phiRef = 0.0001;
	compute PermEffRef = MatProps->perm0_eta_f[0]  *  phiRef*phiRef*phiRef  * ( (1.0-phiRef)*(1.0-phiRef));

	compute perm0;
	SinglePhase* thisPhaseInfo;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		phi = Physics->phi[iCell];

		perm0 = 0.0;
		thisPhaseInfo = Physics->phaseListHead[iCell];
		while (thisPhaseInfo != NULL) {
			perm0 += MatProps->perm0_eta_f[thisPhaseInfo->phase] * thisPhaseInfo->weight;
			thisPhaseInfo = thisPhaseInfo->next;
		}
		perm0 /= Physics->sumOfWeightsCells[iCell];

		Physics->perm_eta_f[iCell] = perm0  *  phi*phi*phi  * ( (1.0-phi)*(1.0-phi));

	}

	Physics_copyValuesToSides(Physics->perm_eta_f, Grid);

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


			Physics->phi[iCell] = Physics->phi0[iCell] + dt*0.5*(    (1.0-Physics->phi0[iCell])*Physics->divV0[iCell] + (1.0-Physics->phi[iCell])*divV   );

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

			sum += Physics->phi[iCell];
		}
	}

	Physics_copyValuesToSides(Physics->phi, Grid);
	Physics_copyValuesToSides(Physics->Dphi, Grid);

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
}





void Physics_computeRho(Physics* Physics, Grid* Grid, MatProps* MatProps)
{

	int iCell;
	SinglePhase* thisPhaseInfo;
#pragma omp parallel for private(iCell, thisPhaseInfo) schedule(dynamic,16)
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->rho[iCell] = 0.0;
		thisPhaseInfo = Physics->phaseListHead[iCell];
		while (thisPhaseInfo != NULL) {
			Physics->rho[iCell] += MatProps->rho0[thisPhaseInfo->phase] * thisPhaseInfo->weight;
			thisPhaseInfo = thisPhaseInfo->next;
		}
		Physics->rho[iCell] /= Physics->sumOfWeightsCells[iCell];

#if (DARCY)

		Physics->rho[iCell] = (1.0 - Physics->phi[iCell])*Physics->rho[iCell] + Physics->phi[iCell]*Physics->rho_f;

#endif

	}

	Physics_copyValuesToSides(Physics->rho, Grid);


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
	int iCell;


	compute scale;

#pragma omp parallel for private(iy, ix, I, iCell, IBC, INeigh, scale) schedule(static,32)
	for (iy = 0; iy<Grid->nyEC; iy++) {
		for (ix = 0; ix<Grid->nxEC; ix++) {
			iCell = ix + iy*Grid->nxEC;
			I = Numbering->map[iCell + INumMap0];
			scale = 1.0;//EqSystem->S[InoDir];

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
	compute contribPhaseAir, contribPhaseWater;
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

			if (phaseAir>-1) {
				contribPhaseAir = contribPhase[phaseAir];
			}else {
				contribPhaseAir = 0.0;
			}

			if (phaseWater>-1) {
				contribPhaseWater = contribPhase[phaseWater];
			}else {
				contribPhaseWater = 0.0;
			}

			if (contribPhaseAir>0) {
				Physics->phase[iCell] = phaseAir;
			} else if (contribPhaseWater>0) {
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
	int nData = 9;
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
			printf("=====  rho  =====\n");
			Data = Physics->rho;
			if (Dim) unit =  kg/m/m/m ;
			break;
		case 5:
			printf("=====  sigma_xx_0  =====\n");
			Data = Physics->sigma_xx_0;
			if (Dim) unit = Pa;
			break;
		case 6:
			printf("=====  Dsigma_xx_0  =====\n");
			Data = Physics->Dsigma_xx_0;
			if (Dim) unit = Pa;
			break;
		case 7:
			printf("=====  sumOfWeightsCells  =====\n");
			Data = Physics->sumOfWeightsCells;
			if (Dim) unit = 1.0;
			break;
		case 8:
			printf("=====  	 P    =====\n");
			Data = Physics->P;
			if (Dim) unit = Pa;
			break;
		case 9:
#if (HEAT)
			printf("=====    T    =====\n");
			Data = Physics->T;
			if (Dim) unit = K;
#endif
			break;
		case 10:
#if (DARCY)
			printf("=====   phi   =====\n");
			Data = Physics->phi;
			if (Dim) unit = 1.0;
#endif
			break;
		case 11:
#if (DARCY)
			printf("=====    Pc    =====\n");
			Data = Physics->Pc;
			if (Dim) unit = Pa;
#endif
			break;
		case 12:
#if (DARCY)
			printf("=====    Pf    =====\n");
			Data = Physics->Pf;
			if (Dim) unit = Pa;
#endif
			break;
		case 13:
#if (DARCY)
			printf("=====    khi_b    =====\n");
			Data = Physics->khi_b;
			if (Dim) unit = Pas;
#endif
			break;
		case 14:
#if (DARCY)
			printf("=====    eta_b    =====\n");
			Data = Physics->eta_b;
			if (Dim) unit = Pas;
#endif
			break;
		case 15:
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

