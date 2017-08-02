/*
 * Interp.c
 *
 *  Created on: Jul 26, 2017
 *      Author: abauville
 */


#include "stokes.h"

#define TEST_SIGMA_INTERP true
#define TEST_SIGMA_INTERP_FROM_PART_TO_CELL true

inline compute Interp_ECVal_Cell2Particle_Local(compute* A, int ix, int iy, int nxEC, compute locX, compute locY)
{
	// Compute a value on particles from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node the particle is attached to
	
	return ( .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxEC]
           + .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxEC]
		   + .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxEC]
		   + .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxEC] );
		   
	/* extra optimization - not so useful
	locX = locX*.25;
	locY = locY*.25;
	compute a = .25-locX;
	compute b = .25+locX;
	compute c = .25-locY;
	compute d = .25+locY;
	return ( a*c*A[ix  +(iy  )*nxEC]
           + a*d*A[ix  +(iy+1)*nxEC]
		   + b*d*A[ix+1+(iy+1)*nxEC]
		   + b*c*A[ix+1+(iy  )*nxEC] );
	*/
}

inline compute Interp_NodeVal_Node2Particle_Local(compute* A, int ix, int iy, int nxS, compute locX, compute locY, int signX, int signY)
{
	// Compute a value on particles from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node the particle is attached to
	return ( .25*(1.0-locX)*(1.0-locY)*A[ix      +(iy  )    *nxS]
		   + .25*(1.0-locX)*(1.0+locY)*A[ix      +(iy+signY)*nxS]
		   + .25*(1.0+locX)*(1.0+locY)*A[ix+signX+(iy+signY)*nxS]
		   + .25*(1.0+locX)*(1.0-locY)*A[ix+signX+(iy  )    *nxS] );
}


inline compute Interp_ECVal_Cell2Node_Local(compute* A, int ix, int iy, int nxEC)
{
	// Compute a value on the shear grid from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy+1)*nxEC] + A[ix+1+(iy+1)*nxEC] + A[ix  +(iy  )*nxEC] + A[ix+1+(iy  )*nxEC])/4;
}

inline compute Interp_NodeVal_Node2Cell_Local(compute* A, int ix, int iy, int nxS)
{
	// Compute a value on an embedded cell center from the A Array of values defined on the shear grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy-1)*nxS] + A[ix-1+(iy-1)*nxS] + A[ix  +(iy  )*nxS] + A[ix-1+(iy  )*nxS])/4;
}

void Interp_All_Particles2Grid_Global(Model* Model)
{


	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	BC* BCStokes 			= &(Model->BCStokes);
	Numbering* NumStokes 	= &(Model->NumStokes);
	Numbering* NumThermal 	= &(Model->NumThermal);
	BC* BCThermal 			= &(Model->BCThermal);





	// Declarations
	// =========================
	int iCell, iNode;
	//int nNeighbours = 4;
	coord locX, locY;


	// Reinitialize Physics array
	bool* changedHead = (bool*) malloc(Grid->nECTot * sizeof(bool));


#pragma omp parallel for private(iCell) OMP_SCHEDULE
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xx_0_ov_G [iCell] = 0.0;
#endif
		Physics->sigma_xx_0 [iCell] = 0.0;
#endif
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


#pragma omp parallel for private(iNode) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
	#if (USE_SIGMA0_OV_G) 
		Physics->sigma_xy_0_ov_G [iNode] = 0.0;
	#endif
		Physics->sigma_xy_0 [iNode] = 0.0;
#endif
		Physics->sumOfWeightsNodes [iNode] = 0.0;
	}


	Physics_PhaseList_reinit(Model);




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
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight, thisPhaseInfo) OMP_SCHEDULE
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

									Physics_Phase_addSingle(&Physics->phaseListHead[iCell],phase);
									thisPhaseInfo = Physics->phaseListHead[iCell];
									break;
								}


							} else {
								thisPhaseInfo = thisPhaseInfo->next;
							}
						}
						thisPhaseInfo->weight += weight;


						// For properties that are stored on the markers, sum contributions
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
						Physics->sigma_xx_0		[iCell] += thisParticle->sigma_xx_0 * weight;
#if (USE_SIGMA0_OV_G)
						Physics->sigma_xx_0_ov_G		[iCell] += (thisParticle->sigma_xx_0/MatProps->G[phase]) * weight;
						//Physics->sigma_xx_0_ov_G		[iCell] += weight / (thisParticle->sigma_xx_0/MatProps->G[phase]);
#endif
						//Physics->sigma_xx_0		[iCell] += weight / thisParticle->sigma_xx_0;
#endif
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
#pragma omp parallel for private(iy, j, iCellS, iCellD) OMP_SCHEDULE
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (j = 0; j<2; ++j) {
				iCellS = j + iy*Grid->nxEC;
				iCellD = Grid->nxEC-2+j + iy*Grid->nxEC;

#if (USE_SIGMA0_OV_G)
				Physics->sigma_xx_0_ov_G		[iCellD] += Physics->sigma_xx_0_ov_G		[iCellS];
				Physics->sigma_xx_0_ov_G		[iCellS]  = Physics->sigma_xx_0_ov_G		[iCellD];
#endif
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
#pragma omp parallel for private(iCell) OMP_SCHEDULE
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {

#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xx_0_ov_G	[iCell] /= Physics->sumOfWeightsCells	[iCell];
		//Physics->sigma_xx_0_ov_G	[iCell] = Physics->sumOfWeightsCells	[iCell] / Physics->sigma_xx_0_ov_G	[iCell];
#endif
		Physics->sigma_xx_0	[iCell] /= Physics->sumOfWeightsCells	[iCell];

#endif
		
		//Physics->sigma_xx_0	[iCell] = Physics->sumOfWeightsCells	[iCell] / Physics->sigma_xx_0	[iCell];
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
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sigma_xx_0, Grid);
#if (HEAT)
	Physics_CellVal_SideValues_getFromBC(Physics->T, Grid, BCThermal, NumThermal);
#endif
#if (DARCY)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->DeltaP0, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->phi0, Grid);
#endif
#if (STRAIN_SOFTENING)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->strain, Grid);
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


	//int iColor; // indexing of the color group for nodes. Nodes of the same color don't collide with each other. i.e. similar to matrix coloring
	int ixStartS[9] = {0,0,0,1,1,1,2,2,2};
	int iyStartS[9] = {0,1,2,0,1,2,0,1,2};
	//SinglePhase* thisPhaseInfo;



	for (iColor = 0; iColor < 9; ++iColor) {
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, signX, signY, phase, i, iNodeNeigh, weight) OMP_SCHEDULE
		for (iy = iyStartS[iColor]; iy < Grid->nyS; iy+=3) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStartS[iColor]; ix < Grid->nxS; ix+=3) { // I don't get why though
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
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
#if (USE_SIGMA0_OV_G)
						Physics->sigma_xy_0_ov_G 		[iNodeNeigh] += (thisParticle->sigma_xy_0 / MatProps->G[phase]) * weight;
						//Physics->sigma_xy_0_ov_G 		[iNodeNeigh] += weight / (thisParticle->sigma_xy_0 / MatProps->G[phase]);
#endif
						Physics->sigma_xy_0 		[iNodeNeigh] += thisParticle->sigma_xy_0 * weight;
#endif
						//Physics->sigma_xy_0 		[iNodeNeigh] += weight / thisParticle->sigma_xy_0;// * weight;
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
#pragma omp parallel for private(iy, iCellS, iCellD,i) OMP_SCHEDULE
		for (iy = 0; iy < Grid->nyS; ++iy) {
			iCellS = 0 + iy*Grid->nxS; // Source
			iCellD = Grid->nxS-1 + iy*Grid->nxS; // destination
#if (USE_SIGMA0_OV_G)
			Physics->sigma_xy_0_ov_G 		[iCellD] += Physics->sigma_xy_0_ov_G  [iCellS];
			Physics->sigma_xy_0_ov_G 		[iCellS]  = Physics->sigma_xy_0_ov_G  [iCellD];
#endif
			Physics->sigma_xy_0 		[iCellD] += Physics->sigma_xy_0  [iCellS];
			Physics->sigma_xy_0 		[iCellS]  = Physics->sigma_xy_0  [iCellD];
			Physics->sumOfWeightsNodes	[iCellD] += Physics->sumOfWeightsNodes[iCellS];
			Physics->sumOfWeightsNodes	[iCellS]  = Physics->sumOfWeightsNodes[iCellD];
		}
	}



#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
	// Dividing by the sum of weights
#pragma omp parallel for private(iNode) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xy_0_ov_G [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic caverage
		//Physics->sigma_xy_0_ov_G [iNode] = Physics->sumOfWeightsNodes[iNode] / Physics->sigma_xy_0_ov_G [iNode]; // harmonic average
#endif
		Physics->sigma_xy_0 [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic caverage
		//Physics->sigma_xy_0 [iNode] = Physics->sumOfWeightsNodes[iNode] / Physics->sigma_xy_0 [iNode]; // harmonic average
	}

#endif


#if (!TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
	Physics_CellVal_advectEulerian(Physics->sigma_xx_0, Model);
#endif



}
































#if (HEAT)
void Interp_Temperature_Grid2Particles_Global(Model* Model)
{
Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps, BC* BCThermal

Grid* Grid 				= &(Model->Grid);
Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	MatProps* MatProps 		= &(Model->MatProps);
	BC* BCThermal 			= &(Model->BCThermal);





	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	compute* DT_sub_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* sumOfWeights_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* DT_rem_OnTheCells = (compute*) malloc( Grid->nECTot *sizeof(compute) );

	int i, iCell;

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i=0;i<4*Grid->nECTot;++i) {
		DT_sub_OnTheCells[i] = 0;
		sumOfWeights_OnTheCells[i] = 0;
	}

	compute TFromCells, DT_sub_OnThisPart, PFromCells;
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
#pragma omp parallel for private(iy, ix, iNode, thisParticle, phase, locX, locY, TFromCells, PFromCells, rhoParticle, dtDiff, DT_sub_OnThisPart, i, iCell, weight) OMP_SCHEDULE
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

				
				TFromCells  = 		Interp_Any_Cell2Particle_Local(Physics->T, ix, iy, Grid->nxEC, locX, locY);
				PFromCells  = 		Interp_Any_Cell2Particle_Local(Physics->P, ix, iy, Grid->nxEC, locX, locY);


				rhoParticle = MatProps->rho0[phase];// * (1+MatProps->beta[phase]*PFromCells) * (1-MatProps->alpha[phase]*thisParticle->T);
				if (rhoParticle<0) {
					printf("error: Negative density on Particles in Physisc_interTempFromCellsParticle\n");
					exit(0);
				}

				dtDiff = (Physics->Cp*rhoParticle)/(  MatProps->k[phase]*( 2.0/(Grid->dx*Grid->dx) + 2.0/(Grid->dy*Grid->dy) )  );


				DT_sub_OnThisPart = ( TFromCells - thisParticle->T ) * ( 1.0 - exp(-d * Physics->dtAdv/dtDiff) );

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
#pragma omp parallel for private(iCell, I, sum, DT_sub_OnThisCell) OMP_SCHEDULE
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

	Physics_CellVal_SideValues_copyNeighbours_Global(DT_rem_OnTheCells, Grid);



	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
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

				
				thisParticle->T  += Interp_ECVal_Cell2Particle_Local(DT_rem_OnTheCells, ix, iy, Grid->nxEC, locX, locY);

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

















void Interp_Phi_Grid2Particles_Global(Model* Model)
{
Grid* Grid 				= &(Model->Grid);
Particles* Particles 	= &(Model->Particles);
Physics* Physics 		= &(Model->Physics);
	

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
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
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

				
				thisParticle->DeltaP0 += Interp_ECVal_Cell2Particle_Local(Physics->DDeltaP, ix, iy, Grid->nxEC, locX, locY);
				}
				thisParticle->phi += Interp_ECVal_Cell2Particle_Local(Physics->DPhi, ix, iy, Grid->nxEC, locX, locY);

#endif

				thisParticle = thisParticle->next;
			}
		}
	}


}



void Interp_Strain_Grid2Particles_Global(Model* Model)
{

Grid* Grid 				= &(Model->Grid);
Particles* Particles 	= &(Model->Particles);
Physics* Physics 		= &(Model->Physics);
	

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

#if (STRAIN_SOFTENING)
	int iCell;
	compute SII;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			Physics_StressInvariant_getLocalCell(Physics, Grid, ix, iy, &SII);
			Physics->Dstrain[iCell] = SII/(2.0*Physics->khi[iCell])*Physics->dtAdv; // Recovering the incremental plastic strain
		}1
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Dstrain, Grid);
#endif

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
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


				thisParticle->strain += Interp_ECVal_Cell2Particle_Local(Physics->DStrain, ix, iy, Grid->nxEC, locX, locY);

#endif

				thisParticle = thisParticle->next;
			}
		}
	}


}





#if (TEST_SIGMA_INTERP)

void Interp_Stresses_Grid2Particles_Global(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	BC* BCThermal 			= &(Model->BCThermal);
	Numbering* NumThermal 	= &(Model->NumThermal);
	Numerics* Numerics 		= &(Model->Numerics);

	int signX, signY;
	compute locX, locY;
	int ix, iy;

	INIT_PARTICLE

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


				thisParticle->sigma_xx_0 += Interp_ECVal_Cell2Particle_Local(Physics->Dsigma_xx_0, ix, iy, Grid->nxEC, locX, locY);
				locX = fabs(locX)-1.0;
				locY = fabs(locY)-1.0;
				thisParticle->sigma_xy_0 += Interp_NodeVal_Node2Particle_Local(Physics->Dsigma_xy_0, ix, iy, Grid->nxS, locX, locY, signX, signY);


				thisParticle = thisParticle->next;
			}
		}
	}

		//END_PARTICLES
	
}





#else
void Interp_Stresses_Grid2Particles_Global(Model* Model)
{
	// General
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	BC* BCThermal 			= &(Model->BCThermal);
	Numbering* NumThermal 	= &(Model->NumThermal);
	Numerics* Numerics 		= &(Model->Numerics);



	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	int signX, signY;


	// See Taras' book pp. 186-187
	compute sigma_xx_0_fromCells;
	compute sigma_xy_0_fromNodes;

	compute d_ve_ini = 0.0;
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
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i=0;i<4*Grid->nSTot;++i) {
		Dsigma_xy_sub_OnTheNodes[i] = 0;
		sumOfWeights_OnTheNodes[i] = 0;
	}
#pragma omp parallel for private(i) OMP_SCHEDULE
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
#pragma omp parallel for private(iy, ix, i, iNode, thisParticle, locX, locY, signX, signY, sigma_xx_0_fromCells, sigma_xy_0_fromNodes, eta, khi, G, dtMaxwell, Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart, iNodeNeigh, weight, iCell) OMP_SCHEDULE
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

				sigma_xx_0_fromCells  = Interp_ECVal_Cell2Particle_Local(Physics->sigma_xx_0, ix, iy, Grid->nxEC, locX, locY);

				eta  				  = Interp_ECVal_Cell2Particle_Local(Physics->eta, ix, iy, Grid->nxEC, locX, locY);
				khi  				  = Interp_ECVal_Cell2Particle_Local(Physics->khi, ix, iy, Grid->nxEC, locX, locY);
				eta_vp = 1.0 / (1.0/eta + 1.0/khi);
				eta_vp = fmax(eta_vp,Numerics->etaMin);
				//printf("eta_vp = %.2e\n",eta_vp);
				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from

				locX = fabs(locX)-1.0;
				locY = fabs(locY)-1.0;


				sigma_xy_0_fromNodes = Interp_NodeVal_Node2Particle_Local(Physics->sigma_xy_0, ix, iy, Grid->nxS, locX, locY, signX, signY);





				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				G = MatProps->G[thisParticle->phase];

				dtMaxwell = eta_vp/G;
				dtMaxwell = fmin(dtm,dtMaxwell);

				// Compute Dsigma sub grid
				Dsigma_xx_sub_OnThisPart =   ( sigma_xx_0_fromCells - thisParticle->sigma_xx_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
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
#pragma omp parallel for private(iNode, I, sum, Dsigma_xy_sub_OnThisNode) OMP_SCHEDULE
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
#pragma omp parallel for private(iCell, I, sum, Dsigma_xx_sub_OnThisCell) OMP_SCHEDULE
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
	Physics_CellVal_SideValues_copyNeighbours_Global(Dsigma_xx_rem_OnTheCells, Grid);







	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, signX, signY) OMP_SCHEDULE
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
					thisParticle->sigma_xx_0  += Interp_ECVal_Cell2Particle_Local(Dsigma_xx_rem_OnTheCells, ix, iy, Grid->nxEC, locX, locY);
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
						
						thisParticle->sigma_xy_0  += Interp_NodeVal_Node2Particle_Local(Dsigma_xy_rem_OnTheNodes, ix, iy, Grid->nxS, locX, locY, signX, signY);
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



#endif
