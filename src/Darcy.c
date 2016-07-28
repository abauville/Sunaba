/*
 * darcy.c
 *
 *  Created on: May 12, 2016
 *      Author: abauville
 */

#include "stokes.h"


/*

void Darcy_setBC(Grid* Grid, Physics* Physics, coord hOcean, PhaseFlag* Phase)
{
	int ix, iy, iCell;
	compute depth, y;

	// Set psi in the air and ocean
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;

			if (Phase[iCell]==Air) {
				Physics->psi[iCell] = 0;
			}

			if (Phase[iCell]==Ocean) {
				y = Grid->ymin - 0.5*Grid->DXEC[ix] + Grid->DYEC[iy]*iy;
				depth = hOcean-y;
				Physics->psi[ix + iy*Grid->nxEC] = depth;
			}


		}
	}

	// Bottom boundary
	iy = 0;
	for (ix = 0; ix < Grid->nxEC; ++ix) {
		iCell = ix+iy*Grid->nxEC;
		Physics->psi[iCell] = Physics->psi[ix   + (iy+1)*Grid->nxEC] + Grid->DYEC[0]; // i.e. gradient = 1
	}

	// Left Boundary
	ix = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		iCell = ix+iy*Grid->nxEC;
		Physics->psi[iCell] = Physics->psi[ix+1+iy*Grid->nxEC]; // i.e. gradient = 0
	}

	// Right Boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		iCell = ix+iy*Grid->nxEC;
		Physics->psi[iCell] = Physics->psi[ix+-1+iy*Grid->nxEC]; // i.e. gradient = 0
	}

}
*/

/*

void Darcy_solve(Darcy* Darcy, Grid* Grid, Physics* Physics, MatProps* MatProps, Particles* Particles)
{
	compute time = 0;

	compute dx;
	compute dy;
	compute dt;

	compute dum1 = max(MatProps->kD, MatProps->nPhase);
	compute kDscale = fmax(dum1,dum1*FAULT_MOD);


	dt =  0.25*fmin(dx*dx,dy*dy)*MatProps->SD[1]/kDscale; // the first 0.25 is a quick and dirty fix for the modified kD



	PhaseFlag* Phase = (PhaseFlag*) malloc(Grid->nECTot * sizeof(PhaseFlag*));

	int ix,iy, iCell;


	int IC, IE, IW, IN, IS;
	compute FluxE, FluxW, FluxN, FluxS;

	Darcy_setPhaseFlag(Phase, Darcy->hOcean, Grid, Particles);



	compute* psiOld = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* psiIni = (compute*) malloc(Grid->nECTot * sizeof(compute));


	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		psiIni[iCell] = Physics->psi[iCell];
	}

	Darcy_setBC(Grid, Physics, Darcy->hOcean, Phase);



	while (time<1*Physics->dt) {
		for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
			psiOld[iCell] = Physics->psi[iCell];
		}



		for (iy = 1; iy < Grid->nyEC-1; ++iy) {
			for (ix = 1; ix < Grid->nxEC-1; ++ix) {
				IC = ix   + (iy)*Grid->nxEC;
				if (Phase[IC]==Solid) {
					dx = Grid->DXEC[ix];
					dy = Grid->DYEC[iy];
					IE = ix+1 + (iy  )*Grid->nxEC;
					IW = ix-1 + (iy  )*Grid->nxEC;
					IN = ix   + (iy+1)*Grid->nxEC;
					IS = ix   + (iy-1)*Grid->nxEC;

					FluxE = 0.5*(Physics->kD[IE]+Physics->kD[IC]) * (psiOld[IE]-psiOld[IC])/dx;
					FluxW = 0.5*(Physics->kD[IW]+Physics->kD[IC]) * (psiOld[IC]-psiOld[IW])/dx;
					FluxN = 0.5*(Physics->kD[IN]+Physics->kD[IC]) * (psiOld[IN]-psiOld[IC]+dy)/dy;
					FluxS = 0.5*(Physics->kD[IS]+Physics->kD[IC]) * (psiOld[IC]-psiOld[IS]+dy)/dy;

					if (iy==1) {
						//printf("GradSouth = %.1e\n", (psiOld[IC]-psiOld[IS]+dy)/dy );
					}

					// Boundary conditions in fluxes at the air/solid and ocean/solid interface
					if (Phase[IE]==Air && psiOld[IC]<0) {
						FluxE = 0;
					}
					if (Phase[IW]==Air && psiOld[IC]<0) {
						FluxW = 0;
					}
					if (Phase[IN]==Air && FluxN>0) {
						if (Darcy->rainFlux<FluxN) {
							FluxN = Darcy->rainFlux;
						}

					}

					Physics->psi[IC] = psiOld[IC] + dt/Physics->SD[IC] * ( (FluxE-FluxW)/dx + (FluxN-FluxS)/dy );
				}
			}
		}





		Darcy_setBC(Grid, Physics, Darcy->hOcean, Phase);

		time += dt;
		//printf("dt = %.2f, time = %.2f, Physics->dt = %.2f\n", dt, time, Physics->dt);


		if (time>Physics->dt-dt) { // update dt on the last iteration to compute up to the correct time
			dt = Physics->dt-time;
		}

	}



	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		//Physics->psi[iCell] = Phase[iCell];
		//printf("Phase[%i] = %.1f\n", iCell, Phase[iCell]);
		Physics->Dpsi[iCell] = Physics->psi[iCell] - psiIni[iCell];

	}




	free(Phase);
	free(psiOld);
	free(psiIni);


}

*/



/*
void Darcy_setPhaseFlag(PhaseFlag* Phase, coord hOcean, Grid* Grid, Particles* Particles)
{
	int ix, iy, iCell, iNode;
	coord depth, y;

	SingleParticle* thisParticle;

	int IxNode[] = {-1,  0, -1, 0};
	int IyNode[] = {-1, -1,  0, 0};
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {

			iCell = ix+iy*Grid->nxEC;



			// Assign Phase
			// ===================

			//Phase[iCell] = Air;
			Phase[iCell] = Solid;


			for (iNode = 0; iNode < 4; ++iNode) {
				thisParticle = Particles->linkHead[ix+IxNode[iNode] + (iy+IyNode[iNode])*Grid->nxS];
				while (thisParticle != NULL && Phase[iCell]==Solid) {

					if (thisParticle->phase==0 || thisParticle->phase==1 ) {
						Phase[iCell] = Air;

					}


					thisParticle = thisParticle->next;
				}


				// updated the ocean particles
				thisParticle = Particles->linkHead[ix+IxNode[iNode] + (iy+IyNode[iNode])*Grid->nxS];
				while (thisParticle != NULL) {

					if (thisParticle->phase==0 || thisParticle->phase==1 ) {
						if (thisParticle->y<hOcean) {
							thisParticle->phase = 1;

						} else {
							thisParticle->phase = 0;
						}
					}
					thisParticle = thisParticle->next;
				}



			}

			//

			// Get Depth
			// ===================
			y = Grid->ymin -0.5*Grid->DYEC[0] + Grid->Y[iy]; // /!\ not so sure about Grid->Y[iy], be careful with node vs EC definition of Y
			depth = -(y-hOcean);
			if (Phase[iCell] == Air && depth>0) {
				Phase[iCell] = Ocean;
			}




		}
	}

	// lower boundary
	iy = 0;
	int I, INeigh;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Phase[I] = Phase[INeigh];
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
		Phase[I] = Phase[INeigh];
	}

	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Phase[I] = Phase[INeigh];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Phase[I] = Phase[INeigh];
	}



}
*/

