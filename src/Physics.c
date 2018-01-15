/*
 * Physics.c
 *
 *  Created on: Feb 24, 2016
 *      Author: abauville
 */

#include "stokes.h"



void Physics_Memory_allocate(Model* Model)
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	


	
	int i;


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

	Physics->lambda 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );

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

	Physics->DeltaP0 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->DDeltaP 		= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
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

#if (USE_SIGMA0_OV_G)
	Physics->sigma_xx_0_ov_G  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0_ov_G	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
#endif
	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->khiShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->GShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->etaShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->ZShear 		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->lambdaShear	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->Eps_pxx = (compute*) malloc(Grid->nECTot*sizeof(compute));
	Physics->Eps_pxy = (compute*) malloc(Grid->nSTot*sizeof(compute));

	Physics->Tau_y = (compute*) malloc(Grid->nECTot*sizeof(compute));
	Physics->Tau_yShear = (compute*) malloc(Grid->nSTot*sizeof(compute));

	Physics->phase 			= (int*) 	malloc( Grid->nECTot * sizeof(int) );


	// Initialize stuff
	//int i;
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0.0;
#if (CRANK_NICHOLSON_VEL || INERTIA)
		Physics->Vx0[i] = 0.0;
#endif
	}
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0.0;
#if (CRANK_NICHOLSON_VEL || INERTIA)
		Physics->Vy0[i] = 0.0;
#endif
	}

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nECTot; ++i) {

		Physics->khi[i] = 1e30;
		Physics->Eps_pxx[i] = 0.0;
#if (STRAIN_SOFTENING)
		Physics->strain[i] = 0.0;
		Physics->Dstrain[i] = 0.0;
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
		Physics->divV0[i] = 0.0;

		Physics->Pf  [i] = 0.0;
		Physics->Pc  [i] = 0.0;
		Physics->DeltaP0 [i] = 0.0;
		Physics->DDeltaP [i] = 0.0;
		Physics->phi [i] = 0.0;
		Physics->phi0[i] = 0.0;



#endif

		Physics->sigma_xx_0[i] = 0.0;
		Physics->Dsigma_xx_0[i] = 0.0;
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xx_0_ov_G[i] = 0.0;
#endif

	}

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->sigma_xy_0[i] = 0.0;
		Physics->Dsigma_xy_0[i] = 0.0;
		Physics->Eps_pxy[i] = 0.0;
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xy_0_ov_G[i] = 0.0;
#endif
	}



	Physics->dtMaxwellMin = 1E+100;
	Physics->dtMaxwellMax = 1E-100;





}


void Physics_Memory_free(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	
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
	free( Physics->lambda );
	free( Physics->khiShear );
	free( Physics->ZShear );
	free( Physics->GShear );
	free(Physics->lambdaShear);

	free( Physics->rho );


	free(Physics->Eps_pxx);
	free(Physics->Eps_pxy);

	free(Physics->Tau_y);
	free(Physics->Tau_yShear);

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

#if (USE_SIGMA0_OV_G)
	free(Physics->sigma_xx_0_ov_G );
	free(Physics->sigma_xy_0_ov_G );
#endif
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





void Physics_Phase_addSingle(SinglePhase** pointerToHead, int phase)
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






void Physics_P_initToLithostatic(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);


	int iy, ix, iCell, iCellS, iCellN, iCellW, iCellE;
	compute rho_g_h = 0.0;

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










void Physics_Velocity_advectEulerian(Model* Model)
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Numbering* NumStokes 	= &(Model->NumStokes);



	int ix, iy;
	compute dVxdx, dVxdy, dVydx, dVydy;
	compute dVxdx0, dVxdy0, dVydx0, dVydy0;
	compute* VxNew = (compute*) malloc(Grid->nVxTot * sizeof(compute));
	compute* VyNew = (compute*) malloc(Grid->nVyTot * sizeof(compute));
	compute Vx, Vy;
	compute dt = Physics->dtAdv;
	for (iy = 1; iy < Grid->nyVx-1; ++iy) {
		for (ix = 1; ix < Grid->nxVx-1; ++ix) {
			dVxdx = (Physics->Vx[ix+1 +  iy   *Grid->nxVx] - Physics->Vx[ix-1 +  iy   *Grid->nxVx])/(2.0*Grid->dx);
			dVxdy = (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   + (iy-1)*Grid->nxVx])/(2.0*Grid->dy);
#if (INERTIA || CRANK_NICHOLSON_VEL)
			dVxdx0 = (Physics->Vx0[ix+1 +  iy   *Grid->nxVx] - Physics->Vx0[ix-1 +  iy   *Grid->nxVx])/(2.0*Grid->dx);
			dVxdy0 = (Physics->Vx0[ix   + (iy+1)*Grid->nxVx] - Physics->Vx0[ix   + (iy-1)*Grid->nxVx])/(2.0*Grid->dy);
#else
			dVxdx0 = dVxdx;
			dVxdy0 = dVxdy;
#endif
			Vy = 0.25* (Physics->Vy[ix   + (iy  )*Grid->nxVy] + Physics->Vy[ix+1 + (iy  )*Grid->nxVy] + Physics->Vy[ix   + (iy-1)*Grid->nxVy] + Physics->Vy[ix+1 + (iy-1)*Grid->nxVy]);
			//VxNew[ix+iy*Grid->nxVx] = Physics->Vx[ix   +  iy   *Grid->nxVx]*(1.0-dt*dVxdx) - dt*Vy*dVxdy;
			VxNew[ix+iy*Grid->nxVx] = Physics->Vx[ix   +  iy   *Grid->nxVx]*(1.0-dt*.5*(dVxdx+dVxdx0)) - dt*Vy*.5*(dVxdy+dVxdy0);
		}
	}

	for (iy = 1; iy < Grid->nyVy-1; ++iy) {
		for (ix = 1; ix < Grid->nxVy-1; ++ix) {
			dVydx = (Physics->Vy[ix+1 +  iy   *Grid->nxVy] - Physics->Vy[ix-1 +  iy   *Grid->nxVy])/(2.0*Grid->dx);
			dVydy = (Physics->Vy[ix   + (iy+1)*Grid->nxVy] - Physics->Vy[ix   + (iy-1)*Grid->nxVy])/(2.0*Grid->dy);
#if (INERTIA || CRANK_NICHOLSON_VEL)
			dVydx0 = (Physics->Vy0[ix+1 +  iy   *Grid->nxVy] - Physics->Vy0[ix-1 +  iy   *Grid->nxVy])/(2.0*Grid->dx);
			dVydy0 = (Physics->Vy0[ix   + (iy+1)*Grid->nxVy] - Physics->Vy0[ix   + (iy-1)*Grid->nxVy])/(2.0*Grid->dy);
#else
			dVydx0 = dVydx;
			dVydy0 = dVydy;
#endif
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
#if (INERTIA || CRANK_NICHOLSON_VEL)
				Physics->Vx0[iVx] = VxNew[iVx];
#endif
			} else {
#if (INERTIA || CRANK_NICHOLSON_VEL)
				Physics->Vx0[iVx] = Physics->Vx[iVx];
#endif
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
#if (INERTIA || CRANK_NICHOLSON_VEL)
				Physics->Vy0[iVy] = VyNew[iVy];
#endif
			} else {
#if (INERTIA || CRANK_NICHOLSON_VEL)
				Physics->Vy0[iVy] = Physics->Vy[iVy];
#endif
			}
		}
	}


	free(VxNew);
	free(VyNew);
}



void Physics_Velocity_retrieveFromSolution(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	BC* BC 					= &(Model->BCStokes);
	Numbering* Numbering 	= &(Model->NumStokes);
	EqSystem* EqSystem		= &(Model->EqStokes);
	


	// Declarations
	// =========================
	int ix, iy, i;
	int I;
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
#pragma omp parallel for private(iy, ix, I, InoDir, IBC, INeigh, scale) OMP_SCHEDULE // maxVx would conflict
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
#pragma omp parallel for private(iy, ix, I, IMap, InoDir, IBC, INeigh, scale) OMP_SCHEDULE // maxVx would conflict
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

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = weight[0]*Physics->Vx[i] + weight[1]*Physics->Vx0[i];
	}
#pragma omp parallel for private(i) OMP_SCHEDULE
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
			maxVx = fmax(maxVx, fabs(Vx));
			maxVy = fmax(maxVy, fabs(Vy));
		}
	}
	Physics->maxVx = maxVx;
	Physics->maxVy = maxVy;



}

#if (CRANK_NICHOLSON_VEL || INERTIA)
void Physics_VelOld_POld_updateGlobal			(Model* Model)
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	
	// A better method would be to intervert the pointers;
	int i;

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx0[i] = Physics->Vx[i];

	}
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy0[i] = Physics->Vy[i];
	}

#if (CRANK_NICHOLSON_P)
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P0[i] = Physics->P[i];
	}
#endif





}
#endif


void Physics_P_retrieveFromSolution(Model* Model)
{
	Physics* Physics 		= &(Model->Physics);
	Grid* Grid 				= &(Model->Grid);
	BC* BCStokes 			= &(Model->BCStokes);
	EqSystem* EqStokes		= &(Model->EqStokes);
	Numbering* NumStokes 	= &(Model->NumStokes);
	


	int iCell;

#if (!DARCY)


	// /!\ For visu it's better if all sides are Neumann
	Physics_CellVal_retrieveFromSolution (Physics->P, 2, Grid, BCStokes, NumStokes, EqStokes);

	// Shift pressure, taking the pressure of the upper left cell (inside) as reference (i.e. 0)
	
	compute RefPressure = 0.0;// = Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];// - 1.0;//Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];//Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];
	
	int ix;
	for (ix=0;ix<Grid->nxEC;++ix) {
		RefPressure += Physics->P[ix+(Grid->nyEC-2)*Grid->nxEC];
	}
	RefPressure/=Grid->nxEC;
	

	
	//compute RefPressure = 0.0;
	//compute RefPressure = Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];// - 1.0;//Physics->P[1 + (Grid->nyEC-2)*Grid->nxEC];//Physics->P[Grid->nxEC/2 + (Grid->nyEC-2)*Grid->nxEC];
	/*
	compute meanP = 0.0;
	compute minP = 1e100;
	compute maxP = -1e100;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		meanP += Physics->P [iCell];
		maxP = fmax(maxP, Physics->P [iCell]);
		minP = fmin(minP, Physics->P [iCell]);
	}
	meanP/= (compute)Grid->nECTot;
	
	RefPressure = meanP;
	//printf("meanP = %.2e, minP = %.2e, maxP = %.2e\n",meanP, minP, maxP);
	*/
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P [iCell] 	= Physics->P [iCell] - RefPressure + Physics->Pback;
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

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P[i] = weight[0]*Physics->P[i] + weight[1]*Physics->P0[i];
	}

	#endif
#endif


#else

	int i;
	Physics_CellVal_retrieveFromSolution (Physics->Pf, 2, Grid, BCStokes, NumStokes, EqStokes);
	Physics_CellVal_retrieveFromSolution (Physics->Pc, 3, Grid, BCStokes, NumStokes, EqStokes);

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

	/*
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
	*/



	// Fill P, the total pressure
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = Physics->Pc[iCell] + Physics->Pf[iCell];
	}



#endif



}


#if (HEAT)
void Physics_T_retrieveFromSolution(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	BC* BCThermal 			= &(Model->BCThermal);
	Numbering* NumThermal 	= &(Model->NumThermal);
	EqSystem* EqThermal  	= &(Model->EqThermal);
	EqSystem* EqThermal  	= &(Model->EqThermal);
	

	Physics_CellVal_retrieveFromSolution (Physics->T, 0, Grid, BCThermal, NumThermal, EqThermal);
}
#endif




void Physics_Sigma0_updateGlobal_fromGrid(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	BC* BC 					= &(Model->BCStokes);
	Numerics* Numerics 		= &(Model->Numerics);


	int ix, iy, iCell, iNode;
#pragma omp parallel for private(iy, ix, iCell) OMP_SCHEDULE
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;
			Physics->sigma_xx_0[iCell] += Physics->Dsigma_xx_0[iCell];
			
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sigma_xx_0, Grid);

#pragma omp parallel for private(iy, ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;
			Physics->sigma_xy_0[iNode] += Physics->Dsigma_xy_0[iNode];
		}
	}








}


void Physics_Dsigma_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	BC* BC 					= &(Model->BCStokes);
	Numerics* Numerics 		= &(Model->Numerics);
	

	// see Taras' book p. 186
	int ix, iy, iCell, iNode;
	compute Z;
	compute Eps_xx, Eps_xy;
	compute dVxdy, dVydx, dVxdx, dVydy;
	compute G;


	compute* EffStrainRate_CellGlobal = (compute*) malloc(Grid->nECTot*sizeof(compute));
	compute* TauII_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));

	Physics_Eta_EffStrainRate_getGlobalCell(Model, EffStrainRate_CellGlobal);

	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			int iCell = ix + iy*Grid->nxEC;
			TauII_CellGlobal[iCell] = 2.0*Physics->Z[iCell] *  EffStrainRate_CellGlobal[iCell];
			//TauII[iCell] -= Physics->Z[iCell]*Physics->Eps_p[iCell];
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(TauII_CellGlobal, Grid);


	compute dt = Physics->dt;
	printf("dt = %.2e, dtaAdv= %.2e\n", Physics->dt, Physics->dtAdv);
	//#pragma omp parallel for private(iy, ix, iCell, dVxdx, dVydy, Eps_xx) OMP_SCHEDULE
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;



			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Eps_xx = 0.5*(dVxdx-dVydy);

			compute Ds0_old = Physics->Dsigma_xx_0[iCell];

#if (USE_SIGMA0_OV_G)
			Physics->Dsigma_xx_0[iCell] = Physics->Z[iCell]*(2.0*Eps_xx + Physics->sigma_xx_0_ov_G[iCell]/(dt)) - Physics->sigma_xx_0[iCell];
#else
			//Physics->Dsigma_xx_0[iCell] = 2.0 * Physics->Z[iCell]*(Eps_xx + Physics->sigma_xx_0[iCell]/(2.0*Physics->G[iCell]*dt)) - Physics->sigma_xx_0[iCell];


			compute SxxVE = 2.0 * Physics->Z[iCell]*(Eps_xx + Physics->sigma_xx_0[iCell]/(2.0*Physics->G[iCell]*dt));
#if (PLASTIC_CORR_RHS)

			//if (Physics->Eps_pxx[iCell]>0.0) { // if yielded
				compute Eps_pxx, SxxVEP;
				compute sign;
				//if (SxxVE>0) { sign = 1.0; } else { sign = -1.0; } // sign of plastic eps should be opposite
				Eps_pxx = Physics->Eps_pxx[iCell];
				SxxVEP = SxxVE - 2.0 * Physics->Z[iCell]*Eps_pxx;


				Physics->Dsigma_xx_0[iCell] = SxxVEP - Physics->sigma_xx_0[iCell];
				//printf("SxxVE = %.2e, SxxVEP = %.2e, Tau_y = %.2e, SIIVE = %.2e Eps_p = %.2e, Epx_xx = %.2e, Eps_pxx = %.2e, SxxVE/SIIVE = %.2e\n", SxxVE, SxxVEP, Physics->Tau_y[iCell], SIIVE, Physics->Eps_p[iCell], Eps_xx, Eps_pxx, SxxVE/SIIVE);
			//} else {
			//	Physics->Dsigma_xx_0[iCell] = SxxVE - Physics->sigma_xx_0[iCell];
			//}
			
#else
			Physics->Dsigma_xx_0[iCell] = SxxVE - Physics->sigma_xx_0[iCell];
#endif


			//Physics->Dsigma_xx_0[iCell] = SxxVE - Physics->sigma_xx_0[iCell];
#endif

#if (USE_UPPER_CONVECTED)
			/*
			// upper convected correction for the rotation of stresses
			compute sigma_xy_0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0,ix,iy,Grid->nxS);
			// Anton's trick
			dVxdy = 0.0;
			compute Sxy_x_Dvxdy = 0.0;
			int iN, Ix, Iy;
			int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
			int IyMod[4] = {0,0,1,1};
			for (iN = 0; iN < 4; ++iN) {
				Ix = (ix-1)+IxMod[iN];
				Iy = (iy-1)+IyMod[iN];
				dVxdy += 0.25*( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
							  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;
				Sxy_x_Dvxdy  +=		  0.25*Physics->sigma_xy_0[Ix+Iy*Grid->nxS]*( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
				- Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;
			}
			
			//Physics->Dsigma_xx_0[iCell] += 2.0 * Physics->Z[iCell]/(Physics->G[iCell])*(Physics->sigma_xx_0[iCell]*dVxdx +  sigma_xy_0*dVxdy );
			Physics->Dsigma_xx_0[iCell] += 2.0 * Physics->Z[iCell]/(Physics->G[iCell])*(Physics->sigma_xx_0[iCell]*dVxdx +  Sxy_x_Dvxdy );
			*/
#endif



			//Physics->Dsigma_xx_0[iCell] *= Physics->dtAdv/Physics->dt; // To update by the right amount according to the time step

			if (isnan(Physics->Dsigma_xx_0[iCell])) {
				printf("isnan Physics->Dsigma_xx_0[iCell]\n");
			}

			if (Numerics->timeStep>0) {
				//Physics->Dsigma_xx_0[iCell] = 0.5*Physics->dtAdv* (Physics->Dsigma_xx_0[iCell]/Physics->dtAdv + Ds0_old/Physics->dtAdv0); // Crank-Nicolson, buggy!!
			}

		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Dsigma_xx_0, Grid);




//#pragma omp parallel for private(iy, ix, iNode,dVxdy, dVydx, Eps_xy, G, Z) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;


			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx] - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy] - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			
			Eps_xy = 0.5*(dVxdy+dVydx);

			//G 	 	= Interp_ECVal_Cell2Node_Local(Physics->G, ix, iy, Grid->nxEC);
			G 		= Physics->GShear[iNode];
			Z 	 	= Physics->ZShear[iNode];
			
			

			compute Ds0_old = Physics->Dsigma_xy_0[iNode];
#if (USE_SIGMA0_OV_G)
			Physics->Dsigma_xy_0[iNode] = Z * (2.0*Eps_xy + Physics->sigma_xy_0_ov_G[iNode]/(dt)) - Physics->sigma_xy_0[iNode];
#else
			//Physics->Dsigma_xy_0[iNode] = 2.0*Z * (Eps_xy + Physics->sigma_xy_0[iNode]/(2.0*G*dt)) - Physics->sigma_xy_0[iNode];


			compute SxyVE = 2.0*Z * (Eps_xy + Physics->sigma_xy_0[iNode]/(2.0*G*dt));
#if (PLASTIC_CORR_RHS)

			//if (Physics->Eps_pxy[iNode]>0.0) { // if yielded

				compute SIIVE =  Interp_ECVal_Cell2Node_Local( TauII_CellGlobal, ix, iy, Grid->nxEC);
				compute Eps_pxy, SxyVEP;
				compute sign;
				//if (SxyVE>0) { sign = 1.0; } else { sign = -1.0; } // sign of plastic eps should be opposite
				Eps_pxy = Physics->Eps_pxy[iNode];
				SxyVEP = SxyVE - 2.0 * Physics->ZShear[iNode]*Eps_pxy;

				Physics->Dsigma_xy_0[iNode] = SxyVEP - Physics->sigma_xy_0[iNode];
				//printf("SxyVE = %.2e, SxyVEP = %.2e, Tau_y = %.2e, SIIVE = %.2e Eps_p = %.2e, Epx_xy = %.2e, Eps_pxy = %.2e, SxyVE/SIIVE = %.2e\n", SxyVE, SxyVEP, Physics->Tau_yShear[iNode], SIIVE, Physics->Eps_pShear[iNode], Eps_xy, Eps_pxy, SxyVE/SIIVE);
			//} else {
			//	Physics->Dsigma_xy_0[iNode] = SxyVE - Physics->sigma_xy_0[iNode];
			//}
			
#else
			Physics->Dsigma_xy_0[iNode] = SxyVE - Physics->sigma_xy_0[iNode];
#endif

#endif	
			






#if (USE_UPPER_CONVECTED)
/*
			compute sigma_xx_0 = Interp_ECVal_Cell2Node_Local(Physics->sigma_xx_0,ix,iy,Grid->nxEC);
			Physics->Dsigma_xy_0[iNode] += 1.0*Z/G * (sigma_xx_0*(dVydx-dVxdy));
*/
#endif	
			

				

			Physics->Dsigma_xy_0[iNode] *= Physics->dtAdv/Physics->dt;


			if (isnan(Physics->Dsigma_xy_0[iNode])) {
				printf("isnan Physics->Dsigma_xy_0[iNode]\n");
			}

			if (Numerics->timeStep>0) {
				//Physics->Dsigma_xy_0[iNode] = 0.5*Physics->dtAdv* (Physics->Dsigma_xy_0[iNode]/Physics->dtAdv + Ds0_old/Physics->dtAdv0); // Crank-Nicolson
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
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;

			phi = Physics->phi[iCell];
			Bulk = Physics->G[iCell]/sqrt(phi);

			divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
			divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
			DeltaP0 = Physics->DeltaP0[iCell];


			Zb 	= Physics->Zb[iCell];

			DeltaP = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

			Physics->DDeltaP[iCell] = DeltaP - Physics->DeltaP0[iCell];
			Physics->DDeltaP[iCell] *= Physics->dtAdv/Physics->dt;
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->DDeltaP, Grid);
#endif



	free(EffStrainRate_CellGlobal);
	free(TauII_CellGlobal);

}



compute Physics_sigma_xxVE_getLocalCell(Model* Model, int ix, int iy) {
	// Where ix and iy are the indices of a Cell
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);

	int iCell = ix + iy*Grid->nxEC;
	compute dt = Physics->dt;

	compute dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
	compute dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

	compute Eps_xx = 0.5*(dVxdx-dVydy);



	return 2.0 * Physics->Z[iCell]*(Eps_xx + Physics->sigma_xx_0[iCell]/(2.0*Physics->G[iCell]*dt));

}


compute Physics_sigma_xyVE_getLocalNode(Model* Model, int ix, int iy) {
	// Where ix and iy are the indices of a Node
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);

	int iNode = ix + iy*Grid->nxS;
	compute dt = Physics->dt;

	compute dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx] - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

	compute dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy] - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;

	compute Eps_xy = 0.5*(dVxdy+dVydx);

	//G 	 	= Interp_ECVal_Cell2Node_Local(Physics->G, ix, iy, Grid->nxEC);
	compute G 		= Physics->GShear[iNode];
	compute Z 	 	= Physics->ZShear[iNode];



	return 2.0*Z * (Eps_xy + Physics->sigma_xy_0[iNode]/(2.0*G*dt)) - Physics->sigma_xy_0[iNode];

}





void Physics_StrainRateInvariant_getLocalCell(Model* Model, int ix, int iy, compute* EII)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	


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


void Physics_StrainRateInvariant_getLocalNode(Model* Model, int ix, int iy, compute* EII)
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	

	// Be careful, Anton's trick not in!!

	compute dVxdy, dVydx, dVxdx, dVydy;

	dVxdy = (Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx]
						 - Physics->Vx[(ix  ) + (iy  )*Grid->nxVx])/Grid->dy;

	dVydx = (Physics->Vy[(ix+1) + (iy  )*Grid->nxVy]
						 - Physics->Vy[(ix  ) + (iy  )*Grid->nxVy])/Grid->dx;


	compute dVxdxCell[4], dVydyCell[4]; // order: NE, NW, SW, SE

	// use Anton's trick for the inner nodes
	if (ix>0 && ix<Grid->nxS-1 && iy>0 && iy<Grid->nyS-1) {
		dVxdxCell[0] = (Physics->Vx[(ix+1)+(iy+1)*Grid->nxVx] - Physics->Vx[(ix  )+(iy+1)*Grid->nxVx])/Grid->dx;
		dVxdxCell[1] = (Physics->Vx[(ix  )+(iy+1)*Grid->nxVx] - Physics->Vx[(ix-1)+(iy+1)*Grid->nxVx])/Grid->dx;
		dVxdxCell[2] = (Physics->Vx[(ix  )+(iy  )*Grid->nxVx] - Physics->Vx[(ix-1)+(iy  )*Grid->nxVx])/Grid->dx;
		dVxdxCell[3] = (Physics->Vx[(ix+1)+(iy  )*Grid->nxVx] - Physics->Vx[(ix  )+(iy  )*Grid->nxVx])/Grid->dx;

		dVydyCell[0] = (Physics->Vy[(ix+1)+(iy+1)*Grid->nxVy] - Physics->Vy[(ix+1)+(iy  )*Grid->nxVy])/Grid->dy;
		dVydyCell[1] = (Physics->Vy[(ix  )+(iy+1)*Grid->nxVy] - Physics->Vy[(ix  )+(iy  )*Grid->nxVy])/Grid->dy;
		dVydyCell[2] = (Physics->Vy[(ix  )+(iy  )*Grid->nxVy] - Physics->Vy[(ix  )+(iy-1)*Grid->nxVy])/Grid->dy;
		dVydyCell[3] = (Physics->Vy[(ix+1)+(iy  )*Grid->nxVy] - Physics->Vy[(ix+1)+(iy-1)*Grid->nxVy])/Grid->dy;
		compute NormalComp_sqr = 0.0;
		int iCell;
		for (iCell = 0; iCell < 4; ++iCell) {

			dVxdx = dVxdxCell[iCell];
			dVydy = dVydyCell[iCell];
			NormalComp_sqr += 0.25*(0.5*(dVxdx-dVydy))*(0.5*(dVxdx-dVydy)) ;

		}
		*EII = sqrt( NormalComp_sqr +  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))   );

	} else {
		if (Grid->isPeriodic) {
			if (ix == 0 || ix == Grid->nxS-1) {
				dVxdx = ( Physics->Vx[(1)+(iy+1)*Grid->nxVx] - Physics->Vx[(Grid->nxVx-1 -1)+(iy+1)*Grid->nxVx] +
						Physics->Vx[(1)+(iy  )*Grid->nxVx] - Physics->Vx[(Grid->nxVx-1 -1)+(iy  )*Grid->nxVx] )/4./Grid->dx;
			}
			else {
				dVxdx = 0.0;
				printf("error in Physics_StrainRateInvariant_getLocalNode. Shouldn't come to this condition");
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

void Physics_StressInvariant_getLocalCell(Model* Model, int ix, int iy, compute* SII) 
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);


	int iCell = ix + iy*Grid->nxEC;


#if (PLASTIC_CORR_RHS)
	int Method = 2;
#else
	int Method = 0; // 0 compute from strain invariant, 1 compute from Dsigma
#endif


	//sigma_xy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0, ix, iy, Grid->nxS);
	if (Method == 0) {
		compute EII;
		compute sq_sigma_xy0,sigma_xx0, sigmaII0;

		compute G, dt, Z;
		compute Eff_strainRate;


		Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);

		// Old stress
		sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
		sigma_xx0     = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];

		sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*sq_sigma_xy0);

		compute dVxdy, dVydx, dVxdx, dVydy, Eps_xx;
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
		compute Exy = 0.0;
		compute dVxdy_av = 0.0;
		compute dVydx_av = 0.0;
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

			dVxdy_av += 0.25*dVxdy;
			dVydx_av += 0.25*dVydx;

			Exy += 0.25*(0.5*(dVxdy+dVydx));

			Exy_x_Sxy0 += 0.25*(0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0[Ix+Iy*Grid->nxS];
		}



		G 		    = Physics->G[iCell];
		Z 			= Physics->Z[iCell];
		dt 			= Physics->dt;
#if (DARCY)
		phi = Physics->phi[iCell];
#endif


		//Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
		
/*
#if (USE_UPPER_CONVECTED)

		compute Exx = Eps_xx;
		compute Txx0 = Physics->sigma_xx_0[iCell];
		compute Txy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0,ix,iy,Grid->nxS);
		dVxdy = dVxdy_av;
		dVydx = dVydx_av;

		//Eff_strainRate = 1.0/(2.0*G*dt) * sqrt(pow((2.0*Exx*G*dt + Txx0 + 2.0*dt*(Txx0*Exx + Txy0*dVxdy)),2.0) + pow((2.0*Exy*G*dt - Txx0*dt*(dVxdy - dVydx) + Txy0),2.0));
		Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(G*dt) + Exy_x_Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
#else
		Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(G*dt) + Exy_x_Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
#endif
*/

		*SII = 2.0*Z*(Eff_strainRate);
	} else if (Method == 1) {
		compute sq_sigma_xy,sigma_xx;
		sq_sigma_xy  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];

		sq_sigma_xy += Physics->Dsigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->Dsigma_xy_0[ix-1+(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->Dsigma_xy_0[ix  +(iy-1)*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->Dsigma_xy_0[ix-1+(iy  )*Grid->nxS];
		sq_sigma_xy += Physics->Dsigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->Dsigma_xy_0[ix  +(iy  )*Grid->nxS];

		sigma_xx     = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];

		*SII = sqrt((sigma_xx)*(sigma_xx)    + 0.25*sq_sigma_xy);


	} else if (Method == 2) {
			


			// Anton's trick
			compute dVxdy = 0;
			compute dVydx = 0;
			compute Exy = 0.0;

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


				Exy += 0.25*(0.5*(dVxdy+dVydx));

			
			}
			// Get Exx
			compute dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			compute dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			compute Exx = 0.5*(dVxdx-dVydy);

			// Get Txx0, Txy0
			compute Txx0 = Physics->sigma_xx_0[iCell];
			compute Txy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0, ix, iy, Grid->nxS);




			// Get friction angle and cohesion
			compute sumOfWeights 	= Physics->sumOfWeightsCells[iCell];
			SinglePhase* thisPhaseInfo;
			int phase;
			compute weight;
			compute cohesion, frictionAngle;
			cohesion = 0.0;
			frictionAngle = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				cohesion 		+= MatProps->cohesion[phase] * weight;
				frictionAngle 	+= MatProps->frictionAngle[phase] * weight;
				thisPhaseInfo = thisPhaseInfo->next;
			}
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;

			compute G = Physics->G[iCell];
			compute Z = Physics->Z[iCell];
			compute Pe = Physics->P[iCell];
			if (Pe<0.0) { // fail safe, avoids  negative yeild stress
				Pe = 0.0;
			}


			compute dt = Physics->dt;

			compute Epxx = Physics->Eps_pxx[iCell];
			compute Epxy = Interp_NodeVal_Node2Cell_Local(Physics->Eps_pxy, ix, iy, Grid->nxS);

			compute Txx = 2.0 * Z*(Exx + Txx0/(2.0*G*dt) - Epxx);
			compute Txy = 2.0 * Z*(Exy + Txy0/(2.0*G*dt) - Epxy);
			
			//if (iCell==150) {
			//printf("Txy_VE0 = %.5e, Txy_VE = %.5e\n", Txy_VE0, Txy_VE);
			//}
			//compute sqr_Txx_VE = Txx_VE*Txx_VE;
			//compute sqr_Txy_VE = Interp_Product_NodeVal_Node2Cell_Local(Txy_VE_CellGlobal, Txy_VE_CellGlobal, ix, iy, nxS);
			*SII = sqrt(Txx*Txx + Txy*Txy);
			//compute TII_VE = 2.0*Physics->Z[iCell]*EffStrainRate_CellGlobal[iCell];
			//compute TII_VE = sqrt(sqr_Txx_VE+sqr_Txy_VE);


	}

}







#if (0)
void Physics_dt_update(Model* Model)
{
	Physics* Physics 		= &(Model->Physics);
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);
	Char* Char 				= &(Model->Char);


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


	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt/2.0); // dtAdv<=dtVep

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
	printf("lastRes = %.2e, absTol = %.2e\n",Numerics->lsLastRes,Numerics->absoluteTolerance);

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
				// Limit the correction to only a small fraction of dt
				// To avoid big jump in time step size from one iteration to the next which potentially causes the residual to increase
				/*
				compute Corr = Numerics->dtAlphaCorr * Numerics->dtCorr;
				compute maxCorr = 0.2*dtOld; 
				if (fabs((Corr))>maxCorr) {
					compute signCorr = Corr/fabs(Corr) ;
					Physics->dt = dtOld + signCorr*maxCorr;
				} else {
					Physics->dt = dtOld + Corr;
				}
				*/
				
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
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt/2.0);
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
			/*
			if ((Physics->dtAdv>1.01*min_dtMaxwell_EP_ov_E)) {
				Physics->dt = Physics->dtAdv;
			}
			*/


#endif
		}
	}
#endif


	/*
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt); // dtAdv<=dtVep

	// Limit according to dtMin, dtMac
	Physics->dtAdv = fmin(Numerics->dtMax,  Physics->dtAdv);
	Physics->dtAdv = fmax(Numerics->dtMin,  Physics->dtAdv);

	Physics->dt = Physics->dtAdv;
	
	
	Physics->dtAdv /= 2.0;

	*/
	Physics->dt = fmin(Numerics->dtMax,  Physics->dt);
	Physics->dt = fmax(Numerics->dtMin,  Physics->dt);
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Physics->dt/2.0); // dtAdv<=dtVep


	/*
	if (Numerics->timeStep<1) {
		Physics->dtAdv	= 1.0*(3600*24*365.25)/Char->time;
	} else {

		
		compute dtCFL 	= Numerics->CFL_fac_Stokes*Grid->dx/(Physics->maxVx); // note: the min(dx,dy) is the char length, so = 1
		dtCFL 	= fmin(dtCFL,  Numerics->CFL_fac_Stokes*Grid->dy/(Physics->maxVy));
		Physics->dtAdv	= fmin(dtCFL,  Physics->dtAdv);
		//Physics->dtAdv	= dtCFL;
		Physics->dtAdv	= fmin(dtCFL, 100000*(3600*24*365.25)/Char->time);
		
		//Physics->dtAdv	= fmin(Physics->dtAdv, 50000*(3600*24*365.25)/Char->time);
		
	
	}
	Physics->dt = 2.0*Physics->dtAdv;
	//Physics->dtT = Physics->dt;
	*/
	

	if (Numerics->use_dtMaxwellLimit) {
		printf("min_dtExp = %.2e, Numerics->dtAlphaCorr = %.2e, dtAdv0 = %.2e, min_dtMaxwell_EP_ov_E = %.2e, min_dtMaxwell_VP_ov_EP = %.2e, dtAdv = %.2e, dt = %.2e\n", min_dtExp, Numerics->dtAlphaCorr, dtAdv0, min_dtMaxwell_EP_ov_E, min_dtMaxwell_VP_ov_EP, Physics->dtAdv, Physics->dt);
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
#else
void Physics_dt_update(Model* Model) {
	// Time step selection based on the analytical solution of the visco-elastic build up

	// Stress build up equation:
	// Sxx = 2*eta*Exx * ( 1 - exp(-G/eta*t) ) 				[1]

	// Derivative wrt time:
	// dSxx/dt = 2*G*Exx * exp(-G/eta*t) 					[2]

	// Solution for the time at a given stress Sxx0:
	// t = eta/G * ln(2*eta*Exx / (2*eta*Exx - Sxx0))
	// can be rewritten:
	// t = tM * ln(1/(1-Sxx0/SxxV_max)
	// with the maxwell time tM = eta/G
	// and the maximum viscous stress (when the viscous strain rate is equal to the total strain rate): SxxV_max = (2*eta*Exx)

	// Let's limit the time step size based on the increment of strain. 
	// The increment of strain is given by:
	// DeltaSxx = dSxx/dt * dt // where dt is the time step size
	// Then from eq. 2:
	// dt = DeltaSxx / (2*G*Exx * exp(-G/eta*t)) 			[3]

	// DeltaSxx can be chosen, for example as a fraction of the maximum viscous stress or of the yield stress.
	// SxxLimit = min(SxxV_max,Syield)
	// DeltaSxx = SxxLimit/n, where n is a non dimensional number representing the fraction of stress

#if (DARCY)
	printf("Time step size selection method not yet adapted to Darcy\n");
	exit(0);
#endif

	

	// Here comes the implementation
	Physics* Physics 		= &(Model->Physics);
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);
	Char* Char 				= &(Model->Char);
	EqSystem* EqStokes 		= &(Model->EqStokes);

	SinglePhase* thisPhaseInfo;
	compute weight, sumOfWeights;
	compute cohesion, frictionAngle;
	compute P;

	compute t;
	compute Sigma_v_max; // maximum viscous stress (if total strain rate = viscous strain rate)
	compute Sigma_yield;
	compute Sigma_limit;
	int phase;

	compute sq_sigma_xy0, sigma_xx0, sigmaII0;
	compute DeltaSigma;
	compute new_dt = 1e200;
	
	
	/*
	if (Physics->time<=2.3 * 1e6 * (3600*24*365)/Char->time) {
		n = 5.0;
	} else {
		n = 100.0;
	}
	*/
	compute dtOld = Physics->dt;
	
	compute EII, sigmaII;


	compute smallest_dt = 1e100;
	compute dt = Physics->dt;
	int ix, iy, iCell;

	compute eta, G;
	compute dAlphaMax = 0.0;

	compute DeltaSigma_Max = 0.0;

	compute DeltaSigma_min = Numerics->deltaSigmaMin;

	//compute stressFac = 1.0;//fmax(0.0,Numerics->dt_stressFac-Numerics->deltaSigmaMin);
	compute stressFac = Numerics->dt_stressFac;


	compute khiLim = 1e29;

	if (Numerics->timeStep<=0) {
		Numerics->dt_DeltaSigma_min_stallFac = 1.0;
	} else {
		if (!fmod(Numerics->stallingCounter+1,5) && EqStokes->normResidual>10.0*Numerics->absoluteTolerance) {
			Numerics->dt_DeltaSigma_min_stallFac/=2.0;
		} else {
			if(Numerics->itNonLin==0) {
				Numerics->dt_DeltaSigma_min_stallFac *= 1.25; // slowly recovers
				Numerics->dt_DeltaSigma_min_stallFac = fmin(1.0,Numerics->dt_DeltaSigma_min_stallFac);
			}
		} 
	}
	Numerics->dt_DeltaSigma_min_stallFac = fmax(Numerics->dt_DeltaSigma_min_stallFac, 1e-3);

	//DeltaSigma_min *= Numerics->dt_DeltaSigma_min_stallFac;

	int it;
	compute dtOld_iter = Physics->dt;

// not simply parallelizable because of smallest_dt
//#pragma omp parallel for private(iy,ix, iCell, eta, G, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, Sigma_limit, cohesion, frictionAngle, thisPhaseInfo, sumOfWeights, phase, weight, P, Sigma_v_max, Sigma_yield, Sigma_limit, sigmaII, DeltaSigma, dt, smallest_dt) OMP_SCHEDULE collapse(2)
	/*
	bool* faultFlag = (bool*) malloc(Grid->nECTot * sizeof(bool));
	for (iy=1;iy<Grid->nyEC-1; ++iy) {
		for (ix=1;ix<Grid->nxEC-1; ++ix) {
			iCell = ix +iy*Grid->nxEC;
			faultFlag[iCell] = false;
		}
	}
	for (iy=1;iy<Grid->nyEC-1; ++iy) {
		for (ix=1;ix<Grid->nxEC-1; ++ix) {
			iCell = ix +iy*Grid->nxEC;
			if (Physics->khi[iCell]<1e29) {
				faultFlag[ix-1 + (iy-1)*Grid->nxEC] = true;
				faultFlag[ix   + (iy-1)*Grid->nxEC] = true;
				faultFlag[ix+1 + (iy-1)*Grid->nxEC] = true;

				faultFlag[ix-1 + (iy  )*Grid->nxEC] = true;
				faultFlag[ix   + (iy  )*Grid->nxEC] = true;
				faultFlag[ix+1 + (iy  )*Grid->nxEC] = true;

				faultFlag[ix-1 + (iy+1)*Grid->nxEC] = true;
				faultFlag[ix   + (iy+1)*Grid->nxEC] = true;
				faultFlag[ix+1 + (iy+1)*Grid->nxEC] = true;
			}
		}
	}
	*/

	int ixLim, iyLim;
	compute P_E, EP_E, V_E, VP_E;
	compute counter = 0;
	compute av_EP_E = 0.0;
	compute minP_E = 1e100;
	compute minEP_E = 1e100;
	compute minV_E = 1e100;
	compute minVP_E = 1e100;
	compute maxEP_E = 0.0;

	bool somethingIsPlastic = false;

	for (iy=1;iy<Grid->nyEC-1; ++iy) {
		for (ix=1;ix<Grid->nxEC-1; ++ix) {
			iCell = ix +iy*Grid->nxEC;
#if (1)
#if (PLASTIC_CORR_RHS)				
			if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]] && Physics->Eps_pxx[iCell] == 0.0) {
#else
			if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]] && Physics->khi[iCell] > khiLim) {
#endif
			//if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]] && !faultFlag[iCell]) {	
			//if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]]) {	
				eta 	= Physics->eta[iCell];
				G 		= Physics->G  [iCell];
				//  Compute sigmaII0
				sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
				sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
				sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
				sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
				sigma_xx0  = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
				sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

				//  Compute sigmaII
				Physics_StressInvariant_getLocalCell(Model, ix, iy, &sigmaII);

				// Get cohesion and frictionAngle
				if (Numerics->timeStep<=0 && Numerics->itNonLin<1) {
					EII = 1.0; // The reference strain in this case is (1/Char.time) / Char.time = 1.0
					Sigma_limit = 2.0*eta*EII/1000.0; 
					
						//printf("Svmax = %.2e, Syield = %.2e, Slimit = %.2e, cohesion = %.2e, frictionAngle = %.2e, P = %.2e\n", Sigma_v_max, Sigma_yield, Sigma_limit, cohesion, frictionAngle, P);
				} else {
					cohesion = 0.0;
					frictionAngle = 0.0;
					thisPhaseInfo = Physics->phaseListHead[iCell];
					sumOfWeights = Physics->sumOfWeightsCells[iCell];
					while (thisPhaseInfo != NULL) {
						phase = thisPhaseInfo->phase;
						weight = thisPhaseInfo->weight;
						cohesion 		+= MatProps->cohesion[phase] * weight;
						frictionAngle 	+= MatProps->frictionAngle[phase] * weight;
						thisPhaseInfo = thisPhaseInfo->next;
					}
					cohesion 		/= sumOfWeights;
					frictionAngle 	/= sumOfWeights;

					P = Physics->P[iCell];

					//  Compute EII
					Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);



					// Get stress limit
					if (0) {
						if (Physics->khi[iCell]<1e29) {
							Sigma_v_max = 2.0*eta*EII;
							Sigma_yield = cohesion*cos(frictionAngle) + P*sin(frictionAngle);
							if (P<0) {
								Sigma_yield  = 0.0;//cohesion * cos(frictionAngle);
							}
							Sigma_limit = fmin(Sigma_v_max,Sigma_yield);
						} else {
							Sigma_limit = 2.0*eta*EII*1.0;
						}
					} else {
						Sigma_v_max = 2.0*eta*EII;
						Sigma_yield = cohesion*cos(frictionAngle) + P*sin(frictionAngle);
						if (P<0) {
							Sigma_yield  = cohesion * cos(frictionAngle);
						}
						
						//Sigma_limit = fmin(Sigma_v_max,Sigma_yield);
						Sigma_limit = Sigma_yield;
						
					}

				}
				
				compute dSigma = fabs(sigmaII - sigmaII0);

				if (sigmaII>Sigma_limit) {
					//printf("SII>Slim!! sigmaII = %.2e, sigma_limit = %.2e, P = %.2e, Sigma_v_max = %.2e,Sigma_yield = %.2e\n", sigmaII, Sigma_limit, P, Sigma_v_max,Sigma_yield);
					sigmaII = Sigma_limit; // because the time step is updated before the viscosity, so stress can be a bit higher than the yield at that moment.
				}
				// Get DeltaSigma
				//DeltaSigma = Sigma_limit*stressFac;
				
				DeltaSigma = DeltaSigma_min;//stressFac * (Sigma_limit-sigmaII)/Sigma_limit   + DeltaSigma_min;
				//DeltaSigma *= Numerics->dt_DeltaSigma_min_stallFac;
				
				new_dt = dtOld * (DeltaSigma/dSigma);
				
				if (new_dt<0) {
					printf("DeltaSigma = %.2e, Sigma_limit = %.2e, sigmaII = %.2e, dSigma = %.2e\n", DeltaSigma, Sigma_limit, sigmaII, dSigma);
					exit(0);
				}

				// compute the corresponding time in the analytical solution
				//t = eta/G * log(2*eta*EII / (2*eta*EII - sigmaII0 ));

				// compute dt using eq. [3]
				//dt = DeltaSigma / (2*G*EII * exp(-G/eta*t));

				if (new_dt<smallest_dt) {
					DeltaSigma_Max = dSigma;
					ixLim = ix;
					iyLim = iy;
					//printf("DeltaSigma = %.2e, dSigma = %.2e, new_dt = %.2e, smallest_dt = %.2e, Physics->dt = %.2e\n",DeltapSigma, dSigma, new_dt, smallest_dt, Physics->dt);
				}
				smallest_dt = fmin(smallest_dt, new_dt);
				
			
			
			}
			
#if (PLASTIC_CORR_RHS)				
			if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]] && Physics->Eps_pxx[iCell] != 0.0) {
#else
			if (MatProps->use_dtMaxwellLimit[Physics->phase[iCell]] && Physics->khi[iCell] < khiLim) {
#endif
				somethingIsPlastic = true;

				EP_E = (1.0/(1.0/(Physics->G[iCell]*Physics->dt) + 1.0/Physics->khi[iCell])) / (Physics->G[iCell]);
				minEP_E = fmin(minEP_E ,EP_E);
				maxEP_E = fmax(maxEP_E ,EP_E);
				av_EP_E += EP_E;
				counter += 1.0;

				P_E = (Physics->khi[iCell]) / (Physics->G[iCell]);
				minP_E = fmin(minP_E ,P_E);

				V_E = (Physics->eta[iCell]) / (Physics->G[iCell]);
				minV_E = fmin(minV_E ,V_E);

				VP_E = (1.0/(1.0/(Physics->eta[iCell]) + 1.0/Physics->khi[iCell])) / (Physics->G[iCell]);
				minVP_E = fmin(minVP_E ,VP_E);
				//printf("VP_E = %.2e, EP_E = %.2e\n",VP_E, EP_E);

			}


#else

			// Limit the stress rotation
				// ==============================================
				compute Sxx, Sxy, Sxx0, Sxy0;
				compute Eps_xx, dVxdx, dVydy, dVxdy, dVydx;

				dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
				dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;
	
				Eps_xx = 0.5*(dVxdx-dVydy);
	
				


				compute Eps_xy = 0.0;
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
					Eps_xy += (0.5*(dVxdy+dVydx))/4.0;
				}
	


				Sxx0 = Physics->sigma_xx_0[iCell];
				Sxx = 2.0 * Physics->Z[iCell]*(Eps_xx + Sxx0/(2.0*Physics->G[iCell]*dtOld));
				
				Sxy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0,ix,iy,Grid->nxS);
				Sxy = 2.0 * Physics->Z[iCell]*(Eps_xy + Sxy0/(2.0*Physics->G[iCell]*dtOld));

				compute alpha, alpha0, dAlpha;

				if (Sxx!=0 && Sxx0!=0) {
					alpha  = 0.5*atan(Sxy /Sxx );
					alpha0 = 0.5*atan(Sxy0/Sxx0);
					dAlpha = alpha-alpha0;
					dAlphaMax = fmax(dAlphaMax,fabs(dAlpha));
					// compute a new dt based by limiting alpha rotation
					if (Numerics->timeStep>2) {
						compute alpha_limit = .25* PI/180.0; // 1 degree max per time step
						//printf("alpha = %.2e deg, alpha0 = %.2e deg, dAlpha = %.2e deg, new_dtAlpha - %.2e, dtOld = %.2e\n", alpha*180.0/PI, alpha0*180.0/PI, dAlpha*180.0/PI, new_dt, dtOld);
						//if (fabs(dAlpha)>alpha_limit) {
							new_dt = dtOld * alpha_limit/fabs(dAlpha);
							if (new_dt<smallest_dt) {
								printf("alpha = %.2e deg, alpha0 = %.2e deg, new_dtAlpha = %.2e, dtOld = %.2e, dAlpha = %.2e, dAlphaMax = %.2e, alpha_limit/fabs(dAlpha) = %.2e\n", alpha*180.0/PI, alpha0*180.0/PI, new_dt, dtOld, dAlpha*180.0/PI, dAlphaMax*180.0/PI, alpha_limit/fabs(dAlpha));
							}
							smallest_dt = fmin(smallest_dt, new_dt);
						//}
					}
				}
				// Limit the stress rotation
				// ==============================================
#endif
		}
	}
	av_EP_E /= counter;
	//Physics->dt = (smallest_dt+Physics->dt)/2.0;
	
	compute dtStress = smallest_dt;
	if (smallest_dt==1e100) { // unlikely case where everything is breaking
		smallest_dt = dtOld;
		printf("The unlikely happened\n");
	}
	
	/*
	if (Numerics->timeStep <= 0) {
		Numerics->dtCorr = dtOld;
		Numerics->dtPrevCorr = Numerics->dtCorr;
		Numerics->dtAlphaCorr = Numerics->dtAlphaCorrIni;
		Physics->dt = dtOld;
	} else {
		
		
		Numerics->dtCorr = Numerics->dtAlphaCorr * (smallest_dt-dtOld);

		//if (fabs(Numerics->dtCorr)/dtOld<0.05) { // avoids small changes 
		//	Numerics->dtCorr = 0.0; 	
		//}
		//printf("Numerics->dtCorr = %.2e, Numerics->dtPrevCorr = %.2e, Ratio = %.2e\n", Numerics->dtCorr, Numerics->dtPrevCorr, Numerics->dtCorr/Numerics->dtPrevCorr);
		if (Numerics->dtCorr/Numerics->dtPrevCorr<-0.9) {
			Numerics->dtAlphaCorr /= 2.0;
		} else {
			Numerics->dtAlphaCorr *= 1.25;
		}
		Numerics->dtAlphaCorr = fmin(Numerics->dtAlphaCorrIni, Numerics->dtAlphaCorr);


		Physics->dt = dtOld + Numerics->dtCorr;
		Numerics->dtPrevCorr = Numerics->dtCorr;
	}

	printf("dtNow = %.2e, 	Numerics->dtCorr = %.2e, smallest_dt = %2e., dtOld = %.2e\n", Physics->dt, 	Numerics->dtCorr, smallest_dt, dtOld);

	//Physics->dt = dtOld;


	*/

	Numerics->lsGoingDown = false;
	Numerics->lsGoingUp = false;
	
	
	
	

	
	



	//if (Numerics->timeStep>40) {
		//Numerics->dtMin = 5e-4;
		//Numerics->dtMax = 5e-4;
	//}
	
	//Physics->dt = fmin(1.01*Numerics->dtPrevTimeStep,Physics->dt);

	compute tol = 0.001;
	printf("(Physics->dt-dtOld)/dtOld = %.2e, dt = %.2e, dtOld = %.2e\n", (Physics->dt-dtOld)/Physics->dt, Physics->dt, dtOld);
	if ((Physics->dt-dtOld)/dtOld<-tol) { 	// going down
		Numerics->lsGoingDown = true;
		printf("going down0\n");
	} else { 						// going up
		Numerics->lsGoingUp = true;
	}

	
	/*
	if (EqStokes->normResidual<10.0*Numerics->absoluteTolerance) { //  don't change the value if the residuals are close to the acceptable solution
		Physics->dt = dtOld;
	}
	*/

	
	Physics->dt = fmin(Numerics->dtMax,  Physics->dt);
	Physics->dt = fmax(Numerics->dtMin,  Physics->dt);

	// dtAdv
	Physics->dtAdv 	= Numerics->CFL_fac_Stokes*Grid->dx/(Physics->maxVx); // note: the min(dx,dy) is the char length, so = 1
	Physics->dtAdv 	= fmin(Physics->dtAdv,  Numerics->CFL_fac_Stokes*Grid->dy/(Physics->maxVy));
	compute dtAdvAlone = Physics->dtAdv;
	Physics->dtAdv 	= fmin(Physics->dtAdv, Physics->dt);
	//Physics->dtAdv 	= fmax(Physics->dtAdv, 0.001*dtAdvAlone);




	compute alpha_lim = 5.0*PI/180.0;
	int iNode;
	compute dtRot;
	compute dtRotMin = 1e100;
	compute omega;
	if (Numerics->timeStep>0) {	 
		// Compute the Alpha array
		// add a condi	ztion with signX signY to avoid recomputing alpha if not necessary
		
	#pragma omp parallel for private(iy, ix, iNode) OMP_SCHEDULE
		for (iy=0; iy<Grid->nyS; iy++) {
			for (ix=0; ix<Grid->nxS; ix++) {
				iNode = ix + iy*Grid->nxS;
				omega  = .5*((Physics->Vy[ix+1 + (iy  )*Grid->nxVy] - Physics->Vy[ix   +(iy  )*Grid->nxVy])/Grid->DXEC[ix]
							- (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   +(iy  )*Grid->nxVx])/Grid->DYEC[iy]);
				dtRot = alpha_lim/fabs(omega);
				dtRotMin = fmin(dtRotMin,dtRot);
			}
		}
		
		Physics->dtAdv = fmin(dtRotMin,Physics->dtAdv);
	}

//#if (!PLASTIC_CORR_RHS)
	if (Numerics->timeStep>5) {
		//if (EP_E<1e100) {
		if (counter>0.0) {
			//Physics->dtAdv = fmax(Physics->dtAdv,1.0*EP_E); // Avoid entering the dominantly elastic domain
			//Physics->dtAdv = fmax(Physics->dtAdv,1.0*maxEP_E); // Avoid entering the dominantly elastic domain
		}
	}
//#endif


	Physics->dtAdv = fmin(2.0*dtOld,  Physics->dtAdv);
	Physics->dtAdv = fmax(0.5*dtOld,  Physics->dtAdv);

	if (somethingIsPlastic) {
		compute dtPFac = 0.5;
		compute dtPlastic = dtPFac*minEP_E+(1.0-dtPFac)*minVP_E;
		Physics->dtAdv = fmin(Physics->dtAdv,dtPlastic);
	}

	Physics->dtAdv = fmin(Numerics->dtMax,  Physics->dtAdv);
	Physics->dtAdv = fmax(Numerics->dtMin,  Physics->dtAdv);

	

#if (ADV_INTERP) 
	
	Physics->dt = Physics->dtAdv;
	
#else
	Physics->dtAdv = Physics->dt;
#endif


	

	
	//Physics->dt = 10.0*Physics->dtAdv;
	//printf("limiting cell: ix = %i, iy = %i \n", ixLim, iyLim);
	//printf("scaled_dt = %.2e yr, dtMin = %.2e, dtMax = %.2e, DeltaSigma_min = %.2e MPa, DeltaSigma_Max = %.2e MPa,  dt_DeltaSigma_min_stallFac = %.2e, Numerics->dtAlphaCorr = %.2e, dAlphaMax = %.1f deg, dtStress = %.2e, dtAdvAlone = %.2e, Physics->dt = %.2e\n", Physics->dt*Char->time/(3600*24*365.25), Numerics->dtMin, Numerics->dtMax, Numerics->dt_DeltaSigma_min_stallFac*DeltaSigma_min *Char->stress/1e6 , DeltaSigma_Max*Char->stress/1e6,  Numerics->dt_DeltaSigma_min_stallFac, Numerics->dtAlphaCorr, dAlphaMax*180.0/PI , dtStress, dtAdvAlone, Physics->dt);
	printf("scaled_dt = %.2e yr, dtMin = %.2e, dtMax = %.2e, DeltaSigma_min = %.2e MPa, DeltaSigma_Max = %.2e MPa,  dt_DeltaSigma_min_stallFac = %.2e, Numerics->dtAlphaCorr = %.2e, dtStress = %.2e, dtAdvAlone = %.2e, dtRotMin = %.2e, Physics->dt = %.2e\n", Physics->dt*Char->time/(3600*24*365.25), Numerics->dtMin, Numerics->dtMax, Numerics->dt_DeltaSigma_min_stallFac*DeltaSigma_min *Char->stress/1e6 , DeltaSigma_Max*Char->stress/1e6,  Numerics->dt_DeltaSigma_min_stallFac, Numerics->dtAlphaCorr, dtStress, dtAdvAlone, dtRotMin, Physics->dt);

	printf("minEP/E = %.2e yr, maxEP/E = %.2e yr, avEP_E = %.2e, P/E = %.2e yr, V/E = %.2e yr, VP/E = %.2e yr\n", minEP_E*Char->time/(3600*24*365.25), maxEP_E*Char->time/(3600*24*365.25), av_EP_E*Char->time/(3600*24*365.25), minP_E*Char->time/(3600*24*365.25), minV_E*Char->time/(3600*24*365.25), minVP_E*Char->time/(3600*24*365.25));

	//free(faultFlag);


}
#endif

#if (DARCY)
void Physics_Perm_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);
	MatProps* MatProps 		= &(Model->MatProps);
	


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

	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->perm_eta_f, Grid);

}


void Physics_Phi_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);
	


	int iy, ix;
	int iCell;
	compute dt = Physics->dtAdv;
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

	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->phi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Dphi, Grid);

}







#endif





void Physics_Rho_updateGlobal(Model* Model)
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);

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

	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->rho, Grid);


}


/*
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
*/








void Physics_Phase_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Particles* Particles 	= &(Model->Particles);
	MatProps* MatProps 		= &(Model->MatProps);
	


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


	Physics_CellVal_SideValues_copyNeighbours_Global_i(Physics->phase,Grid);

}



void Physics_PhaseList_reinit(Model* Model) 
{

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	

	int iCell;
	SinglePhase* temp;
#pragma omp parallel for private(iCell, temp) OMP_SCHEDULE
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


void Physics_check(Model* Model) 
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);

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
#if (HEAT)
	compute K 	= Char->temperature; 	// Kelvin

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s); 		// Watt
#endif
	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	//compute mol = 1.0;



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

void Physics_NodeVal_advectEulerian(compute *A, Model* Model)
{
	
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);



	compute* Anew = (compute*) malloc(Grid->nSTot * sizeof(compute));

	int ix, iy;

	int iC, iN, iS, iW, iE, iVxN, iVxS, iVyW, iVyE;
	compute dAdx_W, dAdx_E, dAdy_S, dAdy_N; 

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	compute dt = Physics->dt;

	compute Vx, Vy;

	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			// Cell indices
			iC   = ix   + (iy  )*Grid->nxS;
			iN   = ix   + (iy+1)*Grid->nxS;
			iS   = ix   + (iy-1)*Grid->nxS;
			iW   = ix-1 + (iy  )*Grid->nxS;
			iE   = ix+1 + (iy  )*Grid->nxS;

			iVxS = ix   + (iy  )+Grid->nxVx;
			iVxN = ix   + (iy+1)+Grid->nxVx;

			iVyW = ix   + (iy  )*Grid->nxVy;
			iVyE = ix+1 + (iy  )*Grid->nxVy;

			dAdx_W = (A[iC] - A[iW])/dx;
			dAdx_E = (A[iE] - A[iC])/dx;

			dAdy_S = (A[iC] - A[iS])/dy;
			dAdy_N = (A[iN] - A[iC])/dy;

			Vx = .5*( Physics->Vx[iVxS] + Physics->Vx[iVxN]);
			Vy = .5*( Physics->Vy[iVyW] + Physics->Vy[iVyE]);

			Anew[iC] = A[iC] + dt* ( - Vx * .5*(dAdx_W + dAdx_E) - Vy*.5*(dAdy_S + dAdy_N) );

		}
	}

	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			iC = ix   + (iy  )*Grid->nxS;
			A[iC] = Anew[iC];
		}
	}

	// Values should be copied as well
	//Physics_CellVal_SideValues_copyNeighbours_Global(A, Grid);



	free(Anew);
	
}
