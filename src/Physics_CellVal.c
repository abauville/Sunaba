/*
 * SideValues.c
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */


#include "stokes.h"


void Physics_CellVal_retrieveFromSolution (compute* Val, int ISub, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{

	

	// Where Val is the value to extract from the solution, and DVal the increment since the last time step, IStep is the index of the subsystem of equations
	int I, IBC, INeigh, iy, ix;
	int INumMap0 = Numbering->subEqSystem0Dir[ISub];
	int iCell;


	compute scale;

#pragma omp parallel for private(iy, ix, I, iCell, IBC, INeigh, scale) OMP_SCHEDULE
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

				Val[iCell] = Physics_CellVal_SideValues_getFromBC_Local(EqSystem->x[INeigh]*scale, BC, IBC, ix, iy, Grid);

			}
		}
	}

}

compute Physics_CellVal_SideValues_getFromBC_Local(compute neighValue, BC* BC, int IBC, int ix, int iy, Grid* Grid)
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
		printf("error in Physics_CellVal_retrieveFromSolution: unknown boundary type\n");
		exit(0);
	}
	return sideValue;


}



void Physics_CellVal_SideValues_copyNeighbours_Global(compute* ECValues, Grid* Grid)
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

void Physics_CellVal_SideValues_copyNeighbours_Global_i(int* ECValues, Grid* Grid)
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




void Physics_CellVal_SideValues_getFromBC_Global(compute* ECValues, Grid* Grid, BC* BC, Numbering* Numbering) {
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
		ECValues[I] = Physics_CellVal_SideValues_getFromBC_Local(ECValues[INeigh], BC, IBC, ix, iy, Grid);
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
		ECValues[I] = Physics_CellVal_SideValues_getFromBC_Local(ECValues[INeigh], BC, IBC, ix, iy, Grid);
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
			ECValues[I] = Physics_CellVal_SideValues_getFromBC_Local(ECValues[INeigh], BC, IBC, ix, iy, Grid);

		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			I = ix + iy*Grid->nxEC;

			INeigh =   ix-1 + (iy)*Grid->nxEC  ;
			IBC = abs(Numbering->map[I])-1; // BC nodes are numbered -1 to -n
			ECValues[I] = Physics_CellVal_SideValues_getFromBC_Local(ECValues[INeigh], BC, IBC, ix, iy, Grid);

		}
	}
}





void Physics_CellVal_advectEulerian(compute *A, Model* Model)
{
	
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);



	compute* Anew = (compute*) malloc(Grid->nECTot * sizeof(compute));

	int ix, iy;

	int iC, iN, iS, iW, iE, iVxW, iVxE, iVyS, iVyN;
	compute dAdx_W, dAdx_E, dAdy_S, dAdy_N; 

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	compute dt = Physics->dt;

	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			// Cell indices
			iC = ix   + (iy  )*Grid->nxEC;
			iN = ix   + (iy+1)*Grid->nxEC;
			iS = ix   + (iy-1)*Grid->nxEC;
			iW = ix-1 + (iy  )*Grid->nxEC;
			iE = ix+1 + (iy  )*Grid->nxEC;

			iVxW = ix-1 + iy+Grid->nxVx;
			iVxE = ix   + iy+Grid->nxVx;

			iVyS = ix + (iy-1)*Grid->nxVy;
			iVyN = ix + (iy  )*Grid->nxVy;

			dAdx_W = (A[iC] - A[iW])/dx;
			dAdx_E = (A[iE] - A[iC])/dx;

			dAdy_S = (A[iC] - A[iS])/dy;
			dAdy_N = (A[iN] - A[iC])/dy;

			Anew[iC] = A[iC] + dt* ( - .5*(Physics->Vx[iVxW]*dAdx_W + Physics->Vx[iVxE]*dAdx_E) - .5*(Physics->Vy[iVyS]*dAdy_S + Physics->Vy[iVyN]*dAdy_N) );

		}
	}

	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iC = ix   + (iy  )*Grid->nxEC;
			A[iC] = Anew[iC];
		}
	}

	Physics_CellVal_SideValues_copyNeighbours_Global(A, Grid);



	free(Anew);
	
}
