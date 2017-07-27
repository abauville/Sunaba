/*
 * SideValues.c
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */


#include "stokes.h"



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
