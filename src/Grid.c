/*
 * grid.c
 *
 *  Created on: Feb 26, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Grid_allocateMemory(Grid* Grid)
{
	Grid->X = (compute*) malloc(Grid->nxS * sizeof(compute) );
	Grid->Y = (compute*) malloc(Grid->nyS * sizeof(compute) );

	Grid->DXS = (compute*) malloc( (Grid->nxS-1) * sizeof(compute) );
	Grid->DYS = (compute*) malloc( (Grid->nyS-1) * sizeof(compute) );

	Grid->DXEC =  (compute*) malloc( (Grid->nxEC-1) * sizeof(compute) );
	Grid->DYEC =  (compute*) malloc( (Grid->nyEC-1) * sizeof(compute) );

}

void Grid_freeMemory(Grid* Grid)
{
	free(Grid->X);
	free(Grid->Y);

	free(Grid->DXS);
	free(Grid->DYS);

	free(Grid->DXEC);
	free(Grid->DYEC);

	free(Grid->nxPerSeg);
	free(Grid->nyPerSeg);

}


void Grid_init(Grid* Grid, Input* Input, Numerics* Numerics)
{
	Grid->userDefined = false;



	// Case for regular grid
	if (Grid->userDefined==false) {

		int ix, iy;
		compute dx = (Grid->xmax-Grid->xmin)/Grid->nxC;
		compute dy = (Grid->ymax-Grid->ymin)/Grid->nyC;
		Numerics->dLmin = fmin(dx,dy);


		Grid->Y[0] = Grid->ymin;

		for (iy = 0; iy < Grid->nyS-1; ++iy) {
			Grid->DYS[iy] = dy;
			Grid->Y[iy+1] = Grid->Y[iy] + dy;
		}

		Grid->DYEC[0] = dy;
		for (iy = 0; iy < Grid->nyS-2; ++iy) {
			Grid->DYEC[iy+1] = Grid->DYS[iy]/2.0 + Grid->DYS[iy+1]/2.0;
		}
		Grid->DYEC[Grid->nyEC-2] = Grid->DYS[Grid->nyS-2];




		Grid->X[0] = Grid->xmin;
		for (ix = 0; ix < Grid->nxS-1; ++ix) {
			Grid->DXS[ix] = dx;
			Grid->X[ix+1] = Grid->X[ix] + dx;
		}

		Grid->DXEC[0] = dx;
		for (ix = 0; ix < Grid->nxS-2; ++ix) {
			Grid->DXEC[ix+1] = Grid->DXS[ix]/2.0 + Grid->DXS[ix+1]/2.0;
		}
		Grid->DXEC[Grid->nxEC-2] = Grid->DXS[Grid->nxS-2];



		/*
		printf("Grid->X\n");
		for (ix = 0; ix < Grid->nxS; ++ix) {
			printf("%.2f  ",Grid->X[ix] );
		}
		printf("\n");

		printf("Grid->DXS\n");
		for (ix = 0; ix < Grid->nxS-1; ++ix) {
			printf("%.2f  ",Grid->DXS[ix] );
		}
		printf("\n");

		printf("Grid->DYS\n");
		for (iy = 0; iy < Grid->nyS-1; ++iy) {
			printf("%.2f  ",Grid->DYS[iy] );
		}
		printf("\n");

		printf("Grid->DXEC\n");
		for (ix = 0; ix < Grid->nxS; ++ix) {
			printf("%.2f  ",Grid->DXEC[ix] );
		}
		printf("\n");

		printf("Grid->Y\n");
		for (iy = 0; iy < Grid->nyS; ++iy) {
			printf("%.2f  ",Grid->Y[iy] );
		}
		printf("\n");

		printf("Grid->DYEC\n");
		for (iy = 0; iy < Grid->nyS; ++iy) {
			printf("%.2f  ",Grid->DYEC[iy] );
		}
		printf("\n");
	*/
		//exit(0);


	}
}


void Grid_updatePureShear(Grid* Grid, BC* BC, compute dt)
{
	// update xmin, xmax, ymin, ymax, dx, dy
	// to take into account boundary conditions
	// useful for pure shear

	compute VxL =  BC->backStrainRate*Grid->xmin;
	compute VxR =  BC->backStrainRate*Grid->xmax;
	compute VyB = -BC->backStrainRate*Grid->ymin;
	compute VyT = -BC->backStrainRate*Grid->ymax;

	int iy, ix;
	compute locY, locX;

	Grid->Y[0] += VyB * dt;

	for (iy = 0; iy < Grid->nyS-1; ++iy) {
		locY = (Grid->Y[iy+1]+Grid->ymin)/(Grid->ymax-Grid->ymin);
		Grid->Y[iy+1] += (1.0-locY)*VyB + locY*VyT;
		Grid->DYS[iy] = Grid->Y[iy+1] - Grid->Y[iy];
	}

	Grid->DYEC[0] = Grid->DYS[0];
	for (iy = 0; iy < Grid->nyS-2; ++iy) {
		Grid->DYEC[iy+1] = Grid->DYS[iy]/2.0 + Grid->DYS[iy+1]/2.0;
	}
	Grid->DYEC[Grid->nyEC-2] = Grid->DYS[Grid->nyS-2];




	Grid->X[0] += VxL * dt;
	for (ix = 0; ix < Grid->nxS-1; ++ix) {
		locX = (Grid->X[ix+1]+Grid->xmin)/(Grid->xmax-Grid->xmin);
		Grid->X[ix+1] += (1.0-locX)*VxL + locX*VxR;
		Grid->DXS[ix] = Grid->X[ix+1] - Grid->X[ix];
	}

	Grid->DXEC[0] = Grid->DXS[0];
	for (ix = 0; ix < Grid->nxS-2; ++ix) {
		Grid->DXEC[ix+1] = Grid->DXS[ix]/2.0 + Grid->DXS[ix+1]/2.0;
	}
	Grid->DXEC[Grid->nxEC-2] = Grid->DXS[Grid->nxS-2];


	Grid->xmin += VxL * dt;
	Grid->xmax += VxR * dt;

	Grid->ymin += VyB * dt;
	Grid->ymax += VyT * dt;

}


