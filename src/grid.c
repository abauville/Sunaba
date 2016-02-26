/*
 * grid.c
 *
 *  Created on: Feb 26, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Grid_updatePureShear(Grid* Grid, BC* BC, compute dt)
{
	// update xmin, xmax, ymin, ymax, dx, dy
	// to take into account boundary conditions
	// useful for pure shear

	Grid->xmin += BC->VxL * dt;
	Grid->xmax += BC->VxR * dt;

	Grid->ymin += BC->VyB * dt;
	Grid->ymax += BC->VyT * dt;

	Grid->dx = (Grid->xmax-Grid->xmin)/Grid->nxC;
	Grid->dy = (Grid->ymax-Grid->ymin)/Grid->nyC;
}


