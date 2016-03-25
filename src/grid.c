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

	compute VxL =  BC->backStrainRate*Grid->xmin;
	compute VxR =  BC->backStrainRate*Grid->xmax;
	compute VyB = -BC->backStrainRate*Grid->ymin;
	compute VyT = -BC->backStrainRate*Grid->ymax;

	Grid->xmin += VxL * dt;
	Grid->xmax += VxR * dt;

	Grid->ymin += VyB * dt;
	Grid->ymax += VyT * dt;

	Grid->dx = (Grid->xmax-Grid->xmin)/Grid->nxC;
	Grid->dy = (Grid->ymax-Grid->ymin)/Grid->nyC;

}


