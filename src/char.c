/*
 * char.c
 *
 *  Created on: Feb 29, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BC)
{
	int i;

	// Grid
	// ======================
	Grid->xmin 		/= Char->length;
	Grid->xmax 		/= Char->length;
	Grid->ymin 		/= Char->length;
	Grid->ymax 		/= Char->length;

	Grid->dx   		/= Char->length;
	Grid->dy   		/= Char->length;


	// Material properties
	// ======================
	for (i = 0; i < MatProps->nPhase; ++i) {
		MatProps->eta0[i] 	/= Char->viscosity;
		MatProps->rho0[i] 	/= Char->density;
	}


	// BC
	// ======================
	BC->VxL    		/= Char->velocity;
	BC->VxR    		/= Char->velocity;
	BC->VyB    		/= Char->velocity;
	BC->VyT    		/= Char->velocity;

	BC->VyL    		/= Char->velocity;
	BC->VyR    		/= Char->velocity;
	BC->VxB    		/= Char->velocity;
	BC->VxT    		/= Char->velocity;


	BC->backStrainRate /= 1.0/Char->time;


	// Physics
	// ======================
	Physics->dt		/= Char->time;
	Physics->g[0] 	/= Char->acceleration;
	Physics->g[1] 	/= Char->acceleration;

	Physics->epsRef /= Char->strainrate;


}



