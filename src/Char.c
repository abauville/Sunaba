/*
 * char.c
 *
 *  Created on: Feb 29, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BCStokes, BC* BCThermal)
{
	// SI units
	compute s 	= Char->time;			// second
	compute m 	= Char->length; 		// meter
	compute kg 	= Char->mass; 			// kilogram
	compute K 	= Char->temperature; 	// Kelvin

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s);


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
		MatProps->eta0 [i] 	/= Char->viscosity;
		MatProps->rho0 [i] 	/= Char->density;
		MatProps->n    [i] 	/= 1.0;
		MatProps->alpha[i]  /= 1/Char->temperature;
		MatProps->beta [i]  /= 1/Char->stress;
		MatProps->k    [i]  /= W/m/K;
		MatProps->G    [i]  /= Char->stress;
		MatProps->cohesion[i] /= Char->stress;
		MatProps->frictionAngle[i] /= 1.0;

		MatProps->kD[i] 	/= m/s;
		MatProps->SD[i] 	/= 1/m;
	}




	BCStokes->backStrainRate /= 1.0/Char->time;

	BCThermal->TT /= Char->temperature;
	BCThermal->TB /= Char->temperature;

	// Physics
	// ======================
	Physics->dt		/= Char->time;
	Physics->g[0] 	/= Char->acceleration;
	Physics->g[1] 	/= Char->acceleration;

	Physics->epsRef /= Char->strainrate;

	Physics->Cp 	/= J/kg/K;


}



