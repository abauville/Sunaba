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
	compute W = kg*m*m/(s*s*s); 		// Watt

	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	int i;

	// Grid
	// ======================
	Grid->xmin 		/= m;
	Grid->xmax 		/= m;
	Grid->ymin 		/= m;
	Grid->ymax 		/= m;



	// Material properties
	// ======================
	for (i = 0; i < MatProps->nPhase; ++i) {
		MatProps->eta0 [i] 	/= Pas;
		MatProps->rho0 [i] 	/= kg/(m*m*m);
		MatProps->n    [i] 	/= 1.0;
		MatProps->alpha[i]  /= 1.0/K;
		MatProps->beta [i]  /= 1.0/Pa;
		MatProps->k    [i]  /= W/m/K;
		MatProps->G    [i]  /= Pa;
		MatProps->cohesion[i] /= Pa;
		MatProps->frictionAngle[i] /= 1.0;

		MatProps->perm0[i] 	/= m*m;
		MatProps->eta_b[i] 	/= Pas;
		MatProps->B	   [i] 	/= Pa;
	}




	BCStokes->backStrainRate /= 1.0/s;

	BCThermal->TT /= K;
	BCThermal->TB /= K;

	// Physics
	// ======================
	Physics->dt		/= s;
	Physics->g[0] 	/= m/(s*s);
	Physics->g[1] 	/= m/(s*s);

	Physics->epsRef /= 1.0/s;

	Physics->Cp 	/= J/kg/K;


}



