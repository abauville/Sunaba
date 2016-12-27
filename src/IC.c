/*
 * IC.c
 *
 *  Created on: Dec 27, 2016
 *      Author: abauville
 */
#include "stokes.h"

#if (HEAT)
void IC_T(Physics* Physics, Grid* Grid, IC* ICThermal, BC* BCThermal)
{
	srand(time(NULL));
	int iCell, iy, ix;
	compute y, Tm, Kappa, age, noise;
	noise 	= ICThermal->data[0];
	Tm 		= ICThermal->data[1] - BCThermal->TT ;
	age 	= ICThermal->data[2];
	printf("Tm = %.2e, age = %.2e, noise = %.2e\n",Tm, age, noise);
	printf("Grid->ymin = %.2e, Grid->ymax = %.2e\n", Grid->ymin, Grid->ymax);
	compute norm_g = sqrt(Physics->g[0]*Physics->g[0] + Physics->g[1]*Physics->g[1]);
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		y = -(Grid->Y[iy] + Grid->DYEC[0]/2.0);

		if (y<0.0) {
			y = 0.0;
		}

		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;

			Kappa = Physics->k[iCell]/(Physics->rho_g[iCell]/norm_g * Physics->Cp);

			Physics->T[iCell] = Tm * erf( y/(2*sqrt(Kappa*age))   ) + BCThermal->TT;

			Physics->T[iCell] += noise*(0.5 - (rand() % 1000)/1000.0);

			Physics->DT[iCell] = Physics->T[iCell];

		}
		//printf("y = %.2e, Physics->T[iCell] = %.2e\n",y, Physics->T[iCell]);
	}

	Physics_copyValuesToSides(Physics->T, Grid);

}
#endif
