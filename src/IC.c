/*
 * IC.c
 *
 *  Created on: Dec 27, 2016
 *      Author: abauville
 */
#include "stokes.h"


void applyGaussian(compute* Val, IC* IC, Grid* Grid);


#if (HEAT)
void IC_T(Physics* Physics, Grid* Grid, IC* ICThermal, BC* BCThermal)
{

	if (ICThermal->SetupType == IC_HSC) {
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
	} else if (ICThermal->SetupType == IC_Gaussian) {
		applyGaussian(Physics->DT, ICThermal, Grid);
	} else {
		printf("error: unknwon ICThermal->SetupType: %d\n", ICThermal->SetupType);
	}
}
#endif


#if (DARCY)
void IC_phi(Physics* Physics, Grid* Grid, Numerics* Numerics, IC* ICDarcy)
{
	if (ICDarcy->SetupType == IC_Gaussian) {
		applyGaussian(Physics->phi, ICDarcy, Grid);
		int iCell;
		int iy, ix;
		for (iCell=0; iCell<Grid->nECTot; ++iCell) {
			if (Physics->phase[iCell] == Physics->phaseAir || Physics->phase[iCell] == Physics->phaseWater) {
				Physics->phi [iCell] = Numerics->phiMax;
			}
		}

		memcpy(Physics->Dphi, Physics->phi, Grid->nECTot * sizeof(compute));

	} else {
		printf("error in Physics_initPhi: unknwon type\n");
		exit(0);
	}

	/*
	printf("Check phi init\n");
	int iy, ix, iCell;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			printf("%.2e  ", Physics->phi[iCell]);
		}
		printf("\n");
	}
	*/

	//exit(0);



}
#endif


void applyGaussian(compute* Val, IC* IC, Grid* Grid)
{
	// Val is a pointer to data stored on EC, e.g. Physics->T or Physics->phi
	int ix, iy, iCell;
	srand(time(NULL));
	compute noise 		= IC->data[0];
	compute background 	= IC->data[1];
	compute A 			= IC->data[2];
	compute xc 			= IC->data[3];
	compute yc 			= IC->data[4];
	compute wx 			= IC->data[5];
	compute wy 			= IC->data[6];

	compute x 			= Grid->xmin-Grid->DXEC[0]/2.0;
	compute y 			= Grid->ymin-Grid->DYEC[0]/2.0;
	compute XFac = 1.0;
	compute YFac = 1.0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;

			Val[iCell] += background + A*exp(   - XFac* (x-xc)*(x-xc)/(2*wx*wx) - YFac* (y-yc)*(y-yc)/(2*wy*wy)      );
			Val[iCell] += noise*(0.5 - (rand() % 1000)/1000.0);

			if (y==yc) {
				//printf("Physics->Dphi [iCell] = %.2e, x = %.2e, y = %.2e, xc, = %.2e, yc = %.2e, w = %.2e\n",Physics->Dphi [iCell], x, y, xc, yc, w);
			}

			if (ix<Grid->nxEC-1) {
				x += Grid->DXEC[ix];
			} else {
				x = Grid->xmin-Grid->DXEC[0]/2.0;
			}

		}
		if (iy<Grid->nyEC-1) {
			y+=Grid->DYEC[iy];
		}
	}

}

