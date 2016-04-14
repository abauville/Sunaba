/*
 * Physics.c
 *
 *  Created on: Feb 24, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal, BC* BCThermal)
{

	// Declarations
	// =========================
	int iCell, i;
	int nNeighbours = 4;
	coord locX, locY;

	coord dx = Grid->dx;
	coord dy = Grid->dy;
	compute* sumOfWeights 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* eta0 			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* n    			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* rho  			= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* k    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* G    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* T    		  	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* cohesion 		= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* frictionAngle 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* sigma_xx_0   	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* sigma_xy_0   	= (compute*) malloc(nNeighbours * Grid->nSTot * sizeof(compute));




	// Reinitialize Physics array
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		eta0[i] = 0;
		n[i] = 0;
		rho[i] = 0;
		T  [i] = 0;
		k  [i] = 0;
		G [i] = 0;
		sigma_xx_0 [i] = 0;
		cohesion[i] = 0;
		frictionAngle[i] = 0;


		sumOfWeights[i] = 0;
	}

	for (i = 0; i < nNeighbours * Grid->nSTot; ++i) {
		sigma_xy_0 [i] = 0;
	}

	// For debugging purpose only
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->sigma_xx_0[i] = -1;
	}
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->sigma_xy_0[i] = -1;
	}

	printf("Init ok\n");


	compute weight;


	//int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	//int IyMod[4] = {0,0,1,1};

	int phase;

	int nxEC = Grid->nxEC;
	int xMod[4], yMod[4], ix,  iy;

	int I, C;

	xMod[0] = -1; yMod[0] = -1;
	xMod[1] =  1; yMod[1] = -1;
	xMod[2] = -1; yMod[2] =  1;
	xMod[3] =  1; yMod[3] =  1;



	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right




	SingleParticle* thisParticle = NULL;
	// Loop through inner cells
	// ========================
	int iNode = 0;
	//int Count = 0;
//#pragma omp parallel for private(ix, iNode, thisParticle, locX, locY, phase, Ix, Iy, i, locEta, iNode, dVxdy, dVydx, dVxdx, dVydy, locEps_II) schedule(static,32)







	compute* StrainRateInvariant = (compute*) malloc(Grid->nECTot * sizeof(compute));



	Physics_computeStrainRateInvariant(Physics, Grid, StrainRateInvariant);








	//printf("=== Part Temp ===\n");
	// Loop through inner nodes
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight) schedule(static,32)
	for (iy = 1; iy < Grid->nyS-1; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 1; ix < Grid->nxS-1; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				phase = thisParticle->phase;
				//printf("phase = %i, locX= %.2f, locY=%.2f  \n", thisParticle->phase,locX,locY);
				// Loop through neighbours
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);

					weight = fabs((locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5));
					weight = sqrt(weight);
					//weight = 1-sqrt((1-fabs(locX))*(1-fabs(locX))    +      (1-fabs(locY))*(1-fabs(locY))) ;
					//if (weight<0)
					//	weight = 0;
					//weight = (1-weight)*(1-weight);

					eta0			[iCell*4+i] += 1/MatProps->eta0[phase] * weight;
					n				[iCell*4+i] += MatProps->n   [phase] * weight;
					rho				[iCell*4+i] += MatProps->rho0[phase] * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
					k				[iCell*4+i] += MatProps->k   [phase] * weight;
					G				[iCell*4+i] += 1/MatProps->G [phase] * weight; // harmonic average
					cohesion		[iCell*4+i] += 1/MatProps->cohesion[phase] * weight;
					frictionAngle	[iCell*4+i] += MatProps->frictionAngle[phase] * weight;
					T 				[iCell*4+i] += thisParticle->T * weight;
					sigma_xx_0 		[iCell*4+i] += thisParticle->sigma_xx_0 * weight;
					sumOfWeights	[iCell*4+i] += weight;


				}
				thisParticle = thisParticle->next;
			}
		}
	}




	printf("=== eta fill ===\n");

	compute sum;
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			I = 4*iCell;
			sum = sumOfWeights[I+0] + sumOfWeights[I+1] + sumOfWeights[I+2] + sumOfWeights[I+3];
			//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
			if (sum==0) {
				printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
				exit(0);
			}

			Physics->T  [iCell] =     (   T[I+0] +    T[I+1] +    T[I+2] +    T[I+3]) / sum;
			Physics->eta0[iCell]= sum/(eta0[I+0] + eta0[I+1] + eta0[I+2] + eta0[I+3]) ; // harmonic average
			Physics->n  [iCell] =     (   n[I+0] +    n[I+1] +    n[I+2] +    n[I+3]) / sum;
			Physics->rho[iCell] =     ( rho[I+0] +  rho[I+1] +  rho[I+2] +  rho[I+3]) / sum;
			Physics->k  [iCell] =     (   k[I+0] +    k[I+1] +    k[I+2] +    k[I+3]) / sum;
			Physics->G  [iCell] = sum/(   G[I+0] +    G[I+1] +    G[I+2] +    G[I+3]) ; // harmonic average
			Physics->cohesion   [iCell] = sum/( cohesion  [I+0] + cohesion  [I+1] + cohesion  [I+2] + cohesion  [I+3]);
			Physics->frictionAngle[iCell] = ( frictionAngle[I+0] + frictionAngle[I+1] + frictionAngle[I+2] + frictionAngle[I+3]) / sum;
			Physics->sigma_xx_0 [iCell] = ( sigma_xx_0[I+0] + sigma_xx_0[I+1] + sigma_xx_0[I+2] + sigma_xx_0[I+3]) / sum ; // harmonic average

			//printf("Physics->simga_xx_0[%i] %.2e, sigma_xx_0 = %.2f %.2f %.2f %.2f, sum = %.2f\n", iCell, Physics->sigma_xx_0[iCell], sigma_xx_0[I+0], sigma_xx_0[I+1], sigma_xx_0[I+2], sigma_xx_0[I+3], sum);


		}
	}
	printf("=== end fill ===\n");

	/*
	printf("=== Check sigma_xx 0 ===\n");
	C = 0;
	//int ix, iy;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.3f  ", Physics->sigma_xx_0[C]);
			C++;

		}
		printf("\n");
	}
	*/



	// Replace boundary values by their neighbours
	int INeigh, IBC;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];
		Physics->k  [I] = Physics->k  [INeigh];
		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}
		}

	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];
		Physics->k  [I] = Physics->k  [INeigh];
		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];

		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}
		}



	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {

		I = ix + iy*Grid->nxEC;
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];
		Physics->k  [I] = Physics->k  [INeigh];
		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];


		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}
		}


	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		IBC = abs(NumThermal->map[I])-1; // BC nodes are numbered -1 to -n
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Physics->eta0[I] = Physics->eta0[INeigh];
		Physics->n   [I] = Physics->n   [INeigh];
		Physics->rho[I] = Physics->rho[INeigh];
		Physics->k  [I] = Physics->k  [INeigh];
		Physics->G  [I] = Physics->G  [INeigh];
		Physics->sigma_xx_0 [I] = Physics->sigma_xx_0 [INeigh];
		Physics->cohesion     [I] = Physics->cohesion     [INeigh];
		Physics->frictionAngle[I] = Physics->frictionAngle[INeigh];


		if (BCThermal->type[IBC]==DirichletGhost) { // Dirichlet
			Physics->T[I] = 2.0*BCThermal->value[IBC] - Physics->T[INeigh];
		}
		else if (BCThermal->type[IBC]==NeumannGhost) { // Neumann
			if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
				Physics->T[I] = Physics->T[INeigh] - BCThermal->value[IBC]*Grid->dx;
			}
			if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
				Physics->T[I] = Physics->T[INeigh] + BCThermal->value[IBC]*Grid->dy;
			}
		}
	}

	printf("end neighbour stuff");

	/*
	printf("=== Check sigma_xx 0b ===\n");
	C = 0;
	//int ix, iy;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.3f  ", Physics->sigma_xx_0[C]);
			C++;

		}
		printf("\n");
	}
*/

	// ==================================
	// Interpolate to nodes
	// ==================================

	// Reinitialize sum of weights
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		sumOfWeights[i] = 0;
	}

	int signX, signY, iNodeNeigh;
	xMod[0] =  1; yMod[0] =  1;
	xMod[1] =  0; yMod[1] =  1;
	xMod[2] =  1; yMod[2] =  0;
	xMod[3] =  0; yMod[3] =  0;
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight, signX, signY, iNodeNeigh) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) { // Gives better result not to give contribution from the boundaries
		for (ix = 0; ix < Grid->nxS; ++ix) { // I don't get why though
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {
				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				phase = thisParticle->phase;

				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}



				locX = fabs(locX);
				locY = fabs(locY);

				for (i=0; i<4; i++) {
					iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

					if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
						printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
						printf("IX = %i, IY = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY);
						exit(0);
					}


					//printf("iNodeNeigh = %i, signX = %i, signY = %i\n", iNodeNeigh, signX, signY);
					weight = (locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5);

					sigma_xy_0 	[iNodeNeigh*4+i] += thisParticle->sigma_xy_0 * weight;
					sumOfWeights[iNodeNeigh*4+i] += weight; // using the same arrays



				}

				thisParticle = thisParticle->next;
			}
		}
	}

	printf("end first loop for sigma_xy\n");

	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		I = 4*iNode;
		sum = sumOfWeights[I+0] + sumOfWeights[I+1] + sumOfWeights[I+2] + sumOfWeights[I+3];
		//printf("%.2f %.2f %.2f %.2f\n", sumOfWeights[I+0], sumOfWeights[I+1], sumOfWeights[I+2], sumOfWeights[I+3]);
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		Physics->sigma_xy_0 [iNode] = ( sigma_xy_0[I+0] +  sigma_xy_0[I+1] +  sigma_xy_0[I+2] +  sigma_xy_0[I+3]) / sum ; // harmonic average


	}

	printf("end filling loop for sigma_xy");









	if (DEBUG) {

		printf("=== Check eta0 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta0[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check eta 1 ===\n");
				C = 0;
				//int ix, iy;
				for (iy = 0; iy < Grid->nyEC; ++iy) {
					for (ix = 0; ix < Grid->nxEC; ++ix) {
						printf("%.3f  ", Physics->eta[C]);
						C++;
					}
					printf("\n");
				}
				/*
		printf("=== Check rho 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->rho[C]);
				C++;
			}
			printf("\n");
		}
		printf("=== Check k 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->k[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== Check T 1 ===\n");
				C = 0;
				//int ix, iy;
				for (iy = 0; iy < Grid->nyEC; ++iy) {
					for (ix = 0; ix < Grid->nxEC; ++ix) {
						printf("%.3f  ", Physics->T[C]);
						C++;
						if (isnan((float) Physics->T[C])!=0) {
							exit(0);
						}
					}
					printf("\n");
				}
		printf("=== Check G 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->G[C]);
				C++;
			}
			printf("\n");
		}
		 */

		/*
		printf("=== Check sigma_xx 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->sigma_xx_0[C]);
				C++;
				if (isnan((float) Physics->sigma_xx_0[C])!=0) {
					//exit(0);
				}
			}
			printf("\n");
		}
		printf("=== Check sigma_xy 1 ===\n");
		C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.3e  ", Physics->sigma_xy_0[C]);
				C++;
			}
			printf("\n");
		}
		*/
	}

	free(sumOfWeights);
	free(eta0);
	free(n);
	free(rho);
	free(k);
	free(T);
	free(G);
	free(sigma_xx_0);
	free(sigma_xy_0);
	free(StrainRateInvariant);
	free(cohesion);
	free(frictionAngle);






}



















void Physics_interpTempFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	int signX, signY;
	compute dum;
	int i;
	/*
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->Dsigma_xx_0[i] = 1;
		Physics->DT[i] = 1;
	}
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->Dsigma_xy_0[i] = -1;
	}
	*/
	/*
	FOR_PARTICLES
		thisParticle->sigma_xx_0 = 1.0;
	END_PARTICLES
	*/


	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, signX, signY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;


				thisParticle->T  += ( .25*(1.0-locX)*(1.0-locY)*Physics->DT[ix  +(iy  )*Grid->nxEC]
								    + .25*(1.0-locX)*(1.0+locY)*Physics->DT[ix  +(iy+1)*Grid->nxEC]
								    + .25*(1.0+locX)*(1.0+locY)*Physics->DT[ix+1+(iy+1)*Grid->nxEC]
								    + .25*(1.0+locX)*(1.0-locY)*Physics->DT[ix+1+(iy  )*Grid->nxEC] );


				thisParticle = thisParticle->next;
			}
		}
	}


}











void Physics_interpStressesFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal)
{

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	int signX, signY;
	compute dum;
	int i;
	/*
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->Dsigma_xx_0[i] = 1;
		Physics->DT[i] = 1;
	}
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->Dsigma_xy_0[i] = -1;
	}
	*/
	/*
	FOR_PARTICLES
		thisParticle->sigma_xx_0 = 1.0;
	END_PARTICLES
	*/


	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, signX, signY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
				locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;




				thisParticle->sigma_xx_0  += ( .25*(1.0-locX)*(1.0-locY)*Physics->Dsigma_xx_0[ix  +(iy  )*Grid->nxEC]
										     + .25*(1.0-locX)*(1.0+locY)*Physics->Dsigma_xx_0[ix  +(iy+1)*Grid->nxEC]
											 + .25*(1.0+locX)*(1.0+locY)*Physics->Dsigma_xx_0[ix+1+(iy+1)*Grid->nxEC]
											 + .25*(1.0+locX)*(1.0-locY)*Physics->Dsigma_xx_0[ix+1+(iy  )*Grid->nxEC] );


				// Sigma_xy is stored on the node, therefore there are 4 possible squares to interpolate from
				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}


				locX = fabs(locX)*2-1;
				locY = fabs(locY)*2-1;


				thisParticle->sigma_xy_0  += ( .25*(1.0-locX)*(1.0-locY)*Physics->Dsigma_xy_0[ix      +(iy  )    *Grid->nxS]
										     + .25*(1.0-locX)*(1.0+locY)*Physics->Dsigma_xy_0[ix      +(iy+signY)*Grid->nxS]
										     + .25*(1.0+locX)*(1.0+locY)*Physics->Dsigma_xy_0[ix+signX+(iy+signY)*Grid->nxS]
										     + .25*(1.0+locX)*(1.0-locY)*Physics->Dsigma_xy_0[ix+signX+(iy  )    *Grid->nxS] );

			 	//printf("thisParticle->sigma_xx_0 = %.3f, thisParticle->sigma_xy_0 = %.3f, ix = %i, iy = %i,signX = %i, signY = %i\n", thisParticle->sigma_xx_0, thisParticle->sigma_xy_0, ix, iy , signX, signY);

				//printf("ix = %i, iy = %i, locX = %.3f, locY = %.3f, T[0] = %.3f, T[1]=%.3f, T[2]=%.3f, T[3]=%.3f, Tpart= %.3f\n",ix, iy, locX, locY, Physics->T[ix  +(iy  )*Grid->nxEC], Physics->T[ix  +(iy+1)*Grid->nxEC], Physics->T[ix+1+(iy+1)*Grid->nxEC], Physics->T[ix+1+(iy)*Grid->nxEC], thisParticle->T);

				thisParticle = thisParticle->next;
			}
		}
	}


}























void Physics_interpFromCellToNode(Grid* Grid, compute* CellValue, compute* NodeValue)
{
	// UC is a scalar CellValue defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	int iNW, iNE, iSW, iSE;



	// CellValue interpolated on the center nodes
	// ======================================
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			I = ix + iy*Grid->nxS;
			iNW = (ix)+ (iy+1)   *Grid->nxEC;
			iNE = ix+1    + (iy+1)   *Grid->nxEC;
			iSW = (ix)+(iy)*Grid->nxEC;
			iSE = ix+1    +(iy)*Grid->nxEC;
			NodeValue[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}


	if (DEBUG) {
		int C = 0;
		printf("=== Check CellValue ===\n");
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", CellValue[C]);
				C++;
			}
			printf("\n");
		}


		C = 0;
		printf("=== Check NodeValue ===\n");
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.3e  ", NodeValue[C]);
				C++;
			}
			printf("\n");
		}

	}

}


void Physics_set_VxVyP_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int InoDir, INeigh;
	compute maxVx = 0;
	compute maxVy = 0;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = -1;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = -1;
	}
	for (i = 0; i < Grid->nCTot; ++i) {
		Physics->P[i] = -1;
	}


	// Set Vx
	// =========================
	C = 0;
	int IBC;
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			I = ix + iy*Grid->nxVx;
			InoDir = Numbering->map[I];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vx[C] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1; // BC nodes are numbered -1 to -n
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vx[C] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (iy==0)  // lower boundary
						INeigh = Numbering->map[  ix + (iy+1)*Grid->nxVx  ];
					if (iy==Grid->nyVx-1)  // lower boundary
						INeigh = Numbering->map[  ix + (iy-1)*Grid->nxVx  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vx[C] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (iy==0)  // lower boundary
							Physics->Vx[C] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dy;
						if (iy==Grid->nyVx-1)  // lower boundary
							Physics->Vx[C] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dy;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}

			// Get maxVx
			if (Physics->Vx[C]*Physics->Vx[C] > maxVx)
				maxVx = Physics->Vx[C]*Physics->Vx[C];


			C++;
		}
	}

	// Set Vy
	// =========================

	C = 0;
	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			I = ix + iy*Grid->nxVy + Grid->nVxTot;

			InoDir = Numbering->map[I];

			if (InoDir>=0) { // Not a Dirichlet node
				Physics->Vy[C] = EqSystem->x[InoDir];
			}
			// Deal with boundary conditions
			else  { // Dirichlet or Neumann
				IBC = abs(InoDir)-1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					Physics->Vy[C] = BC->value[IBC];
				}
				else { // on a ghost node

					// Get neighbours index
					if (ix==0)  // lower boundary
						INeigh = Numbering->map[  ix+1 + (iy)*Grid->nxVy + Grid->nVxTot ];
					if (ix==Grid->nxVy-1)  // lower boundary
						INeigh = Numbering->map[  ix-1 + (iy)*Grid->nxVy + Grid->nVxTot  ];


					if (BC->type[IBC]==DirichletGhost) { // Dirichlet
						Physics->Vy[C] = 2.0*BC->value[IBC] - EqSystem->x[INeigh];
					}
					else if (BC->type[IBC]==NeumannGhost) { // Neumann
						if (ix==0)  // lower boundary
							Physics->Vy[C] = EqSystem->x[INeigh] - BC->value[IBC]*Grid->dx;
						if (ix==Grid->nxVy-1)  // lower boundary
							Physics->Vy[C] = EqSystem->x[INeigh] + BC->value[IBC]*Grid->dx;
					}
					else {
						printf("error: unknown boundary type\n");
						exit(0);
					}
				}
			}

			if (Physics->Vy[C]*Physics->Vy[C] > maxVy)
				maxVy = Physics->Vy[C]*Physics->Vy[C];
			C++;
		}
	}


	int IE; // Index of embedded nodes
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			IE = ix + iy*Grid->nxEC; // Index embedded
			I = (ix-1) + (iy-1)*Grid->nxC + Grid->nVxTot + Grid->nVyTot; // Index in NumStokes
			InoDir = Numbering->map[I];
			if (InoDir>=0) { // Not a Dirichlet node
				Physics->P[IE] = EqSystem->x[InoDir];
			} else {
				IBC = abs(InoDir)-1;
				Physics->P[IE] = BC->value[ IBC ];
			}
		}
	}
	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->P[ix + iy*Grid->nxEC] = Physics->P[ix-1 + (iy)*Grid->nxEC];
	}





	Physics->maxV = sqrt(maxVx+maxVy);





	if (DEBUG) {
		// Check Vx
		// =========================
		printf("=== Vx ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyVx; ++iy) {
			for (ix = 0; ix < Grid->nxVx; ++ix) {
				printf("%.3e  ", Physics->Vx[C]);
				C++;
			}
			printf("\n");
		}


		// Check Vy
		// =========================
		C = 0;
		printf("=== Vy ===\n");
		for (iy = 0; iy < Grid->nyVy; ++iy) {
			for (ix = 0; ix < Grid->nxVy; ++ix) {
				printf("%.3e  ", Physics->Vy[C]);
				C++;
			}
			printf("\n");
		}

		// Check P
		// =========================
		printf("=== P ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3e  ", Physics->P[C]);
				C++;
			}
			printf("\n");
		}
	}






}
















void Physics_set_T_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumThermal, EqSystem* EqThermal)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int INeigh, IBC;
	//int InoDir, INeigh;
	//compute maxVx = 0;
	//compute maxVy = 0;
	// Init Vx, Vy, P to -1, for debugging purposes
	// =========================

	/*
	for (i = 0; i < Grid->nCTot; ++i) {
		Physics->T[i] = EqThermal->x[i];
	}
	*/







	printf("== Check T filling ==\n");
	C = 0;
	i = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		for (ix = 0; ix<Grid->nxEC; ix++) {
			I = NumThermal->map[C];
			if (I>=0) {
				Physics->DT[C] = EqThermal->x[I] - Physics->T[C];
				Physics->T[C] = EqThermal->x[I];
			}
			else {

				IBC = abs(I)-1; // BC nodes are numbered -1 to -n

				// Get neighbours index
				if (iy==0) { // lower boundary
					if (ix==0) {
						INeigh = NumThermal->map[  ix+1 + (iy+1)*Grid->nxEC  ];
					} else if (ix==Grid->nxEC-1) {
						INeigh = NumThermal->map[  ix-1 + (iy+1)*Grid->nxEC  ];
					} else {
						INeigh = NumThermal->map[  ix + (iy+1)*Grid->nxEC  ];
					}

				} else if (iy==Grid->nyEC-1)  { //  upper boundary
					if (ix==0) {
						INeigh = NumThermal->map[  ix+1 + (iy-1)*Grid->nxEC  ];
					} else if (ix==Grid->nxEC-1) {
						INeigh = NumThermal->map[  ix-1 + (iy-1)*Grid->nxEC  ];
					} else {
						INeigh = NumThermal->map[  ix + (iy-1)*Grid->nxEC  ];
					}
				} else if (ix==0) { // left boundary
					INeigh = NumThermal->map[  ix+1 + (iy)*Grid->nxEC  ];
				} else if (ix==Grid->nxEC-1) { // right boundary
					INeigh = NumThermal->map[  ix-1 + (iy)*Grid->nxEC  ];
				}





				if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Physics->DT[C] = 2.0*BC->value[IBC] - EqThermal->x[INeigh] - Physics->T[C];
					Physics->T[C] = 2.0*BC->value[IBC] - EqThermal->x[INeigh];
				}
				else if (BC->type[IBC]==NeumannGhost) { // Neumann
					if (ix==0 || ix==Grid->nxEC-1)  {// left or right boundary
						Physics->DT[C] = EqThermal->x[INeigh] - BC->value[IBC]*Grid->dx - Physics->T[C];
						Physics->T[C] = EqThermal->x[INeigh] - BC->value[IBC]*Grid->dx;
					}
					if (iy==0 || iy==Grid->nyEC-1) { // top or bottom boundary
						Physics->DT[C] = EqThermal->x[INeigh] + BC->value[IBC]*Grid->dy - Physics->T[C];
						Physics->T[C] = EqThermal->x[INeigh] + BC->value[IBC]*Grid->dy;
					}
				}
				else {
					printf("error: unknown boundary type\n");
					exit(0);
				}
				//printf("C=%i, IBC=%i, Type=%i, value=%.3f, valueNeigh=%.3f, FinalValue=%.3f\n",C, IBC,BC->type[IBC], BC->value[IBC], EqThermal->x[INeigh], Physics->T[C]);


				//Physics->T[C] = BC->value[abs(I)];
			}
			C++;
		}
	}



	if (DEBUG) {
		// Check T
		// =========================
		printf("=== T ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->T[C]);
				C++;
			}
			printf("\n");
		}
	}




}






void Physics_computeStressChanges(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes)
{

	// see Taras' book p. 186
	int ix, iy, iCell, iNode;
	compute Z;
	compute Eps_xx, Eps_xy, Eps_yy;
	compute dVxdy, dVydx;
	compute GShear;
	// compute stress
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix + iy*Grid->nxEC;
			Eps_xx 	=  (Physics->Vx[ix + iy*Grid->nxVx] - Physics->Vx[ix-1 + iy*Grid->nxVx])/Grid->dx;


			//printf("Eps_xx = %.3e, Eps_yy = % .3e, Eps_ref = %.3e\n", Eps_xx, Eps_yy, Physics->epsRef);

			//Eps_xx 	=  (Physics->Vx[ix + iy*Grid->nxVx] - Physics->Vx[ix-1 + iy*Grid->nxVx])/Grid->dx;

			Z 		= (Physics->G[iCell]*Physics->dt)  /  (Physics->eta[iCell] + Physics->G[iCell]*Physics->dt);

			Physics->Dsigma_xx_0[iCell] = ( 2*Physics->eta[iCell] * Eps_xx  -  Physics->sigma_xx_0[iCell] ) * Z;
			//printf("eta = %.2e, G = %.2e, Eps_xx = %.2e, Z = %.2e, sigma_xx_0 = %.2e, dt = %.2e\n", Physics->eta[iCell], Physics->G[iCell], Eps_xx, Z, Physics->sigma_xx_0[iCell], Physics->dt);

		}
	}



	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->Dsigma_xx_0[ix + iy*Grid->nxEC] = Physics->Dsigma_xx_0[ix-1 + (iy)*Grid->nxEC];
	}






	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]
					- Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]
				    - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Eps_xy = 0.5*(dVxdy+dVydx);

			//printf("iNode = %i, dVxdy = %.2e, dVydx = %.2e, VxS = %.3e, VxN = %.3e, VyW = %.3e, VyE = %.3e\n", iNode, dVxdy, dVydx, Physics->Vx[ix  + (iy  )*Grid->nxVx] , Physics->Vx[ix  + (iy+1)*Grid->nxVx], Physics->Vy[ix  + iy*Grid->nxVy] , Physics->Vy[ix+1+ iy*Grid->nxVy] );

			//GShear = 0.25 * ( Physics->G[ix+iy*Grid->nxEC] + Physics->G[ix+1+iy*Grid->nxEC] + Physics->G[ix+(iy+1)*Grid->nxEC] + Physics->G[ix+1+(iy+1)*Grid->nxEC] );

			GShear = Physics->GShear[iNode];

			Z 		= (GShear*Physics->dt)  /  (Physics->etaShear[iNode] + GShear*Physics->dt);

			Physics->Dsigma_xy_0[iNode] = ( 2*Physics->etaShear[iNode] * Eps_xy   -   Physics->sigma_xy_0[iNode] ) * Z;
			//printf("Dsigma_xy_0[%i] = %.3e, ix = %i, iy = %i\n", iNode, Physics->Dsigma_xy_0[iNode], ix, iy);
			//printf("eta = %.2e, G = %.2e, Eps_xy = %.2e, Z = %.2e, sigma_xy_0 = %.2e\n", Physics->etaShear[iNode], Physics->GShear[iNode], Eps_xy, Physics->dt, Physics->sigma_xy_0[iNode]);
		}
	}





/*
	// Check T
	// =========================
	printf("=== DSigma_xx_0 ===\n");
	int C = 0;
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			printf("%.3e  ", Physics->Dsigma_xx_0[C]);
			C++;
		}
		printf("\n");
	}
	printf("=== DSigma_xy_0 ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			printf("%.3e  ", Physics->Dsigma_xy_0[C]);
			C++;
		}
		printf("\n");
	}
*/


}








void Physics_computeStrainRateInvariant(Physics* Physics, Grid* Grid, compute* StrainRateInvariant)
{
	// Definition of second invariant: // E_II = sqrt( Eps_xx^2 + Eps_xy^2  );
	// Declarations
	// =========================
	int ix, iy, I, iNode, Ix, Iy, IE;
	compute dVxdy, dVydx, dVxdx, dVydy;
	// ix, iy modifiers
	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			IE = ix+iy*Grid->nxEC;
			//I = (ix-1)+(iy-1)*Grid->nxC;

			// Compute Eps_xy at the four nodes of the cell
			// 1. Sum contributions
			dVxdy = 0;
			dVydx = 0;
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];
				dVxdy += ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
					     - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;

				dVydx += ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
						 - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;
			}
			// 2. Average
			dVxdy /= 4;
			dVydx /= 4;

			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
								 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
								 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;


			StrainRateInvariant[IE] = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +  0.5*dVxdx*dVxdx  +  0.5*dVydy*dVydy);
			if (StrainRateInvariant[IE]<1E-12) { // tolerance for accuracy of 0
				StrainRateInvariant[IE] = Physics->epsRef;
			}

		}
	}

	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		StrainRateInvariant[ix + iy*Grid->nxEC] = StrainRateInvariant[ix-1 + (iy)*Grid->nxEC];
	}
}


/*
It's actually not the vorticity but the rotation rate
void Physics_computeVorticity(Physics* Physics, Grid* Grid, compute* Vorticity) {
	// The vorticity array must have size nxEC*nyEX

	// In 2D the vorticity is defined as: omega = dVy/dx - dVx/dy

	int iy, ix, iCell;



	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell = ix + iy*Grid->nxEC;
			Vorticity[iCell]  = (Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix+(iy-1)*Grid->nxVy])/Grid->dx;
			Vorticity[iCell] += (Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+(iy)*Grid->nxVy])/Grid->dy;
		}
	}

	// Loop through the sides and put the same values
	// Replace boundary values by their neighbours
		// lower boundary
		iy = 0;
		for (ix = 0; ix<Grid->nxEC; ix++) {
			Vorticity[ix + iy*Grid->nxEC] = Vorticity[ix + (iy+1)*Grid->nxEC];
		}
		// upper boundary
		iy = Grid->nyEC-1;
		for (ix = 0; ix<Grid->nxEC; ix++) {
			Vorticity[ix + iy*Grid->nxEC] = Vorticity[ix + (iy-1)*Grid->nxEC];
		}
		// left boundary
		ix = 0;
		for (iy = 0; iy<Grid->nyEC; iy++) {
			Vorticity[ix + iy*Grid->nxEC] = Vorticity[ix+1 + (iy)*Grid->nxEC];
		}
		// right boundary
		ix = Grid->nxEC-1;
		for (iy = 0; iy<Grid->nyEC; iy++) {
			Vorticity[ix + iy*Grid->nxEC] = Vorticity[ix-1 + (iy)*Grid->nxEC];
		}



}

*/

void Physics_computeEta(Physics* Physics, Grid* Grid)
{
	int iCell, C, iy, ix, i;
	compute eta0;
	compute sigma_y, sigmaII, eta_corr;
	compute* EII = (compute*) malloc(Grid->nECTot*sizeof(compute));
	int nLocIt = 10000;
	int it = 0;
	compute eta_y;
	compute EII_visc;
	Physics_computeStrainRateInvariant(Physics, Grid, EII);
	compute eta;
	compute tolerance;
	compute EIILoc;
	if (Physics->itNonLin<1) {
		tolerance = 1E-6;
	} else {
		tolerance = 1E-6;
	}

	if ( Physics->time < 0  && Physics->itNonLin==0) {
		for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
			Physics->eta[iCell] = Physics->eta0[iCell];
		}

		if (Physics->eta[iCell]<Physics->etaMin) {
			Physics->eta[iCell] = Physics->etaMin;
		}
		else if (Physics->eta[iCell]>Physics->etaMax) {
			Physics->eta[iCell] = Physics->etaMax;
		}

	}

	else {
#pragma omp parallel for private(iCell, sigma_y, sigmaII) schedule(static,32)
		for (iCell = 0; iCell < Grid->nECTot; ++iCell) {

			// Compute powerlaw rheology
			Physics->eta[iCell] = Physics->eta0[iCell] * pow(EII[iCell]/Physics->epsRef     ,    1.0/Physics->n[iCell] - 1.0);


			// Compute the yield stress
			sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Physics->P[iCell] * sin(Physics->frictionAngle[iCell]);

			sigmaII = 2*Physics->eta[iCell] * EII[iCell];



			if (sigmaII>sigma_y) {

				Physics->eta[iCell] = (sigma_y /(2*EII[iCell]));

			}





			if (Physics->eta[iCell]<Physics->etaMin) {
				Physics->eta[iCell] = Physics->etaMin;
			}
			else if (Physics->eta[iCell]>Physics->etaMax) {
				Physics->eta[iCell] = Physics->etaMax;
			}

		}
	}


















	if (DEBUG) {
		printf("===== Compute Eta =====\n");
		printf("=== Check eta0  ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta0[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== Check EII ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.2e  ", EII[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== Check eta ===\n");
		C = 0;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
		}


		printf("=== Check n ===\n");
			C = 0;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					printf("%.3f  ", Physics->n[C]);
					C++;
				}
				printf("\n");
			}
	}


	free(EII);

}




