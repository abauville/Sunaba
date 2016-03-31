/*
 * Physics.c
 *
 *  Created on: Feb 24, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal)
{

	// Declarations
	// =========================
	int iCell, i;
	int nNeighbours = 4;
	coord locX, locY;

	coord dx = Grid->dx;
	coord dy = Grid->dy;
	compute* sumOfWeights 	= (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	compute* eta = (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* rho = (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute*   k = (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));
	compute* T   = (compute*) malloc(nNeighbours * Grid->nECTot * sizeof(compute));

	// Reinitialize Physics array
	for (i = 0; i < nNeighbours * Grid->nECTot; ++i) {
		eta[i] = 0;
		rho[i] = 0;
		T  [i] = 0;
		k  [i] = 0;
		sumOfWeights[i] = 0;
	}

	compute Tsum;
	compute Tweight;
	compute weight;
	compute locEta, locEps_II, dVxdx, dVxdy, dVydx, dVydy;

	int iNode, IX, IY;

	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};

	int phase;

	int nxEC = Grid->nxEC;
	int xMod[4], yMod[4], Ix[4], ix, Iy[4], iy;

	int I;

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


	int quadrant = 0;
	int oldQuadrant = 0;

	SingleParticle* thisParticle = NULL;
	double tocTot = 0;
	// Loop through inner cells
	// ========================
	iNode = 0;
	//int Count = 0;
//#pragma omp parallel for private(ix, iNode, thisParticle, locX, locY, phase, Ix, Iy, i, locEta, iNode, dVxdy, dVydx, dVxdx, dVydy, locEps_II) schedule(static,32)







	compute* StrainRateInvariant = (compute*) malloc(Grid->nECTot * sizeof(compute));



	Physics_computeStrainRateInvariant(Physics, Grid, StrainRateInvariant);









	//printf("=== Part Temp ===\n");
	// Loop through inner nodes
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
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

					//Tweight = weight;

					// code to avoid giving contribution to boundary cells in an asymmetric way
					// hopefully can be optimized


					if (i<2 && iy+IyN[i]==1 && locY>0) {
						weight = 0;
					}
					else if (i>1 && iy+IyN[i]==Grid->nyEC-2 && locY<0) {
						weight = 0;
					}
					else if ((i==0 || i==2) && ix+IxN[i]==1 && locX>0) {
						weight = 0;
					}
					else if ((i==1 || i==3) && ix+IxN[i]==Grid->nxEC-2 && locX<0) {
						weight = 0;
					}



					eta			[iCell*4+i] += MatProps->eta0[phase] * pow(  (StrainRateInvariant[iCell]/Physics->epsRef)    ,   1.0/MatProps->n[phase]-1.0 ) * weight;
					rho			[iCell*4+i] += MatProps->rho0[phase] * (1+MatProps->beta[phase]*Physics->P[iCell]) * (1-MatProps->alpha[phase]*Physics->T[iCell])   *  weight;
					k			[iCell*4+i] += MatProps->k   [phase] * weight;
					T 			[iCell*4+i] += thisParticle->T * weight;
					sumOfWeights[iCell*4+i] += weight;
				}
				thisParticle = thisParticle->next;
			}
			//printf("\n");
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

			Physics->T  [iCell] = (  T[I+0] +   T[I+1] +   T[I+2] +   T[I+3]) / sum;
			Physics->eta[iCell] = (eta[I+0] + eta[I+1] + eta[I+2] + eta[I+3]) / sum;
			Physics->rho[iCell] = (rho[I+0] + rho[I+1] + rho[I+2] + rho[I+3]) / sum;
			Physics->k  [iCell] = (  k[I+0] +   k[I+1] +   k[I+2] +   k[I+3]) / sum;

			//printf("Physics->T[%i] %.2e, T = %.2f %.2f %.2f %.2f, sum = %.2f\n", iCell, Physics->T[iCell], T[I+0], T[I+1], T[I+2], T[I+3], sum);

			if (Physics->eta[iCell]<Physics->etaMin) {
				Physics->eta[iCell] = Physics->etaMin;
			}
			else if (Physics->eta[iCell]>Physics->etaMax) {
				Physics->eta[iCell] = Physics->etaMax;
			}
		}
	}


	// Replace boundary values by their neighbours
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->eta[ix + iy*Grid->nxEC] = Physics->eta[ix + (iy+1)*Grid->nxEC];
		Physics->rho[ix + iy*Grid->nxEC] = Physics->rho[ix + (iy+1)*Grid->nxEC];
		Physics->k  [ix + iy*Grid->nxEC] = Physics->k  [ix + (iy+1)*Grid->nxEC];
	}
	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		Physics->eta[ix + iy*Grid->nxEC] = Physics->eta[ix + (iy-1)*Grid->nxEC];
		Physics->rho[ix + iy*Grid->nxEC] = Physics->rho[ix + (iy-1)*Grid->nxEC];
		Physics->k  [ix + iy*Grid->nxEC] = Physics->k  [ix + (iy-1)*Grid->nxEC];
	}
	// left boundary
	ix = 0;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->eta[ix + iy*Grid->nxEC] = Physics->eta[ix+1 + (iy)*Grid->nxEC];
		Physics->rho[ix + iy*Grid->nxEC] = Physics->rho[ix+1 + (iy)*Grid->nxEC];
		Physics->k  [ix + iy*Grid->nxEC] = Physics->k  [ix+1 + (iy)*Grid->nxEC];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 0; iy<Grid->nyEC; iy++) {
		Physics->eta[ix + iy*Grid->nxEC] = Physics->eta[ix-1 + (iy)*Grid->nxEC];
		Physics->rho[ix + iy*Grid->nxEC] = Physics->rho[ix-1 + (iy)*Grid->nxEC];
		Physics->k  [ix + iy*Grid->nxEC] = Physics->k  [ix-1 + (iy)*Grid->nxEC];
	}







	if (DEBUG) {
		printf("=== Check eta 1 ===\n");
		int C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (ix = 0; ix < Grid->nxEC; ++ix) {
				printf("%.3f  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
		}
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
					}
					printf("\n");
				}
	}

	free(sumOfWeights);
	free(eta);
	free(rho);
	free(k);
	free(T);
	free(StrainRateInvariant);




	/*

		// Loop through the inner bottom and upper boundaries
		iy = 0;
		int C = 2;
		int iOut;
		for (iOut = 0; iOut<2; ++iOut) {
			for (ix = 1; ix < Grid->nxS-1; ++ix) {
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle!=NULL) {
					locX = (thisParticle->x-Grid->xmin)/dx - ix;
					locY = (thisParticle->y-Grid->ymin)/dy - iy;
					phase = thisParticle->phase;

					for (i=0;i<2;++i) {
						I = (ix+IxN[C+i] + (iy+IyN[C+i]) * nxC) * nNeighbours;

						weight = fabs((locX + xMod[C+i]*0.5)   *   (locY + yMod[C+i]*0.5));

						eta			[I+C+i] += MatProps->eta0[phase] * weight;
						rho			[I+C+i] += MatProps->rho0[phase] * weight;
						T 			[I+C+i] += thisParticle->T * weight;
						sumOfWeights[I+C+i] += weight;
					}
					thisParticle = thisParticle->next;
				}
			}
			C = 0;
			iy = Grid->nyS-1;
		}








		// Loop through the inner left and right boundaries
		ix = 0;
		C  = 1;
		int Cperiodic = 0;
		int periodicMod = 0;
		if (BC->SetupType==SimpleShearPeriodic) {
			periodicMod = Grid->nxC;
		}
		else {
			periodicMod = 0;
		}
		for (iOut = 0; iOut<2; ++iOut) {
			for (iy = 1; iy < Grid->nyS-1; ++iy) {
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle!=NULL) {
					locX = (thisParticle->x-Grid->xmin)/dx - ix;
					locY = (thisParticle->y-Grid->ymin)/dy - iy;
					phase = thisParticle->phase;
					for (i=0;i<3;i+=2) {
						I = ix+IxN[C+i] + (iy+IyN[C+i]) * nxC * nNeighbours;
						weight = fabs((locX + xMod[C+i]*0.5)   *   (locY + yMod[C+i]*0.5));
						eta			[I+C+i] += MatProps->eta0[phase] * weight;
						rho			[I+C+i] += MatProps->rho0[phase] * weight;
						T 			[I+C+i] += thisParticle->T * weight;
						sumOfWeights[I+C+i] += weight;
					}


					if (BC->SetupType==SimpleShearPeriodic) {
							for (i=0;i<3;i+=2) {
								I = (ix+IxN[Cperiodic+i] + (iy+IyN[Cperiodic+i]) * nxC + periodicMod) * nNeighbours ;
								weight = fabs((locX + xMod[C+i]*0.5)   *   (locY + yMod[C+i]*0.5));
								eta			[I+Cperiodic+i] += MatProps->eta0[phase] * weight;
								rho			[I+Cperiodic+i] += MatProps->rho0[phase] * weight;
								T 			[I+Cperiodic+i] += thisParticle->T * weight;
								sumOfWeights[I+Cperiodic+i] += weight;
							}
						}

					thisParticle = thisParticle->next;
				}




			}
			C = 0;
			Cperiodic = 1;
			ix = Grid->nxS-1;
			if (BC->SetupType==SimpleShearPeriodic) {
				periodicMod = - Grid->nxC;
			}
		}


		// Corners
		C = 3;
		for (iy = 0; iy < Grid->nyS; iy+=Grid->nyS-1) {
			for (ix = 0; ix < Grid->nxS; ix+=Grid->nxS-1) {
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle!=NULL) {
					locX = (thisParticle->x-Grid->xmin)/dx - ix;
					locY = (thisParticle->y-Grid->ymin)/dy - iy;
					phase = thisParticle->phase;

					I = (ix+IxN[C] + (iy+IyN[C]) * nxC) * nNeighbours;
					weight = fabs((locX + xMod[C]*0.5)   *   (locY + yMod[C]*0.5));
					eta			[I+C] += MatProps->eta0[phase] * weight;
					rho			[I+C] += MatProps->rho0[phase] * weight;
					T 			[I+C] += thisParticle->T * weight;
					sumOfWeights[I+C] += weight;
					thisParticle = thisParticle->next;
				}
				C--;
			}
		}

	*/





}




void Physics_interpFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal)
{



	int ix, iy, iNode;
	SingleParticle* thisParticle;
	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;


	// Loop through nodes
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

				//printf("ix = %i, iy = %i, locX = %.3f, locY = %.3f, T[0] = %.3f, T[1]=%.3f, T[2]=%.3f, T[3]=%.3f, Tpart= %.3f\n",ix, iy, locX, locY, Physics->T[ix  +(iy  )*Grid->nxEC], Physics->T[ix  +(iy+1)*Grid->nxEC], Physics->T[ix+1+(iy+1)*Grid->nxEC], Physics->T[ix+1+(iy)*Grid->nxEC], thisParticle->T);

				thisParticle = thisParticle->next;
			}
		}
	}




	/*
	// Loop through inner nodes of the bottom boundary
	int INeigh, IBC;
	iy = 0;
	for (ix = 1; ix < Grid->nxS-1; ix++) {
		iNode = ix  + (iy  )*Grid->nxS;
		thisParticle = Particles->linkHead[iNode];

		// Loop through the particles in the shifted cell
		// ======================================
		while (thisParticle!=NULL) {
			locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
			locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;


			// Get neighbours index




			thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*BCThermal->value[(int) abs(NumThermal->map[ix-1+1+(iy-1+1 )*(Grid->nxC+2)]) ] // the +1 is because the Numbering map contains the ghost nodes but Physics->T does not
							   + .25*(1.0-locX)*(1.0+locY)*Physics->T[ix-1+(iy)*Grid->nxC]
							   + .25*(1.0+locX)*(1.0+locY)*Physics->T[ix  +(iy)*Grid->nxC]
							   + .25*(1.0+locX)*(1.0-locY)*BCThermal->value[(int) abs(NumThermal->map[ix+1  +(iy-1+1 )*(Grid->nxC+2)]) ];

			thisParticle = thisParticle->next;
		}
	}





	// Loop through inner nodes of the top boundary
	iy = Grid->nyS-1;
	for (ix = 1; ix < Grid->nxS-1; ix++) {
		iNode = ix  + (iy  )*Grid->nxS;
		thisParticle = Particles->linkHead[iNode];



		// Loop through the particles in the shifted cell
		// ======================================
		while (thisParticle!=NULL) {
			locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
			locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;



			thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*Physics->T[ix-1+(iy-1 )*Grid->nxC]
							   + .25*(1.0-locX)*(1.0+locY)*BCThermal->value[(int) abs(NumThermal->map[ix-1+1+(iy+1 )*(Grid->nxC+2)]) ]
							   + .25*(1.0+locX)*(1.0+locY)*BCThermal->value[(int) abs(NumThermal->map[ix+1  +(iy+1 )*(Grid->nxC+2)]) ]
							   + .25*(1.0+locX)*(1.0-locY)*Physics->T[ix  +(iy-1)*Grid->nxC];

			thisParticle = thisParticle->next;
		}
	}

	// Loop through inner nodes of the left
	ix = 0;
	for (iy = 1; iy < Grid->nyS-1; iy++) {
		iNode = ix  + (iy  )*Grid->nxS;
		thisParticle = Particles->linkHead[iNode];
		// Loop through the particles in the shifted cell
		// ======================================
		while (thisParticle!=NULL) {
			locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
			locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;

			if (BCStokes->SetupType==SimpleShearPeriodic) {
				thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*Physics->T[ix-1+(iy-1)*Grid->nxC]
								   + .25*(1.0-locX)*(1.0+locY)*Physics->T[ix-1+(iy)*Grid->nxC]
								   + .25*(1.0+locX)*(1.0+locY)*Physics->T[ix  +(iy)*Grid->nxC]
								   + .25*(1.0+locX)*(1.0-locY)*Physics->T[ix  +(iy-1)*Grid->nxC];

			}
			else {
				thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*BCThermal->value[(int) abs(NumThermal->map[ix-1+1+(iy-1+1)*(Grid->nxC+2)]) ]
								   + .25*(1.0-locX)*(1.0+locY)*BCThermal->value[(int) abs(NumThermal->map[ix-1+1+(iy+1 )*(Grid->nxC+2)]) ]
								   + .25*(1.0+locX)*(1.0+locY)*Physics->T[ix  +(iy)*Grid->nxC]
								   + .25*(1.0+locX)*(1.0-locY)*Physics->T[ix  +(iy-1)*Grid->nxC];

			}

			thisParticle = thisParticle->next;
		}
	}

	// Loop through inner nodes of the right
	ix = Grid->nxS-1;
	for (iy = 1; iy < Grid->nyS-1; iy++) {
		iNode = ix  + (iy  )*Grid->nxS;
		thisParticle = Particles->linkHead[iNode];
		// Loop through the particles in the shifted cell
		// ======================================
		while (thisParticle!=NULL) {
			locX = ((thisParticle->x-Grid->xmin)/dx - ix)*2.0;
			locY = ((thisParticle->y-Grid->ymin)/dy - iy)*2.0;

			if (BCStokes->SetupType==SimpleShearPeriodic) {
				ix = 0;
				thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*Physics->T[ix-1+(iy-1)*Grid->nxC]
								   + .25*(1.0-locX)*(1.0+locY)*Physics->T[ix-1+(iy)*Grid->nxC]
								   + .25*(1.0+locX)*(1.0+locY)*Physics->T[ix  +(iy)*Grid->nxC]
								   + .25*(1.0+locX)*(1.0-locY)*Physics->T[ix  +(iy-1)*Grid->nxC];
			}
			else {
				thisParticle->T +=   .25*(1.0-locX)*(1.0-locY)*Physics->T[ix-1+(iy-1)*Grid->nxC]
								   + .25*(1.0-locX)*(1.0+locY)*Physics->T[ix-1+(iy  )*Grid->nxC]
								   + .25*(1.0+locX)*(1.0+locY)*BCThermal->value[(int) abs(NumThermal->map[ix+1  +(iy +1 )*(Grid->nxC+2)]) ]
								   + .25*(1.0+locX)*(1.0-locY)*BCThermal->value[(int) abs(NumThermal->map[ix-1+1+(iy-1+1)*(Grid->nxC+2)]) ];
			}

			thisParticle = thisParticle->next;
		}
	}



	// Lower left corner
	// using linear shape functions for a triangular element
	ix = 0;
	iy = 0;
	iNode = ix  + (iy  )*Grid->nxS;// index of the node
	thisParticle = Particles->linkHead[iNode];
	while (thisParticle!=NULL) {
		locX = fabs((thisParticle->x-Grid->xmin)/dx - ix - 0.5);
		locY = fabs((thisParticle->y-Grid->ymin)/dy - iy - 0.5);
		thisParticle->T += (1.0-locX-locY)*Physics->T[ix  +(iy)*Grid->nxC]
						            + locX*BCThermal->value[(int) abs(NumThermal->map[0  +(1  )*(Grid->nxC+2)]) ]
									+ locY*BCThermal->value[(int) abs(NumThermal->map[1  +(0  )*(Grid->nxC+2)]) ];
		thisParticle = thisParticle->next;
	}

	// Lower right corner
	ix = Grid->nxS-1;
	iy = 0;
	iNode = ix  + (iy  )*Grid->nxS;// index of the node
	thisParticle = Particles->linkHead[iNode];
	while (thisParticle!=NULL) {
		locX = fabs((thisParticle->x-Grid->xmin)/dx - ix - 0.5);
		locY = fabs((thisParticle->y-Grid->ymin)/dy - iy - 0.5);
		thisParticle->T += (1.0-locX-locY)*Physics->T[ix  +(iy)*Grid->nxC]
									+ locX*BCThermal->value[(int) abs(NumThermal->map[Grid->nxC+2-1   +(1 )*(Grid->nxC+2)]) ]
									+ locY*BCThermal->value[(int) abs(NumThermal->map[Grid->nxC+2-1-1  +(0  )*(Grid->nxC+2)]) ];
		thisParticle = thisParticle->next;
	}

	// Upper left corner
	ix = 0;
	iy = Grid->nyS-1;
	iNode = ix  + (iy  )*Grid->nxS;// index of the node
	thisParticle = Particles->linkHead[iNode];
	while (thisParticle!=NULL) {
		locX = fabs((thisParticle->x-Grid->xmin)/dx - ix - 0.5);
		locY = fabs((thisParticle->y-Grid->ymin)/dy - iy - 0.5);
		thisParticle->T += (1.0-locX-locY)*Physics->T[ix  +(iy)*Grid->nxC]
								    + locX*BCThermal->value[(int) abs(NumThermal->map[0   +(Grid->nyC+2-1-1 )*(Grid->nxC+2)]) ]
									+ locY*BCThermal->value[(int) abs(NumThermal->map[1   +(Grid->nyC+2-1  )*(Grid->nxC+2)]) ];
		thisParticle = thisParticle->next;
	}

	// Upper right corner
	ix = 0;
	iy = Grid->nyS-1;
	iNode = ix  + (iy  )*Grid->nxS;// index of the node
	thisParticle = Particles->linkHead[iNode];
	while (thisParticle!=NULL) {
		locX = fabs((thisParticle->x-Grid->xmin)/dx - ix - 0.5);
		locY = fabs((thisParticle->y-Grid->ymin)/dy - iy - 0.5);
		thisParticle->T += (1.0-locX-locY)*Physics->T[ix  +(iy)*Grid->nxC]
									+ locX*BCThermal->value[(int) abs(NumThermal->map[Grid->nxC+2-1   +(Grid->nyC+2-1-1 )*(Grid->nxC+2)]) ]
									+ locY*BCThermal->value[(int) abs(NumThermal->map[Grid->nxC+2-1-1   +(Grid->nyC+2-1  )*(Grid->nxC+2)]) ];
		thisParticle = thisParticle->next;
	}


*/

























}

























void Physics_interpFromCellToNode(Grid* Grid, compute* CellValue, compute* NodeValue, int BCType)
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


/*

	// CellValue extrapolated on the lower boundary
	// ======================================
	// o: centered CellValue
	// x: CellValue extrapolated (in 1D) fron the o nodes
	// X: value interpolated between the two x
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1b   and 2b
	// |       |       |
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1a   and 2a
	// |       |       |
	// | - x - X - x - |       nodes tempa and tempb
	//
	iy = 0;
	compute temp1, temp2;
	int i1a, i1b, i2a, i2b;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+ iy   *Grid->nxC;
		i2b =  ix   +(iy+1)*Grid->nxC;
		i2a =  ix   + iy   *Grid->nxC;


		//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		//NodeValue[I] = (temp1+temp2)/2;

		// No extrapolation just interpolation: average of the two closest cell centers
		NodeValue[I] = (CellValue[i1a] + CellValue[i2a])/2;
	}
	// CellValue extrapolated on the upper boundary
	// ======================================
	//   x  X  x
	//  1a    2a
	//  1b    2b
	iy = Grid->nyS-1;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		i2a =  ix   +(iy-1) *Grid->nxC;
		//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		//NodeValue[I] = (temp1+temp2)/2;
		NodeValue[I] = (CellValue[i1a] + CellValue[i2a])/2;
	}


	if (BCType!=SimpleShearPeriodic) { // not periodic
		// CellValue extrapolated on the left boundary
		// ======================================
		//  x 1a   1b
		//  X
		//  x 2a   2b
		ix = 0;
		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			I = ix + iy*Grid->nxS;
			i1a =  ix   +(iy  )*Grid->nxC;
			i2a =  ix   +(iy-1)*Grid->nxC;
			//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			//NodeValue[I] = (temp1+temp2)/2;
			NodeValue[I] = (CellValue[i1a] + CellValue[i2a])/2;
		}

		// CellValue extrapolated on the right boundary
		// ======================================
		//  1b   1a x
		//          X
		//  2b   2a x
		ix = Grid->nxS-1;
		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			I = ix + iy*Grid->nxS;
			i1a = (ix-1)+(iy  )*Grid->nxC;
			i2a = (ix-1)+(iy-1)*Grid->nxC;
			//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			//NodeValue[I] = (temp1+temp2)/2;
			NodeValue[I] = (CellValue[i1a] + CellValue[i2a])/2;
		}

		// Lower left corner
		//          1b
		//      1a
		//   X
		ix = 0; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy+1)*Grid->nxC;
		i1a =  ix   +(iy  )*Grid->nxC;
		//NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		// Value of the closest
		NodeValue[I] = CellValue[i1a];

		// Lower right corner
		//  1b
		//      1a
		//          X
		ix = Grid->nxS-1; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+(iy  )*Grid->nxC;
		//NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		NodeValue[I] = CellValue[i1a];

		// Upper left corner
		//  X
		//      1a
		//          1b
		ix = 0; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy-2)*Grid->nxC;
		i1a =  ix   +(iy-1)*Grid->nxC;
		//NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		NodeValue[I] = CellValue[i1a];

		// Upper right corner
		//          X
		//      1a
		//  1b
		ix = Grid->nxS-1; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		//NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		NodeValue[I] = CellValue[i1a];


	}
	else { // if periodic boundaries

		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			// Left and right boundary
			ix = 0;
			I = ix + iy*Grid->nxS;
			iNW = (ix+Grid->nxC-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix+Grid->nxC-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			NodeValue[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
			NodeValue[I+Grid->nxS-1] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}

		// Upper left and right corners
		// ======================================
		//   x  X  x
		//  1a    2a
		//  1b    2b
		iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (Grid->nxC-1)+(iy-2)*Grid->nxC;
		i1a = (Grid->nxC-1)+(iy-1)*Grid->nxC;
		i2b =  0   +(iy-2)*Grid->nxC;
		i2a =  0   +(iy-1) *Grid->nxC;
		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		NodeValue[I] = (temp1+temp2)/2;
		NodeValue[I+Grid->nxS-1] = (temp1+temp2)/2;


		// Lower left and right corners
		// ======================================
		//  1b    2b
		//  1a    2a
		//   x  X  x
		iy = 0;
		compute temp1, temp2;
		int i1a, i1b, i2a, i2b;
		I = ix + iy*Grid->nxS;
		i1b = (Grid->nxC-1)+(iy+1)*Grid->nxC;
		i1a = (Grid->nxC-1)+ iy   *Grid->nxC;
		i2b =  0   +(iy+1)*Grid->nxC;
		i2a =  0   + iy   *Grid->nxC;
		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		NodeValue[I] = (temp1+temp2)/2;
		NodeValue[I+Grid->nxS-1] = (temp1+temp2)/2;


	}

*/
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


/*
	// Set P
	// =========================
	C = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			I = ix + iy*Grid->nxC + Grid->nVxTot + Grid->nVyTot;
			InoDir = Numbering->map[I];
			if (InoDir>=0) { // Not a Dirichlet node
				Physics->P[C] = EqSystem->x[InoDir];
			} else {
				IBC = abs(InoDir)-1;
				Physics->P[C] = BC->value[ IBC ];
			}
			C++;
		}
	}


*/
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
				printf("%.2f  ", Physics->Vx[C]);
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
				printf("%.2f  ", Physics->Vy[C]);
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
				printf("%.2f  ", Physics->P[C]);
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

	/*
	// Corner values
	compute value;
	// Lower left
	ix = 0; iy = 0;
	value = (Physics->T[1] + Physics->T[Grid->nxEC])/2;
	INeigh = NumThermal->map[  ix+1 + (iy+1)*Grid->nxEC  ];
	Physics->T[ ix + iy*Grid->nxEC ] = 2.0*value - EqThermal->x[INeigh];
	printf("LL Corner: T = %.3f\n", Physics->T[ ix + iy*Grid->nxEC ]);

	// Lower right
	ix = Grid->nxEC-1; iy = 0;
	value = (Physics->T[Grid->nxEC-2] + Physics->T[2*Grid->nxEC-1])/2;
	INeigh = NumThermal->map[  ix-1 + (iy+1)*Grid->nxEC  ];
	Physics->T[ ix + iy*Grid->nxEC ] = 2.0*value- EqThermal->x[INeigh];
	printf("LR Corner: T = %.3f\n", Physics->T[ ix + iy*Grid->nxEC ]);


	// Upper left
	ix = 0;
	iy = Grid->nyEC-1;
	value = (Physics->T[(Grid->nyEC-2)*Grid->nxEC] + Physics->T[(Grid->nyEC-1)*Grid->nxEC+1])/2;
	INeigh = NumThermal->map[  ix+1 + (iy-1)*Grid->nxEC  ];
	Physics->T[ ix + iy*Grid->nxEC ] = 2.0*value- EqThermal->x[INeigh];
	printf("UL Corner: T = %.3f, value = %.3f\n", Physics->T[ ix + iy*Grid->nxEC ], value);

	// Upper right
	ix = Grid->nxEC-1;
	iy = Grid->nyEC-1;
	value = (Physics->T[(Grid->nyEC)*Grid->nxEC-2] + Physics->T[(Grid->nyEC-1)*Grid->nxEC-1])/2;
	INeigh = NumThermal->map[  ix-1 + (iy-1)*Grid->nxEC  ];
	Physics->T[ ix + iy*Grid->nxEC ] = 2.0*value- EqThermal->x[INeigh];
	printf("UR Corner: T = %.3f\n", Physics->T[ ix + iy*Grid->nxEC ]);
	*/



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


void Physics_computeStrainRateInvariant(Physics* Physics, Grid* Grid, compute* StrainRateInvariant)
{
// Definition of second invariant: // E_II = sqrt( Eps_xx^2 + Eps_xy^2  );
// Declarations
// =========================
int ix, iy, I, iNode, Ix, Iy, IE;
compute dVxdy, dVydx, dVxdx;
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


		StrainRateInvariant[IE] = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +    dVxdx*dVxdx  );
		//StrainRateInvariant[IE] = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))  );
		if (StrainRateInvariant[IE]==0) {
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








// Non linear rheology code, to put somewhere else

/*


switch (MatProps->flowLaw[phase]) {
				case LinearViscous:
					locEta = MatProps->eta0[phase];
					break;


				case PowerLawViscous:
					// Compute Eps_xy at the four nodes of the cell
					// 1. Sum contributions
					dVxdy = 0;
					dVydx = 0;
					for (iNode = 0; iNode < 4; ++iNode) {
						IX = ix+IxMod[iNode];
						IY = iy+IyMod[iNode];
						dVxdy += ( Physics->Vx[(IX  )+(IY+1)*Grid->nxVx]
											   - Physics->Vx[(IX  )+(IY  )*Grid->nxVx] )/Grid->dy;

						dVydx += ( Physics->Vy[(IX+1)+(IY  )*Grid->nxVy]
											   - Physics->Vy[(IX  )+(IY  )*Grid->nxVy] )/Grid->dx;
					}
					// 2. Average
					dVxdy /= 4;
					dVydx /= 4;

					dVxdx = (Physics->Vx[(ix+1) + (iy+1)*Grid->nxVx]
										 - Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx])/Grid->dx;
					dVydy = (Physics->Vy[(ix+1) + (iy+1)*Grid->nxVy]
										 - Physics->Vy[(ix+1) + (iy  )*Grid->nxVy])/Grid->dy;


					// Local Strain rate invariant
					locEps_II = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +    0.5*dVxdx*dVxdx + 0.5*dVydy*dVydy );

					//locEps_II = sqrt(   (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))  );
					//printf("ix=%i, iy=%i, locEps_II = %.3e, EpsRef = %.3e, dVxdx2=%.4e, dVydy2=%.4e\n",ix, iy, locEps_II,Physics->epsRef, dVxdx*dVxdx, dVydy*dVydy );
					//locEps_II = dVxdx;
					if (locEps_II==0)
						locEps_II = Physics->epsRef;
					locEta = MatProps->eta0[phase] * pow(  (locEps_II/Physics->epsRef)    ,   1.0/MatProps->n[phase]-1.0 );
					break;


				default:
					printf("Unknown flow law %i, for phase %i", MatProps->flowLaw[phase], phase);
					exit(0);
				}

*/

