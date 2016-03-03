/*
 * Physics.c
 *
 *  Created on: Feb 24, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps)
{
	// Declarations
	// =========================
	int iCell, iP, ix, iy, i;
	coord locX, locY;

	coord dx = Grid->dx;
	coord dy = Grid->dy;
	compute* sumOfWeights = (compute*) malloc(Grid->nCTot * sizeof(compute));

	compute weight;
	int phase;

	int nxC = Grid->nxC;
	int xMod[4], yMod[4], Ix[4], Iy[4];
	int I;

	xMod[0] =  1; yMod[0] =  1;
	xMod[1] = -1; yMod[1] =  1;
	xMod[2] = -1; yMod[2] = -1;
	xMod[3] =  1; yMod[3] = -1;


	// Reinitialize Physics array
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		Physics->eta[iCell] = 0;
		Physics->rho[iCell] = 0;
		sumOfWeights[iCell] = 0;
	}

	printf("Initialized\n");

	//int quadrant = 0;

	// Loop through inner cells
	// ========================
	iCell = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			iCell = ix  + (iy  )*nxC;
			iP = Particles->linkHead[iCell];

			// Loop through the particles in the cell
			// ======================================
			while (iP!=-1) {
				locX = (Particles->xy[2*iP  ]-Grid->xmin)/dx - ix;
				locY = (Particles->xy[2*iP+1]-Grid->ymin)/dy - iy;

				phase = Particles->phase[iP];


				// Get the index of the neighbours
				if (locX<0.5) {
					if (locY<0.5) { // Lower left quadrant
						//quadrant = 0;
						Ix[0] = ix  ; Iy[0] = iy  ;
						Ix[1] = ix-1; Iy[1] = iy  ;
						Ix[2] = ix-1; Iy[2] = iy-1;
						Ix[3] = ix  ; Iy[3] = iy-1;
					} else { 		// Upper left quadrant
						//quadrant = 1;
						locY =  (1-locY);
						Ix[0] = ix  ; Iy[0] = iy  ;
						Ix[1] = ix-1; Iy[1] = iy  ;
						Ix[2] = ix-1; Iy[2] = iy+1;
						Ix[3] = ix  ; Iy[3] = iy+1;
					}
				}
				else {
					if (locY<0.5) { // Lower right quadrant
						//quadrant = 2;
						locX = (1-locX);
						Ix[0] = ix  ; Iy[0] = iy  ;
						Ix[1] = ix+1; Iy[1] = iy  ;
						Ix[2] = ix+1; Iy[2] = iy-1;
						Ix[3] = ix  ; Iy[3] = iy-1;
					} else { 		// Upper right quadrant
						//quadrant = 3;
						locX = (1-locX);
						locY = (1-locY);
						Ix[0] = ix  ; Iy[0] = iy  ;
						Ix[1] = ix+1; Iy[1] = iy  ;
						Ix[2] = ix+1; Iy[2] = iy+1;
						Ix[3] = ix  ; Iy[3] = iy+1;
					}
				}

				// Add contribution of the particle to each of the four cells that its area overlaps
				// the contribution is the non-dimensional area (Total Area of the particle: dx*dy/(dx*dy))
				//if (DEBUG)
				//printf("ix=%i, iy=%i, iP=%i, quadrant=%i,  phase=%i, eta0=%2f, rho0=%.2f =====\n",ix,iy, iP, quadrant, phase, MatProps->eta0[phase], MatProps->rho0[phase]);
				for (i = 0; i < 4; ++i) {
					if (Ix[i]>=0 && Ix[i]<Grid->nxC) { // Check for boundaries
						if (Iy[i]>=0 && Iy[i]<Grid->nyC) {
							I = Ix[i] + Iy[i] * nxC;
							weight = fabs((locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5));
							Physics->eta[I] += MatProps->eta0[phase] * weight;
							Physics->rho[I] += MatProps->rho0[phase] * weight;
							sumOfWeights[I] += weight;
							//if (DEBUG)
							//printf("i=%i, Ix[i]=%i, Iy[i]=%i, weight=%.2f, locX=%.2f, locY=%.2f, A=%.2f, B=%.2f\n",i,Ix[i], Iy[i], weight, locX, locY, (locX + xMod[i]*0.5), (locY + yMod[i]*0.5) );
						}
					}

				}
				//if (DEBUG)
				//printf("\n");

				iP = Particles->linkNext[iP];
			}
		}
	}











	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		Physics->eta[iCell] /= sumOfWeights[iCell];
		Physics->rho[iCell] /= sumOfWeights[iCell];
	}


	if (DEBUG) {
		printf("=== Check eta 1 ===\n");
		int C = 0;
		//int ix, iy;
		for (iy = 0; iy < Grid->nyC; ++iy) {
			for (ix = 0; ix < Grid->nxC; ++ix) {
				printf("%.3f  ", Physics->eta[C]);
				C++;
			}
			printf("\n");
		}
	}

	free(sumOfWeights);



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
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			I = ix + iy*Grid->nxS;
			iNW = (ix-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			NodeValue[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}




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
	// |   o   |   o   |       nodes 1a   and 1b
	// |       |       |
	// | - - - | - - - | -
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


		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		NodeValue[I] = (temp1+temp2)/2;
	}
	// CellValue extrapolated on the upper boundary
	// ======================================
	//   x  X  x
	//  1a    2a
	//  1b    2b
	iy = Grid->nyS-1;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		i2b =  ix   +(iy-2)*Grid->nxC;
		i2a =  ix   +(iy-1) *Grid->nxC;
		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		NodeValue[I] = (temp1+temp2)/2;
	}


	if (BCType!=1) { // not periodic
		// CellValue extrapolated on the left boundary
		// ======================================
		//  x 1a   1b
		//  X
		//  x 2a   2b
		ix = 0;
		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			I = ix + iy*Grid->nxS;
			i1b = (ix+1)+(iy  )*Grid->nxC;
			i1a =  ix   +(iy  )*Grid->nxC;
			i2b = (ix+1)+(iy-1)*Grid->nxC;
			i2a =  ix   +(iy-1)*Grid->nxC;
			temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			NodeValue[I] = (temp1+temp2)/2;
		}

		// CellValue extrapolated on the right boundary
		// ======================================
		//  1b   1a x
		//          X
		//  2b   2a x
		ix = Grid->nxS-1;
		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			I = ix + iy*Grid->nxS;
			i1b = (ix-2)+(iy  )*Grid->nxC;
			i1a = (ix-1)+(iy  )*Grid->nxC;
			i2b = (ix-2)+(iy-1)*Grid->nxC;
			i2a = (ix-1)+(iy-1)*Grid->nxC;
			temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			NodeValue[I] = (temp1+temp2)/2;
		}

		// Lower left corner
		//          1b
		//      1a
		//   X
		ix = 0; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy+1)*Grid->nxC;
		i1a =  ix   +(iy  )*Grid->nxC;
		NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;

		// Lower right corner
		//  1b
		//      1a
		//          X
		ix = Grid->nxS-1; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+(iy  )*Grid->nxC;
		NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;

		// Upper left corner
		//  X
		//      1a
		//          1b
		ix = 0; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy-2)*Grid->nxC;
		i1a =  ix   +(iy-1)*Grid->nxC;
		NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;

		// Upper right corner
		//          X
		//      1a
		//  1b
		ix = Grid->nxS-1; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		NodeValue[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;



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

	if (DEBUG) {
		int C = 0;
		printf("=== Check CellValue ===\n");
		for (iy = 0; iy < Grid->nyC; ++iy) {
			for (ix = 0; ix < Grid->nxC; ++ix) {
				printf("%.2f  ", CellValue[C]);
				C++;
			}
			printf("\n");
		}


		C = 0;
		printf("=== Check U ===\n");
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.2f  ", NodeValue[C]);
				C++;
			}
			printf("\n");
		}

	}

}


void Physics_set_VxVyP_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, compute* sol)
{
	// Declarations
	// =========================
	int ix, iy, i;
	int I, C;
	int InoDir;
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
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			I = ix + iy*Grid->nxVx;
			InoDir = Numbering->map[I];
			if (InoDir!=-1) { // Not a Dirichlet node
				Physics->Vx[C] = sol[InoDir];
			} else {
				Physics->Vx[C] = BC->valueDir[ findi(BC->listDir,BC->nDir,I) ];
			}
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

			if (InoDir!=-1) { // Not a Dirichlet node
				Physics->Vy[C] = sol[InoDir];
			} else {
				Physics->Vy[C] = BC->valueDir[ findi(BC->listDir,BC->nDir,I) ];
			}
			if (Physics->Vy[C]*Physics->Vy[C] > maxVy)
				maxVy = Physics->Vy[C]*Physics->Vy[C];
			C++;
		}
	}


	// Set P
	// =========================
	C = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			I = ix + iy*Grid->nxC + Grid->nVxTot + Grid->nVyTot;
			InoDir = Numbering->map[I];
			if (InoDir!=-1) { // Not a Dirichlet node
				Physics->P[C] = sol[InoDir];
			} else {
				Physics->P[C] = BC->valueDir[ findi(BC->listDir,BC->nDir,I) ];
			}
			C++;
		}
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
		for (iy = 0; iy < Grid->nyC; ++iy) {
			for (ix = 0; ix < Grid->nxC; ++ix) {
				printf("%.2f  ", Physics->P[C]);
				C++;
			}
			printf("\n");
		}
	}












}



