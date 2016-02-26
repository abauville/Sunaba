/*
 * BC.c
 *
 *  Created on: Feb 23, 2016
 *      Author: abauville
 */

#include "stokes.h"

//==========================================================================
//
//                            BOUNDARY INDEXING
//
//==========================================================================
void BC_set(BC* BC, Grid* Grid, EqSystem* EqSystem, Physics* Physics)
{
	// Declarations
	// =========================
	int I, C, C2;
	int i;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	// Set and fill Dirichlet boundary conditions
	// =======================================
	BC->nPDir = 1;
	BC->nDir = 2*Grid->nxVy + 2*Grid->nyVx + BC->nPDir; // Vx eq + Vy Eq + P eq
	EqSystem->nEq = EqSystem->nEqIni - BC->nDir;


	BC->listDir    = (int*)     malloc( BC->nDir * sizeof(  int  ));
	BC->valueDir   = (compute*) malloc( BC->nDir * sizeof(compute));

	// Assign VxL, VyB etc... to Dirichlet values
	BC_updateDir(BC, Grid);

	// Set and fill Dirichlet Pressure boundary conditions
	// =======================================

	BC->listDir[BC->nDir-1]  = EqSystem->nEqIni-1;
	BC->valueDir[BC->nDir-1] = 0.0;


	// Set and fill Neumann boundary conditions (In the numbering with Dirichlet)
	// =======================================
	BC->nNeu = ((Grid->nxVx-2)*2 + (Grid->nyVy-2)*2);
	BC->listNeu       = (int*)     malloc( BC->nNeu * sizeof(  int  ));
	BC->listNeuNeigh  = (int*)     malloc( BC->nNeu * sizeof(  int  ));
	BC->coeffNeu      = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->coeffNeuNeigh = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->valueNeu      = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->isNeu         = (bool*)    malloc( EqSystem->nEqIni   * sizeof(bool)   );

	I = 0;
	C = Grid->nVxTot + Grid->nxVy;
	C2 = Grid->nxS;
	for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
		BC->listNeu[I]          = C;
		BC->listNeuNeigh[I]    = C+1;

		BC->coeffNeu[I]         = -Physics->etaShear[C2]/dx/dx;
		BC->coeffNeuNeigh[I]   =   Physics->etaShear[C2]/dx/dx;

		BC->valueNeu[I]         =  0.0;

		I++;
		C = C+Grid->nxVy;
		C2 += Grid->nxS;
	}

	C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
	C2 = 2*Grid->nxS-1;
	for (i=0;i<Grid->nyVy-2;i++){ // Vy Right

		BC->listNeu[I]          = C;
		BC->listNeuNeigh[I]    = C-1;

		BC->coeffNeu[I]         = -Physics->etaShear[C2]/dx/dx;
		BC->coeffNeuNeigh[I]   =   Physics->etaShear[C2]/dx/dx;

		BC->valueNeu[I]         =  -0.0;

		I++;
		C = C+Grid->nxVy;
		C2+=Grid->nxS;
	}

	C = 1;
	C2 = 1;
	for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
		BC->listNeu[I]          = C;
		BC->listNeuNeigh[I]    = C+Grid->nxVx;

		BC->coeffNeu[I]         = -Physics->eta[C2]/dy/dy;
		BC->coeffNeuNeigh[I]   =  Physics->eta[C2]/dy/dy;

		BC->valueNeu[I]         =  0.0;

		I++;
		C = C+1;
		C2 += 1;
	}

	C = Grid->nxVx*(Grid->nyVx-1)+1;
	C2 = (Grid->nyC-1)*Grid->nxC+1;
	for (i=0;i<Grid->nxVx-2;i++){ // Vx Top

		BC->listNeu[I]          = C;
		BC->listNeuNeigh[I]    = C-Grid->nxVx;

		BC->coeffNeu[I]         = -Physics->eta[C2]/dy/dy;
		BC->coeffNeuNeigh[I]   =   Physics->eta[C2]/dy/dy;

		BC->valueNeu[I]         =  0.0;

		I++;
		C = C+1;
		C2+=1;
	}


	// Fill BC->isNeu
	for (i=0; i<EqSystem->nEqIni; i++) {
		BC->isNeu[i] = false;
	}
	for (i=0; i<BC->nNeu; i++) {
		BC->isNeu[ BC->listNeu[i] ] = true;
	}

	if (DEBUG) {
		printf("======= BC->isNeu =======\n");
		for (i=0; i<EqSystem->nEqIni; i++) {
			printf("%i  ", BC->isNeu[i]);
		}
		printf("\n");
	}

	if (DEBUG) {
		printf("=== List Dir ===\n");
		printListi(BC->listDir,BC->nDir);
		printf("=== List Neu ===\n");
		printListi(BC->listNeu,BC->nNeu);
		printf("=== List Neu Neigh ===\n");
		printListi(BC->listNeuNeigh,BC->nNeu);
	}






}

void BC_updateDir(BC* BC, Grid* Grid)
{
	int I, C, i;
	I = 0;

		C = 0;
		for (i=0; i<Grid->nyVx; i++) { // Vx Left
			BC->listDir[I] = C;

			BC->valueDir[I] = BC->VxL;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			BC->listDir[I] = C;

			BC->valueDir[I] = BC->VxR;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			BC->listDir[I] = C;

			BC->valueDir[I] = BC->VyB;

			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			BC->listDir[I] = C;

			BC->valueDir[I] = BC->VyT;

			I++;
			C += 1;
		}
}



