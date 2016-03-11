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
	int i;

	// Set and fill Dirichlet boundary conditions
	// =======================================
	switch (BC->SetupType) {
	case 0:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		BC->nPDir = 1;
		BC->nDir = 2*Grid->nxVy + 2*Grid->nyVx + BC->nPDir; // Vx eq + Vy Eq + P eq
		BC->nNeu = ((Grid->nxVx-2)*2 + (Grid->nyVy-2)*2);

		EqSystem->nEq = EqSystem->nEqIni - BC->nDir;

		break;
	case 1:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		BC->nPDir = 1;
		BC->nDir = 2*Grid->nxVy + 2*Grid->nxVx + BC->nPDir; // Vx eq + Vy Eq + P eq
		// no Neumann nodes for this setup
		BC->nNeu = 0;

		int nPeriod = Grid->nyVx-2 + 2*(Grid->nyVy-2);

		EqSystem->nEq = EqSystem->nEqIni - BC->nDir - nPeriod;
		//EqSystem->nEq = EqSystem->nEqIni - BC->nDir;
		break;
	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;

	}









	BC->listDir    = (int*)     malloc( BC->nDir * sizeof(  int  ));
	BC->valueDir   = (compute*) malloc( BC->nDir * sizeof(compute));

	BC_updateDir(BC, Grid);





	// Set and fill Dirichlet Pressure boundary conditions
	// =======================================
	BC->listDir[BC->nDir-1]  = EqSystem->nEqIni-1;
	//BC->listDir[BC->nDir-1]  = Grid->nVxTot+Grid->nVyTot;
	BC->valueDir[BC->nDir-1] = 0.0;






	// Set and fill Neumann boundary conditions (In the numbering with Dirichlet)
	// =======================================

	BC->listNeu       = (int*)     malloc( BC->nNeu * sizeof(  int  ));
	BC->listNeuNeigh  = (int*)     malloc( BC->nNeu * sizeof(  int  ));
	BC->coeffNeu      = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->coeffNeuNeigh = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->valueNeu      = (compute*) malloc( BC->nNeu * sizeof(compute));
	BC->isNeu         = (bool*)    malloc( EqSystem->nEqIni   * sizeof(bool)   );

	BC_numberNeu(BC, Grid, EqSystem);
	BC_updateNeuCoeff(BC, Grid, Physics);




	// Check
	// =======================================
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
	int C, i;
	int I = 0;

	switch (BC->SetupType) {
	case 0:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		BC->VxL = -BC->backStrainRate*Grid->xmin;
		BC->VxR = -BC->backStrainRate*Grid->xmax;
		BC->VyB =  BC->backStrainRate*Grid->ymin;
		BC->VyT =  BC->backStrainRate*Grid->ymax;

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
		break;

	case 1:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		BC->VxB =  0*1.0*BC->backStrainRate*Grid->ymin;
		BC->VxT =  2*1.0*BC->backStrainRate*Grid->ymax;
		BC->VyB =  0;
		BC->VyT =  0;

		// Top and bottom Vy
		// =======================================
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

		// Top and bottom Vx
		// =======================================
		C = 0;
		for (i=0; i<Grid->nxVx; i++) { // Vy Bottom
			BC->listDir[I] = C;
			BC->valueDir[I] = 2*BC->VxB; // factor 2 because it's a dirichlet condition on a ghost node
			I++;
			C += 1;
		}


		C = Grid->nxVx*(Grid->nyVx-1);
		for (i=0; i<Grid->nxVx; i++) { // Vy Top
			BC->listDir[I] = C;
			BC->valueDir[I] = 2*BC->VxT;
			I++;
			C += 1;
		}






		break;

	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;
	}
}



void BC_numberNeu(BC* BC, Grid* Grid, EqSystem* EqSystem)
{
	int i, C;
	int I = 0;
	switch (BC->SetupType) {
	case 0:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			BC->listNeu[I]          = C;
			BC->listNeuNeigh[I]    	= C+1;
			I++;
			C = C+Grid->nxVy;
		}

		C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right

			BC->listNeu[I]          = C;
			BC->listNeuNeigh[I]    = C-1;

			I++;
			C = C+Grid->nxVy;
		}

		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
			BC->listNeu[I]          = C;
			BC->listNeuNeigh[I]    = C+Grid->nxVx;

			I++;
			C = C+1;
		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top

			BC->listNeu[I]         = C;
			BC->listNeuNeigh[I]    = C-Grid->nxVx;

			I++;
			C = C+1;
		}

		// Fill BC->isNeu
		for (i=0; i<EqSystem->nEqIni; i++) {
			BC->isNeu[i] = false;
		}
		for (i=0; i<BC->nNeu; i++) {
			BC->isNeu[ BC->listNeu[i] ] = true;
		}
		break;

	case 1:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================

		// no Neumann nodes for this setup

		for (i=0; i<EqSystem->nEqIni; i++) {
			BC->isNeu[i] = false;
		}

		break;

	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;
		break;
	}
}



void BC_updateNeuCoeff(BC* BC, Grid* Grid, Physics* Physics)
{
	compute dx = Grid->dx;
	compute dy = Grid->dy;
	int i, I, C2;

	switch (BC->SetupType) {
	case 0:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		I = 0;
		C2 = Grid->nxS;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			BC->coeffNeu[I]         = -Physics->etaShear[C2]/dx/dx;
			BC->coeffNeuNeigh[I]   =   Physics->etaShear[C2]/dx/dx;

			BC->valueNeu[I]         =  0.0;

			I++;
			C2 += Grid->nxS;
		}

		C2 = 2*Grid->nxS-1;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right

			BC->coeffNeu[I]         = -Physics->etaShear[C2]/dx/dx;
			BC->coeffNeuNeigh[I]   =   Physics->etaShear[C2]/dx/dx;

			BC->valueNeu[I]         =  -0.0;

			I++;
			C2+=Grid->nxS;
		}

		C2 = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom

			BC->coeffNeu[I]         = -Physics->eta[C2]/dy/dy;
			BC->coeffNeuNeigh[I]   =  Physics->eta[C2]/dy/dy;

			BC->valueNeu[I]         =  0.0;

			I++;
			C2 += 1;
		}

		C2 = (Grid->nyC-1)*Grid->nxC+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top

			BC->coeffNeu[I]         = -Physics->eta[C2]/dy/dy;
			BC->coeffNeuNeigh[I]   =   Physics->eta[C2]/dy/dy;

			BC->valueNeu[I]         =  0.0;

			I++;
			C2+=1;
		}

		break;

	case 1:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================

		// no Neumann nodes for this setup

		break;

	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;
	}


}

