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
void BC_init(BC* BC, Grid* Grid, EqSystem* EqSystem, Physics* Physics)
{


	// Set and fill Dirichlet boundary conditions
	// =======================================
	switch (BC->SetupType) {
	case PureShear:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		BC->nPDir 	= 1;
		BC->nDir 	= 2*Grid->nxVy + 2*Grid->nyVx + BC->nPDir; // Vx eq + Vy Eq + P eq
		BC->nNeu 	= ((Grid->nxVx-2)*2 + (Grid->nyVy-2)*2);
		BC->n 		= BC->nDir + BC->nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;
	case SimpleShearPeriodic:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		BC->nPDir = 1;
		BC->nDir = 2*Grid->nxVy + 2*Grid->nxVx + BC->nPDir; // Vx eq + Vy Eq + P eq
		// no Neumann nodes for this setup
		BC->nNeu = 0;
		BC->n 		= BC->nDir + BC->nNeu;

		int nPeriod = Grid->nyVx-2 + 2*(Grid->nyVy-2);

		EqSystem->nEq = EqSystem->nEqIni - BC->nDir - nPeriod;
		//EqSystem->nEq = EqSystem->nEqIni - BC->nDir;
		break;
	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;

	}









	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*) malloc ( BC->n * sizeof(BCType));
	BC_update(BC, Grid);





	// Set and fill Dirichlet Pressure boundary conditions
	// =======================================
	BC->list[BC->n-1]  = EqSystem->nEqIni-1;
	//BC->listDir[BC->nDir-1]  = Grid->nVxTot+Grid->nVyTot;
	BC->value[BC->n-1] = 0.0;
	BC->type[BC->n-1] = 0;





	// Check
	// =======================================

	if (DEBUG) {
		printf("=== BC list ===\n");
		printListi(BC->list,BC->n);
		printf("=== BC value ===\n");
		printListd(BC->value,BC->n);
		printf("=== BC type ===\n");
		int i;
		for (i=0;i<BC->n;i++) {
			printf("%d  ", BC->type[i]);
		}
		printf("\n");
	}

}






void BC_update(BC* BC, Grid* Grid)
{
	int C, i;
	int I = 0;

	switch (BC->SetupType) {
	case PureShear:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		BC->VxL =  BC->backStrainRate*Grid->xmin;
		BC->VxR =  BC->backStrainRate*Grid->xmax;
		BC->VyB = -BC->backStrainRate*Grid->ymin;
		BC->VyT = -BC->backStrainRate*Grid->ymax;

		C = 0;
		for (i=0; i<Grid->nyVx; i++) { // Vx Left
			BC->list[I] = C;

			BC->value[I] = BC->VxL;
			BC->type[I] = Dirichlet;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			BC->list[I] = C;
			BC->value[I] = BC->VxR;
			BC->type[I] = Dirichlet;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			BC->list[I] = C;

			BC->value[I] = BC->VyB;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			BC->list[I] = C;
			BC->value[I] = BC->VyT;
			BC->type[I] = Dirichlet;

			I++;
			C += 1;
		}




		// Neumann
		// =======================================


		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			BC->list[I]          = C;
			BC->value[I]         =  0.0;
			BC->type[I] 		 = NeumannGhost;
			I++;
			C = C+Grid->nxVy;
		}

		C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right

			BC->list[I]          = C;
			BC->value[I]         = 0.0;
			BC->type[I] 		 = NeumannGhost;

			I++;
			C = C+Grid->nxVy;
		}

		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
			BC->list[I]          = C;
			BC->value[I]         = 0.0;
			BC->type[I] = NeumannGhost;

			I++;
			C = C+1;
		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top

			BC->list[I]         = C;
			BC->value[I]         = 0.0;
			BC->type[I] = NeumannGhost;

			I++;
			C = C+1;
		}















		break;

	case SimpleShearPeriodic:
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
			BC->list[I] = C;
			BC->value[I] = BC->VyB;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			BC->list[I] = C;
			BC->value[I] = BC->VyT;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}

		// Top and bottom Vx
		// =======================================
		C = 0;
		for (i=0; i<Grid->nxVx; i++) { // Vx Bottom
			BC->list[I] = C;
			BC->value[I] = 2*BC->VxB; // factor 2 because it's a dirichlet condition on a ghost node
			BC->type[I] = DirichletGhost;
			I++;
			C += 1;
		}


		C = Grid->nxVx*(Grid->nyVx-1);
		for (i=0; i<Grid->nxVx; i++) { // Vx Top
			BC->type[I] = DirichletGhost;
			BC->list[I] = C;
			BC->value[I] = 2*BC->VxT;
			I++;
			C += 1;
		}


		break;

		// no Neumann nodes for this setup


	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;
	}
}




