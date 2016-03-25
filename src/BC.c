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
void BC_initStokes(BC* BC, Grid* Grid, EqSystem* EqSystem)
{

	int nPDir, nDir, nNeu;

	// Set and fill Dirichlet boundary conditions
	// =======================================
	switch (BC->SetupType) {
	case PureShear:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		nPDir 	= 1;
		nDir 	= 2*Grid->nxVy + 2*Grid->nyVx + nPDir; // Vx eq + Vy Eq + P eq
		nNeu 	= ((Grid->nxVx-2)*2 + (Grid->nyVy-2)*2);
		BC->n 		= nDir + nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;
	case SimpleShearPeriodic:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		nPDir = 1;
		nDir = 2*Grid->nxVy + 2*Grid->nxVx + nPDir; // Vx eq + Vy Eq + P eq
		// no Neumann nodes for this setup
		nNeu = 0;
		BC->n 		= nDir + nNeu;

		int nPeriod = Grid->nyVx-2 + 2*(Grid->nyVy-2);

		EqSystem->nEq = EqSystem->nEqIni - nDir - nPeriod;
		//EqSystem->nEq = EqSystem->nEqIni - BC->nDir;
		break;
	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;

	}


	if (UPPER_TRI) {
		EqSystem->nRow = EqSystem->nEq + nPDir - Grid->nCTot;
	}
	else {
		EqSystem->nRow = EqSystem->nEq;
	}




	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*) malloc ( BC->n * sizeof(BCType));
	BC_updateStokes(BC, Grid);


	// Set and fill Dirichlet Pressure boundary conditions
	// =======================================
	BC->list[BC->n-1]  = EqSystem->nEqIni-1;
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


void BC_initThermal(BC* BC, Grid* Grid, EqSystem* EqSystem)
{
	int nDir, nNeu;

	// Set and fill Dirichlet boundary conditions
	// =======================================
	switch (BC->SetupType) {
	case PureShear:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		// Dirichlet on upper and lower
		// Neumann on the sides
		nDir 	= 2*Grid->nxC; // Vx eq + Vy Eq + P eq
		nNeu 	= 2*Grid->nyC;
		BC->n 	= nDir + nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n - 4;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;
	case SimpleShearPeriodic:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		nDir 	= 2*Grid->nxC; // Vx eq + Vy Eq + P eq
		// no Neumann nodes for this setup
		nNeu = 0;
		BC->n 	= nDir + nNeu;

		int nPeriod = 2*(Grid->nyC);

		EqSystem->nEq = EqSystem->nEqIni - nDir - nPeriod - 4; // the -4 corresponds to the corners
		break;
	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;

	}

	EqSystem->nRow = EqSystem->nEq;


	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*) malloc ( BC->n * sizeof(BCType));
	BC_updateThermal(BC, Grid);




	// Check
	// =======================================

	if (DEBUG) {
		printf("=== BC list thermal ===\n");
		printListi(BC->list,BC->n);
		printf("=== BC value thermal ===\n");
		printListd(BC->value,BC->n);
		printf("=== BC type thermal ===\n");
		int i;
		for (i=0;i<BC->n;i++) {
			printf("%d  ", BC->type[i]);
		}
		printf("\n");
	}

}





void BC_updateStokes(BC* BC, Grid* Grid)
{
	int C, i;
	int I = 0;

	if (BC->SetupType==PureShear) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		C = 0;
		for (i=0; i<Grid->nyVx; i++) { // Vx Left
			BC->list[I] = C;

			BC->value[I] = VxL;
			BC->type[I] = Dirichlet;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			BC->list[I] = C;
			BC->value[I] = VxR;
			BC->type[I] = Dirichlet;

			I++;
			C += Grid->nxVx;
		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			BC->list[I] = C;

			BC->value[I] = VyB;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			BC->list[I] = C;
			BC->value[I] = VyT;
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















	} else if (BC->SetupType==SimpleShearPeriodic) {
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		compute VxB =  0*1.0*BC->backStrainRate*Grid->ymin;
		compute VxT =  2*1.0*BC->backStrainRate*Grid->ymax;
		compute VyB =  0;
		compute VyT =  0;

		// Top and bottom Vy
		// =======================================
		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			BC->list[I] = C;
			BC->value[I] = VyB;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			BC->list[I] = C;
			BC->value[I] = VyT;
			BC->type[I] = Dirichlet;
			I++;
			C += 1;
		}

		// Top and bottom Vx
		// =======================================
		C = 0;
		for (i=0; i<Grid->nxVx; i++) { // Vx Bottom
			BC->list[I] = C;
			BC->value[I] = VxB; // factor 2 because it's a dirichlet condition on a ghost node
			BC->type[I] = DirichletGhost;
			I++;
			C += 1;
		}


		C = Grid->nxVx*(Grid->nyVx-1);
		for (i=0; i<Grid->nxVx; i++) { // Vx Top
			BC->type[I] = DirichletGhost;
			BC->list[I] = C;
			BC->value[I] = VxT;
			I++;
			C += 1;
		}




		// no Neumann nodes for this setup

	}

	else {
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);

	}
}






void BC_updateThermal(BC* BC, Grid* Grid)
{
	int C, i;
	int I = 0;

	if (BC->SetupType==PureShear) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute TT =  0.0;
		compute TB =  1.0;

		C = 1; // the first element in the numbering map is a ghost (in the sense of empty, i.e. there are no nodes in the corners)
		for (i=0; i<Grid->nxC; i++) { // Vx Left
			BC->list[I] = C;

			BC->value[I] = TB;
			BC->type[I] = Dirichlet;

			I++;
			C += 1;
		}


		C = (Grid->nxC+2)*(Grid->nyC+2-1)+1;
		for (i=0; i<Grid->nxC; i++) { // Vx Right
			BC->list[I] = C;
			BC->value[I] = TT;
			BC->type[I] = Dirichlet;

			I++;
			C += 1;
		}


		// Neumann
		// =======================================
		C = Grid->nxC+2;
		for (i=0;i<Grid->nyC;i++){ // Vy Left
			BC->list[I]          = C;
			BC->value[I]         = 0.0;
			BC->type[I] 		 = NeumannGhost;
			I++;
			C += Grid->nxC+2;
		}

		C = 2*(Grid->nxC+2)-1;
		for (i=0;i<Grid->nyC;i++){ // Vy Left
			BC->list[I]          = C;
			BC->value[I]         = 0.0;
			BC->type[I] 		 = NeumannGhost;
			I++;
			C += Grid->nxC+2;
		}








	} else if (BC->SetupType==SimpleShearPeriodic) {
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================

		compute TT =  0.0;
		compute TB =  1.0;

		C = 1; // the first element in the numbering map is a ghost (in the sense of empty, i.e. there are no nodes in the corners)
		for (i=0; i<Grid->nxC; i++) { // Vx Left
			BC->list[I] = C;

			BC->value[I] = TB;
			BC->type[I] = Dirichlet;

			I++;
			C += 1;
		}


		C = (Grid->nxC+2)*(Grid->nyC+2-1)+1;
		for (i=0; i<Grid->nxC; i++) { // Vx Right
			BC->list[I] = C;
			BC->value[I] = TT;
			BC->type[I] = Dirichlet;

			I++;
			C += 1;
		}



		// no Neumann nodes for this setup

	}

	else {
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);

	}
}



