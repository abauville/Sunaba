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
void BC_freeMemory(BC* BC) {
	free(BC->list);
	free(BC->value);
	free(BC->type);
}


void BC_initStokes(BC* BC, Grid* Grid, Physics* Physics, EqSystem* EqSystem)
{

	if (!DARCY) {
		BC->counter = 0;
		int nP, nV;
		BC_updateStokes_Vel(BC, Grid, Physics, false);
		nV = BC->counter;
		BC_updateStokes_P(BC, Grid, Physics, false);
		nP = BC->counter-nV;



		BC->n = BC->counter;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		if (BC->SetupType==SimpleShearPeriodic) {
			EqSystem->nEq -= Grid->nyVx-2 + 2*(Grid->nyVy-2) + 2*(Grid->nyEC-2); // because of periodic nodes
		}

		printf("BC->n = %i, nEq = %i, nEqIni = %i\n",BC->n, EqSystem->nEq, EqSystem->nEqIni);
		EqSystem->nRow = EqSystem->nEq;


		if (UPPER_TRI) {
			if (BC->SetupType==SimpleShearPeriodic) {
				EqSystem->nRow = EqSystem->nEq + nP  + 2*(Grid->nyEC-2) - Grid->nECTot;
			} else {
				EqSystem->nRow = EqSystem->nEq + nP  - Grid->nECTot;
			}
		}



	} else if (DARCY) {

		BC->counter = 0;
		BC_updateStokes_Vel(BC, Grid, Physics, false);
		BC_updateStokesDarcy_P(BC, Grid, Physics, false);
		BC->n = BC->counter;
		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		EqSystem->nRow = EqSystem->nEq;
	}





	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*)  malloc( BC->n * sizeof(BCType ));

	BC->counter = 0;

	if (!DARCY) {
		BC_updateStokes_Vel(BC, Grid, Physics, true);
		BC_updateStokes_P(BC, Grid, Physics, true);
	} else if (DARCY) {
		BC_updateStokes_Vel(BC, Grid, Physics, true);
		BC_updateStokesDarcy_P(BC, Grid, Physics, true);
	}






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


void BC_initThermal(BC* BC, Grid* Grid, Physics* Physics, EqSystem* EqSystem)
{
	/*
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
		nDir 	= 2*Grid->nxEC; // Vx eq + Vy Eq + P eq
		nNeu 	= 2*(Grid->nyEC-2);
		BC->n 	= nDir + nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;

	case SimpleShearPeriodic:
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================
		nDir 	= 2*Grid->nxEC; // Vx eq + Vy Eq + P eq
		// no Neumann nodes for this setup
		nNeu = 0;
		BC->n 	= nDir + nNeu;

		int nPeriod = 2*(Grid->nyC);

		EqSystem->nEq = EqSystem->nEqIni - nDir - nPeriod; // the -4 corresponds to the corners
		break;
	case FixedLeftWall:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		// Dirichlet on upper and lower
		// Neumann on the sides
		nDir 	= 2*Grid->nxEC; // Vx eq + Vy Eq + P eq
		nNeu 	= 2*(Grid->nyEC-2);
		BC->n 	= nDir + nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;
	case Sandbox:
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================
		// Dirichlet on upper and lower
		// Neumann on the sides
		nDir 	= 2*Grid->nxEC; // Vx eq + Vy Eq + P eq
		nNeu 	= 2*(Grid->nyEC-2);
		BC->n 	= nDir + nNeu;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		printf("### nEq = %i\n", EqSystem->nEq);
		break;
	default:
		printf("Unknown BC.SetupType %i", BC->SetupType);
		exit(0);
		break;

	}
	 */

	BC->counter = 0;
	BC_updateThermal(BC, Grid, Physics, false);
	BC->n = BC->counter;
	EqSystem->nEq = EqSystem->nEqIni - BC->n;
	printf("EqSystem->nEq = %i\n",EqSystem->nEq);
	if (BC->SetupType==SimpleShearPeriodic) {
		EqSystem->nEq -= 2*(Grid->nyEC-2); // because of periodic nodes
	}
	EqSystem->nRow = EqSystem->nEq;
	printf("EqSystem->nEq = %i\n",EqSystem->nEq);

	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*) malloc ( BC->n * sizeof(BCType));

	BC->counter = 0;
	BC_updateThermal(BC, Grid, Physics, true);




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





void BC_updateStokes_Vel(BC* BC, Grid* Grid, Physics* Physics, bool assigning)
{
	// if assigining == true BC->val, Bc->list and BC->type are filled
	// otherwise the number of BC are just counted

	int C, i;

	int I = BC->counter;




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
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VxL;
				BC->type[I] = Dirichlet;
				C += Grid->nxVx;
			}
			I++;

		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VxR;
				BC->type[I] = Dirichlet;
				C += Grid->nxVx;
			}

			I++;

		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VyB;
				BC->type[I] = Dirichlet;
				C += 1;
			}
			I++;

		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);

		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VyT;
				BC->type[I] = Dirichlet;
				C += 1;
			}
			I++;

		}




		// Neumann
		// =======================================


		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         =  0.0;
				BC->type[I] 		 = NeumannGhost;
				C = C+Grid->nxVy;
			}
			I++;

		}




		C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
				C = C+Grid->nxVy;
			}
			I++;

		}

		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
			if (assigning) {

				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] = NeumannGhost;
				C = C+1;
			}
			I++;

		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]         = 0.0;
				BC->type[I] = NeumannGhost;
				C = C+1;
			}
			I++;

		}
	}



	else if (BC->SetupType==SimpleShearPeriodic) {
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================1
		compute VxB =  2.0*1.0*BC->backStrainRate*Grid->ymin;
		compute VxT =  2.0*1.0*BC->backStrainRate*Grid->ymax;
		compute VyB =  0;
		compute VyT =  0;

		// Top and bottom Vy
		// =======================================
		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VyB;
				BC->type[I] = Dirichlet;
				C += 1;
			}
			I++;

		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VyT;
				BC->type[I] = Dirichlet;
				C += 1;
			}
			I++;
		}

		// Top and bottom Vx
		// =======================================
		C = 0;
		for (i=0; i<Grid->nxVx; i++) { // Vx Bottom
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VxB; // factor 2 because it's a dirichlet condition on a ghost node
				BC->type[I] = DirichletGhost;
				C += 1;
			}
			I++;

		}


		C = Grid->nxVx*(Grid->nyVx-1);
		for (i=0; i<Grid->nxVx; i++) { // Vx Top
			if (assigning) {
				BC->type[I] = DirichletGhost;
				BC->list[I] = C;
				BC->value[I] = VxT;
				C += 1;
			}
			I++;

		}

	}


	else if (BC->SetupType==FixedLeftWall) {
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
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VxL;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += Grid->nxVx;
		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VxR;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += Grid->nxVx;
		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VyB;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VyT;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += 1;
		}







		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = DirichletGhost;
			}
			I++;
			C = C+Grid->nxVy;
		}




		// Neumann
		// =======================================


		C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
			}
			I++;
			C = C+Grid->nxVy;
		}

		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] = NeumannGhost;
			}
			I++;
			C = C+1;
		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]         = 0.0;
				BC->type[I] = NeumannGhost;
			}
			I++;
			C = C+1;
		}

	}
	else if (BC->SetupType==Sandbox) {
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
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VxL;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += Grid->nxVx;
		}




		C = 2*Grid->nxVx-1;
		for (i=0; i<Grid->nyVx-1; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;


				BC->value[I] = VxR;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += Grid->nxVx;
		}




		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VyB;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += 1;
		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VyT;
				BC->type[I] = Dirichlet;
			}
			I++;
			C += 1;
		}



		C = 1;
		for (i=0;i<Grid->nxVx-1;i++){ // Vx Bottom
			if (assigning) {
				BC->list[I]  = C;
				BC->value[I] = VxL;//(VxL+VxR)/2.0;
				BC->type[I]  = DirichletGhost;
			}
			I++;
			C = C+1;
		}




		// Neumann
		// =======================================


		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
			}
			I++;
			C = C+Grid->nxVy;
		}




		//C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		C = Grid->nVxTot + 2*Grid->nxVy-1;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
			}
			I++;
			C = C+Grid->nxVy;
		}



		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]        = 0.0;
				BC->type[I] 		= NeumannGhost;
			}
			I++;
			C = C+1;
		}

	}

	else {
		printf("Unknown Stokes BC.SetupType %i", BC->SetupType);
		exit(0);

	}



	BC->counter = I;


}



void BC_updateStokes_P(BC* BC, Grid* Grid, Physics* Physics, bool assigning)
{
	int C, I, i;

	I = BC->counter;

	// Pressure BC for all setup
	// =======================================
	// in normal stokes there is lagrangian operator on the Pressure, only the gradient
	// and it doesn't have to be build for boundary velocity nodes (since they are BC)
	// therefore these are just dummy values in this case
	// However the following acts as Pf Boundary conditions for the two-phase flow


	C = Grid->nVxTot + Grid->nVyTot;
	for (i=0;i<Grid->nxEC;i++){ // PBottom
		if (assigning) {
			BC->list[I]         = C;
			BC->value[I]        = 0.0;
			BC->type[I] = NeumannGhost;
		}
		I++;
		C = C+1;
	}


	C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-1)*Grid->nxEC;
	for (i=0;i<Grid->nxEC;i++){ // PTop
		if (assigning) {
			BC->list[I]         = C;
			BC->value[I]        = 0.0;
			BC->type[I] = NeumannGhost;
		}
		I++;
		C = C+1;
	}





	if (BC->SetupType!=SimpleShearPeriodic) {
		C =  Grid->nVxTot + Grid->nVyTot + Grid->nxEC;
		for (i=0;i<Grid->nyEC-2;i++){ // Pleft
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]        = 0.0;
				BC->type[I] = NeumannGhost;
			}
			I++;
			C = C+Grid->nxEC;
		}

		C = Grid->nVxTot + Grid->nVyTot + Grid->nxEC-1 + Grid->nxEC;
		for (i=0;i<Grid->nyEC-2;i++){ // Prigth
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]        = 0.0;
				BC->type[I] = NeumannGhost;
			}
			I++;
			C = C+Grid->nxEC;
		}
	}

	BC->counter = I;


}


#if (DARCY)
void BC_updateStokesDarcy_P(BC* BC, Grid* Grid, Physics* Physics, bool assigning) {
	int C, I, i, iP;
	int NumberMod;

	I = BC->counter;


	// Pressure BC for all setup
	// =======================================
	// in normal stokes there is lagrangian operator on the Pressure, only the gradient
	// and it doesn't have to be build for boundary velocity nodes (since they are BC)
	// therefore these are just dummy values in this case
	// However the following acts as Pf Boundary conditions for the two-phase flow

	for (iP = 0; iP < 2; ++iP) {
		if (iP==0) {
			NumberMod = 0;
		} else if (iP == 1) {
			NumberMod = Grid->nECTot;
			//iP = 2;// used only for debugging, to set all Pc to 0
		}


		if (iP==1) { // Pc, i.e. Dummy
			C = Grid->nVxTot + Grid->nVyTot + NumberMod;
			for (i=0;i<Grid->nxEC;i++){ // PBottom
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;//-1.0*Physics->rho[i]*Physics->g[1];
					BC->type[I] = NeumannGhost;
				}
				I++;
				C = C+1;
			}
		} else if (iP == 0 ) {
			C = Grid->nVxTot + Grid->nVyTot + NumberMod;
			for (i=0;i<Grid->nxEC;i++){ // PBottom
				if (assigning) {
					BC->list[I]         = C;
					//BC->value[I]        = 0.0;//
					BC->value[I]        = 1.0*Physics->rho_f*Physics->g[1];//1.0*Physics->rho[i]*Physics->g[1];
					//BC->value[I]        = 1.0*Physics->rho[i]*Physics->g[1];
					//BC->value[I]        = 1.0*Physics->rho_f*Physics->g[1] + 0.5*(1.0*Physics->rho[i]*Physics->g[1]+1.0*Physics->rho_f*Physics->g[1]);
					BC->type[I] 		= NeumannGhost;
				}
				I++;
				C = C+1;
			}
		}





		if (iP==1) { // Pc, i.e. Dummy
			C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-1)*Grid->nxEC + NumberMod;
			for (i=0;i<Grid->nxEC;i++){ // PTop
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					BC->type[I] 		= Dirichlet;
				}
				I++;
				C = C+1;
			}
		} else if (iP == 0 ) {
			C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-1)*Grid->nxEC + NumberMod;
			for (i=0;i<Grid->nxEC;i++){ // PTop
				if (assigning) {
					BC->list[I]         = C;
					if (Physics->y_oceanSurface < 0.0 + 1e-8) {
						// 0.0 is the default value
						BC->value[I]        = 0.0;
					} else {
						BC->value[I] =  Physics->rho_f*Physics->g[1]*(Grid->ymax-(Physics->y_oceanSurface+Grid->ymin));
					}
					BC->type[I] 		= DirichletGhost;
				}
				I++;
				C = C+1;
			}
		}






		if (iP==1) { // Pc, i.e. Dummy
			if (BC->SetupType!=SimpleShearPeriodic) {
						C =  Grid->nVxTot + Grid->nVyTot + Grid->nxEC + NumberMod;
						for (i=0;i<Grid->nyEC-2;i++){ // PLeft
							if (assigning) {
								BC->list[I]         = C;
								BC->value[I]        = 0.0;
								BC->type[I] 		= Dirichlet;
							}
							I++;
							C = C+Grid->nxEC;
						}

						C = Grid->nVxTot + Grid->nVyTot + Grid->nxEC-1 + Grid->nxEC + NumberMod;
						for (i=0;i<Grid->nyEC-2;i++){ // PRight
							if (assigning) {
								BC->list[I]         = C;
								BC->value[I]        = 0.0;
								BC->type[I] 		= Dirichlet;
							}
							I++;
							C = C+Grid->nxEC;
						}
					}
		} else if (iP == 0 ) {
			if (BC->SetupType!=SimpleShearPeriodic) {
						C =  Grid->nVxTot + Grid->nVyTot + Grid->nxEC + NumberMod;
						for (i=0;i<Grid->nyEC-2;i++){ // PLeft
							if (assigning) {
								BC->list[I]         = C;
								BC->value[I]        = 0.0;
								BC->type[I] 		= NeumannGhost;
							}
							I++;
							C = C+Grid->nxEC;
						}

						C = Grid->nVxTot + Grid->nVyTot + Grid->nxEC-1 + Grid->nxEC + NumberMod;
						for (i=0;i<Grid->nyEC-2;i++){ // PRight
							if (assigning) {
								BC->list[I]         = C;
								BC->value[I]        = 0.0;
								BC->type[I] 		= NeumannGhost;
							}
							I++;
							C = C+Grid->nxEC;
						}
					}
		}





		// Second row from the top, set Pc to 0
		if (iP==1) {
			/*
			C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-2)*Grid->nxEC + 1 + NumberMod;
			for (i=0;i<Grid->nxEC-2;i++){ // PTop
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					//BC->value[I]        = Physics->rho[i + (Grid->nyEC-2)*Grid->nxEC + 1] * fabs(Physics->g[1]) * Grid->dy/2  - Physics->Pf[i + (Grid->nyEC-2)*Grid->nxEC + 1];
					BC->type[I] 		= Dirichlet;
				}
				I++;
				C = C+1;
			}
			*/




		} else if (iP == 0 ) {
			/*
			C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-2)*Grid->nxEC + 1 + NumberMod;
			for (i=0;i<Grid->nxEC-2;i++){ // PTop
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					BC->type[I] 		= Dirichlet;
				}
				I++;
				C = C+1;
			}
			*/


		}




		if (iP==2) {
			for (i = 0; i < Grid->nECTot; ++i) {
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					BC->type[I] 		= Dirichlet;
				}
				I++;
				C = C+1;
			}
		}





	}


















	BC->counter = I;

}


#endif

void BC_updateThermal(BC* BC, Grid* Grid, Physics* Physics, bool assigning)
{
	int C, i;
	int I = BC->counter;

	if (BC->SetupType==PureShear || BC->SetupType==FixedLeftWall || BC->SetupType==Sandbox) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================


		C = 0; // the first element in the numbering map is a ghost (in the sense of empty, i.e. there are no nodes in the corners)
		for (i=0; i<Grid->nxEC; i++) { // Bottom boundary
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = BC->TB;
				BC->type[I] = DirichletGhost;
			}
			I++;
			C += 1;
		}


		C = (Grid->nxEC)*(Grid->nyEC-1);
		for (i=0; i<Grid->nxEC; i++) { // Top boundary
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = BC->TT;
				BC->type[I] = DirichletGhost;
			}
			I++;
			C += 1;
		}




		// Neumann
		// =======================================
		C = Grid->nxEC;
		for (i=1;i<Grid->nyEC-1;i++){ // Left boundary
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
			}
			I++;
			C += Grid->nxEC;
		}

		C = 2*(Grid->nxEC)-1;
		for (i=1;i<Grid->nyEC-1;i++){ // Right boundary
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;
			}
			I++;
			C += Grid->nxEC;
		}









	} else if (BC->SetupType==SimpleShearPeriodic) {
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================


		C = 0; // the first element in the numbering map is a ghost (in the sense of empty, i.e. there are no nodes in the corners)
		for (i=0; i<Grid->nxEC; i++) { // Bottom
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = BC->TB;
				BC->type[I] = DirichletGhost;
			}
			I++;
			C += 1;
		}


		C = (Grid->nxEC)*(Grid->nyEC-1);
		for (i=0; i<Grid->nxEC; i++) { // Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = BC->TT;
				BC->type[I] = DirichletGhost;
			}
			I++;
			C += 1;
		}





		// no Neumann nodes for this setup

	}


	else {
		printf("Unknown Temp BC.SetupType %i", BC->SetupType);
		exit(0);

	}


	BC->counter = I;
}



