/*
 * BC.c
 *
 *  Created on: Feb 23, 2016
 *      Author: abauville
 */

#include "stokes.h"
#define SQUARE(x) (x)*(x)


// Functions for the Corner Flow BC
static inline compute VxArc(compute alpha, compute U, compute x, compute y)
{
	compute atanYX = atan(y/x);
	if (atanYX<0) {
		//atanYX += PI;
	}
	return -2*U*(x*(x*sin(alpha)*atan(tan(alpha)) + y*(sin(alpha) - cos(alpha)*atan(tan(alpha)))) + (x*x + y*y)*((sin(alpha) - cos(alpha)*atan(tan(alpha)))*atanYX - sin(alpha)*atan(tan(alpha))))/((x*x + y*y)*(cos(2*alpha) + 2*SQUARE(atan(tan(alpha))) - 1));
}

static inline compute VyArc(compute alpha, compute U, compute x, compute y)
{
	compute atanYX = atan(y/x);
	if (atanYX<0) {
		//atanYX += PI;
	}
	return 2*U*(-y*(x*sin(alpha)*atan(tan(alpha)) + y*(sin(alpha) - cos(alpha)*atan(tan(alpha)))) + (x*x + y*y)*sin(alpha)*atanYX*atan(tan(alpha)))/((x*x + y*y)*(cos(2*alpha) + 2*SQUARE(atan(tan(alpha))) - 1));

}
static inline compute VxOcean(compute alpha, compute U, compute x, compute y)
{
	compute atanYX = atan(y/x);
	if (x<0) {
		atanYX += PI;
	}
	return U*(x*(x*sin(alpha) - y*cos(alpha) + y) + (x*x + y*y)*(-(cos(alpha) - 1)*atanYX + PI*cos(alpha) - atan(tan(alpha))))/((x*x + y*y)*(sin(alpha) - atan(tan(alpha)) + PI));
	//return U*(x*(x*sin(alpha) - y*cos(alpha) + y) - (x*x + y*y)*(atan(tan(alpha)) - PI))/((x*x + y*y)*(sin(alpha) - atan(tan(alpha)) + PI));

	//return U*(x*(x*sin(alpha) - y*cos(alpha) + y) - (x*x + y*y)*((cos(alpha) - 1)*atan(y/x) + atan(tan(alpha))))/((x*x + y*y)*(sin(alpha) - atan(tan(alpha))));
}
static inline compute VyOcean(compute alpha, compute U, compute x, compute y)
{
	compute atanYX = atan(y/x);
	if (x<0) {
		atanYX += PI;
	}

	return U*(-y*(-x*sin(alpha) + y*(cos(alpha) - 1)) + (x*x + y*y)*(-atanYX + PI)*sin(alpha))/((x*x + y*y)*(sin(alpha) - atan(tan(alpha)) + PI));
	//return U*y*(x*sin(alpha) - y*cos(alpha) + y)/((x*x + y*y)*(sin(alpha) - atan(tan(alpha)) + PI));
	//return -U*(-y*(x*sin(alpha) - y*(cos(alpha) - 1)) + (x*x + y*y)*sin(alpha)*atan(y/x))/((x*x + y*y)*(sin(alpha) - atan(tan(alpha))));
}

compute CornerVelocity(Grid* Grid, compute alpha, compute U, compute x, compute y, bool type) {
	// give the Corner Flow velocity at the given point, returns Vx if type = 0, or Vy if type = 1
	compute r, xSlab, Value;
	//x = Grid->xmin + Grid->dx*ix;
	//y = 0.0*Grid->ymax - (Grid->ymin + Grid->dy*iy);

	y = -y;


	//printf("iy = %i, y = %.2e\n",iy, y);
	r = sqrt(x*x + y*y);
	xSlab = r*cos(alpha);
	//printf("ix = %i, iy = %i, x = %.2e, y = %.2e\n", ix, iy, x, y);
	if (x>xSlab) {
		//printf("use Arc\n");
		if (type==0) {
			Value = VxArc(alpha,U,x,y);
		} else {
			Value = - VyArc(alpha,U,x,y);
		}
	} else {
		//printf("use Ocean\n");
		if (type==0) {
			Value = VxOcean(alpha,U,x,y);
		} else {
			Value = - VyOcean(alpha,U,x,y);
		}
	}

	/*
	if (type==0) {
		Value = VxOcean(alpha,U,x,y);
	} else {
		Value = - VyOcean(alpha,U,x,y);
	}
	*/





	//printf("ix, = %i, iy = %i, Value = %.2e\n",ix, iy, Value);
	/*
	if (x==0) {
		printf("WARNING WARNING WARNING: Found x = 0, but the corner flow solution is not defined for x=0 at y=0 (because of division by (x^2+y^2)), if you apply Vx at the surface it will result  in NaN. If you want to do so change your mesh.\n");
		//exit(0);
	}
	*/

	return Value;
}


//==========================================================================
//
//                            BOUNDARY INDEXING
//
//==========================================================================
void BC_Memory_free(BC* BC) {
	free(BC->list);
	free(BC->value);
	free(BC->type);
}


void BC_initStokes(BC* BC, Grid* Grid, Physics* Physics, EqSystem* EqSystem)
{
		BC->reCompute_SymbolicFactorization = false;
		BC->iyTopRow = Grid->nyS;
	if (!DARCY) {
		BC->counter = 0;
		int nP, nV;
		BC_updateStokes_Vel(BC, Grid, Physics, false);
		nV = BC->counter;
		BC_updateStokes_P(BC, Grid, Physics, false);
		nP = BC->counter-nV;



		BC->n = BC->counter;

		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		if (Grid->isPeriodic) {
			EqSystem->nEq -= Grid->nyVx-2 + 2*(Grid->nyVy-2) + 2*(Grid->nyEC-2); // because of periodic nodes
		}

		printf("BC->n = %i, nEq = %i, nEqIni = %i\n",BC->n, EqSystem->nEq, EqSystem->nEqIni);
		EqSystem->nRow = EqSystem->nEq;


		if (UPPER_TRI) {
#if (PENALTY_METHOD)
				EqSystem->nRow = EqSystem->nEq;
#else
			if (Grid->isPeriodic) {
				EqSystem->nRow = EqSystem->nEq + nP  + 2*(Grid->nyEC-2) - Grid->nECTot;
			} else {
				EqSystem->nRow = EqSystem->nEq + nP  - Grid->nECTot;
			}
#endif

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
	if (Grid->isPeriodic) {
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


	BC->IsFreeSlipLeft	= false;
	BC->IsFreeSlipRight = false;
	BC->IsFreeSlipBot 	= false;
	BC->IsFreeSlipTop 	= false;

	if (BC->SetupType==Stokes_PureShear) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= true;
		BC->IsFreeSlipTop 	= true;



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



	else if (BC->SetupType==Stokes_SimpleShear) {
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


	else if (BC->SetupType==Stokes_FixedLeftWall) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		BC->IsFreeSlipLeft	= false;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= true;
		BC->IsFreeSlipTop 	= true;

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
	else if (BC->SetupType==Stokes_Sandbox) {
		// =======================================
		// =======================================
		// 				Sandbox
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		//compute outFlowH = (Grid->ymax-Grid->ymin)/5.0;
		compute integralOutflowVxdy = 0.0;
		compute extraOutFlowVy;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipBot 	= false;
		BC->IsFreeSlipTop 	= true;

		if (BC->Sandbox_NoSlipWall) {
			BC->IsFreeSlipRight = false;
		} else {
			BC->IsFreeSlipRight = true;
		}


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




		C = 1*Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;


				BC->value[I] = VxR;
				BC->type[I] = Dirichlet;






					 // OutFlow
				compute y00 = (BC->Sandbox_TopSeg00 - (Grid->ymin + (i) * Grid->dy))/BC->Sandbox_TopSeg00;
				compute y01 = (BC->Sandbox_TopSeg01 - (Grid->ymin + (i) * Grid->dy))/(BC->Sandbox_TopSeg01-BC->Sandbox_TopSeg00);
				//printf("y00 = %.2e, y01= %.2e y = %.2e,\n",y00, y01, (Grid->ymin + (i) * Grid->dy));
				if (y00>0.0) {
					//BC->value[I] = y00*VxL;
					BC->value[I] = VxL;
				}
				else if (y01>0.0) {
					BC->value[I] = y01*VxL;
				}


				if (y01>0.0) {

					if (i>0) {
						integralOutflowVxdy += BC->value[I]*Grid->dy;
					}

					//printf("y = %.2e, VxL = %.2e, BC->value[I] = %.2e\n",y, VxL, y*VxL);
				}



					//printf("y = %.2e, VxL = %.2e, BC->value[I] = %.2e\n",y, VxL, y*VxL);


			}
			I++;
			C += Grid->nxVx;
		}
		extraOutFlowVy = -integralOutflowVxdy/((Grid->nxVy-2)*Grid->dx);



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


		int iy, ix;
		int highestTopo_iy = 0;
		int lowestTopo_iy = Grid->nyS+2;
		//int nySedIni = (int) round(1.0/Grid->dy); // assumes that the non dimensional thickness of sediment is 1.0
		//int nyAirMin = (int) round(nySedIni/2.0);



		for(iy = 0;iy < Grid->nyEC;iy++) {
			for(ix = 0;ix < Grid->nxEC;ix++) {
				if (Physics->phase[ix+iy*Grid->nxEC]!=Physics->phaseAir && Physics->phase[ix+iy*Grid->nxEC]!=Physics->phaseWater ) {
					if (iy>highestTopo_iy) {
						highestTopo_iy = iy;
					}
				} else {
					if (iy<lowestTopo_iy) {
						lowestTopo_iy = iy;
					}
				}
			}
		}
		int nyAirMin = (int) round(.75*lowestTopo_iy);

		int	nTopRowsWithBC = 1;
		//BC->iyTopRow_tolerance = (int) round(0.15*lowestTopo_iy);
		BC->iyTopRow_tolerance = (int) round(100.0*lowestTopo_iy);
		int iyTopRow_ideal = highestTopo_iy + nyAirMin;
		
		
		
		
		  
		
		if (abs(iyTopRow_ideal-BC->iyTopRow)>BC->iyTopRow_tolerance && BC->iyTopRow!= Grid->nyS) {
			// then update iyTopRow
			BC->iyTopRow = highestTopo_iy + nyAirMin;
			BC->reCompute_SymbolicFactorization = true;
		} else {
			BC->reCompute_SymbolicFactorization = false;
		}
		nTopRowsWithBC = Grid->nyS - BC->iyTopRow+1;
		if (nTopRowsWithBC<1) {
			nTopRowsWithBC = 1;
			BC->iyTopRow = Grid->nyS;
		}
		
		printf("iyTopRow = %i, nTopRowsWithBC=%i, highestTopo_iy = %i, nyAirMin = %i\n", BC->iyTopRow, nTopRowsWithBC, highestTopo_iy, nyAirMin);
		//nTopRowsWithBC = 1;
		
		for (iy=BC->iyTopRow-1;iy<Grid->nyVy;iy++) {
			//C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1-iy);
			for (ix=0; ix<Grid->nxVy; ix++) { // Vy Top
				if (assigning) {
					BC->list[I] = Grid->nVxTot + ix+iy*Grid->nxVy;//C;
					compute y = (Grid->Y[iy]-Grid->dy-Grid->ymin)/(Grid->ymax-Grid->ymin);

					BC->value[I] = (1.0-y)*VyB + y*VyT + extraOutFlowVy;
					BC->type[I] = Dirichlet;

					//BC->value[I] = 0.0;
					//BC->type[I] = Neumann;

				}
				I++;
				//C += 1;
			}
		}



		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
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
		for (i=0;i<Grid->nyVy-2-(nTopRowsWithBC-1);i++){ // Vy Left
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
		for (i=0;i<Grid->nyVy-2-(nTopRowsWithBC-1);i++){ // Vy Right
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;

				//y = Grid->ymin + Grid->dy*(i+1);
				//if (y<Grid->ymin+(Grid->ymax-Grid->ymin)/12.0) {
				//y = (outFlowH - (Grid->ymin + (i) * Grid->dy))/outFlowH;
				//if (i<Grid->nyVy-10) {
					if (BC->Sandbox_NoSlipWall) {
						if (Physics->phase[Grid->nxEC-1 + i * Grid->nxEC] != Physics->phaseAir && Physics->phase[Grid->nxEC-1 + i * Grid->nxEC] != Physics->phaseWater) {
							BC->type[I] 		 = DirichletGhost;
						}
					}
				//}

				//BC->type[I] 		 = DirichletGhost;
			}
			I++;
			C = C+Grid->nxVy;
		}


		
		for (iy=BC->iyTopRow;iy<Grid->nyVx;iy++) {
			//C = Grid->nxVx*(Grid->nyVx-1-iy)+1;
			for (ix=1;ix<Grid->nxVx-1;ix++){ // Vx Top
				if (assigning) {
					BC->list[I]         = ix+iy*Grid->nxVx;//C;
					compute x = (Grid->X[ix]-Grid->dx-Grid->xmin)/(Grid->xmax-Grid->xmin);
					BC->value[I]        = (1.0-x)*VxL + (x)*VxR;
					BC->type[I] 		= Dirichlet;
					
					//BC->value[I]        = 0.0;
					//BC->type[I] 		= NeumannGhost;
					
				}
				I++;
				//C = C+1;
			}
		}
		
		
		/*
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
		*/
		

	} else if (BC->SetupType==Stokes_SandboxWeakBackstop) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= false;
		BC->IsFreeSlipTop 	= true;


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
				if (Physics->phase[i + 0*Grid->nxC] == BC->specialPhase) { // if special phase => dragging down
					BC->list[I]  = C;
					BC->value[I] = VxL;//(VxL+VxR)/2.0;
					BC->type[I]  = DirichletGhost;
				} else { 												   // else free slip
					BC->list[I]          = C;
					BC->value[I]         = 0.0;
					BC->type[I] = NeumannGhost;
				}
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

	} else if (BC->SetupType==Stokes_Sandbox_InternalBC) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyB = -BC->backStrainRate*Grid->ymin;
		compute VyT = -BC->backStrainRate*Grid->ymax;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= false;
		BC->IsFreeSlipTop 	= true;


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
				if (i<Grid->nxVx/2) { // if special phase => dragging down
					BC->list[I]  = C;
					BC->value[I] = VxL;//(VxL+VxR)/2.0;
					BC->type[I]  = DirichletGhost;
				} else { 												   // stick slip
					BC->list[I]          = C;
					BC->value[I]         = 0.0;
					BC->type[I] = DirichletGhost;
				}
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

	else if (BC->SetupType==Stokes_CornerFlow) {
		// =======================================
		// =======================================
		// 				Corner Flow
		// =======================================
		// =======================================
		compute VxL =  BC->backStrainRate*Grid->xmin;
		compute VxR =  BC->backStrainRate*Grid->xmax;
		compute VyT = -BC->backStrainRate*Grid->ymax;
		printf("Ini, VyT = %.2e, Alt = %.2e, VxR = %.2e, VxL = %.2e\n",VyT, -Grid->ymax*0.5*(VxR/Grid->xmax + VxL/Grid->xmin), VxR, VxL);
		//int* CornerBCType = (int*) malloc(BC->n * sizeof(int)); // 0: Vx Arc, 1: Vy Arc, 2: Vx Ocean, 3: VyOcean

		int ix, iy;
		compute alpha = BC->Corner_SubductionAngle;//PI/4;

		compute U = BC->refValue;
		compute y,x;
		compute ySurf = 0.0;
		C = 0;
		printf("VxLeft\n");

		BC->IsFreeSlipLeft	= false;
		BC->IsFreeSlipRight = false;
		BC->IsFreeSlipBot 	= false;
		BC->IsFreeSlipTop 	= true;


		for (i=0; i<Grid->nyVx; i++) { // Vx Left
			if (assigning) {
				BC->list[I] = C;
				BC->type[I] = Dirichlet;

				ix = 0;
				//iy = i;
				x = (Grid->xmin + Grid->dx*ix);
				y = (Grid->ymin + Grid->dy*i) - 0.5*Grid->dy;

				if (y<=ySurf) {
					BC->value[I] = CornerVelocity(Grid, alpha, U, x, y, 0);
				} else {
					BC->value[I] = CornerVelocity(Grid, alpha, U, x, ySurf, 0);
					VxL = BC->value[I];
				}

				C += Grid->nxVx;

			}
			I++;

		}


		C = Grid->nxVx-1;
		printf("VxRight\n");
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;
				BC->type[I] = Dirichlet;

				ix = Grid->nxS-1;
				x = (Grid->xmin + Grid->dx*ix);
				y = (Grid->ymin + Grid->dy*i) - 0.5*Grid->dy;

				if (y<=ySurf) {
					BC->value[I] = CornerVelocity(Grid, alpha, U, x, y, 0);
				} else {
					BC->value[I] = CornerVelocity(Grid, alpha, U, x, ySurf, 0);
					VxR = BC->value[I];
				}


				C += Grid->nxVx;
			}

			I++;

		}


		printf("VyBottom\n");
		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;
				BC->type[I] = Dirichlet;

				ix = i;
				iy = 0;
				x = (Grid->xmin + Grid->dx*i) - 0.5*Grid->dx;
				y = (Grid->ymin + Grid->dy*iy);
				BC->value[I] = CornerVelocity(Grid, alpha, U, x, y, 1);




				C += 1;
			}
			I++;

		}

		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);
		printf("VyTop\n");
		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->type[I] = Dirichlet;

				ix = i;
				iy = Grid->nyS-1; // Boundary are defined on the shear nodes


				BC->value[I] = Grid->ymax*((VxL-VxR)/(Grid->xmax-Grid->xmin));//equivalent to VyT/2.0, but I'm not sure that the right hand side is always at VxR = 0, so I leave it like that
				//printf("VyT = %.2e, Alt = %.2e, VxR = %.2e, VxL = %.2e\n",VyT, Grid->ymax*((VxL-VxR)/(Grid->xmax-Grid->xmin)), VxR, VxL);

				C += 1;
			}
			I++;

		}



		// Ghost nodes
		// =======================================

		printf("VyLeft\n");
		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			if (assigning) {
				BC->list[I]         = C;

				ix = 0;
				//iy = i+1; // Boundary are defined on the shear nodes

				y = (Grid->ymin + Grid->dy*i);
				//if (y<=ySurf) { // it will stop updating iy in the sticky air, so that the stickyair has the surface velocity
					iy = i+1;
				//}

				//BC->value[I] 		= CornerVelocity(Grid, alpha, U, ix, iy, 1);
				//BC->type[I] 		= DirichletGhost;


				BC->value[I] 		= 0.0;
				BC->type[I] 		= NeumannGhost;

				C = C+Grid->nxVy;
			}
			I++;

		}




		C = Grid->nVxTot + Grid->nxVy-1 + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Right
			if (assigning) {
				BC->list[I]         = C;

				ix = Grid->nxS-1;
				//iy = i+1; // Boundaries are defined on the shear nodes
				y = (Grid->ymin + Grid->dy*i);
				//if (y<=ySurf) { // it will stop updating iy in the sticky air, so that the stickyair has the surface velocity
					iy = i+1;

				x = Grid->xmin + Grid->dx*ix  - 0.5*Grid->dx;
				y = (Grid->ymin + Grid->dy*iy);
				BC->value[I] 		= CornerVelocity(Grid, alpha, U, x, y, 0);
				BC->type[I] 		= DirichletGhost;


				BC->value[I] 		= 0.0;
				BC->type[I] 		= NeumannGhost;

				C = C+Grid->nxVy;
			}
			I++;

		}

		C = 1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Bottom
			if (assigning) {

				BC->list[I]         = C;

				ix = i+1;
				iy = 0; // Boundary are defined on the shear nodes
				x = Grid->xmin + Grid->dx*ix;
				y = (Grid->ymin + Grid->dy*iy);
				BC->value[I] 		= CornerVelocity(Grid, alpha, U, x, y, 0);
				BC->type[I] 		= DirichletGhost;

				//BC->value[I] 		= 0.0;
				//BC->type[I] 		= NeumannGhost;

				C = C+1;
			}
			I++;

		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top
			if (assigning) {
				BC->list[I]         = C;
				ix = i+1;
				iy = Grid->nyS-1; // Boundary are defined on the shear nodes

				//BC->value[I] 		= CornerVelocity(Grid, alpha, U, ix, iy, 0);
				//BC->type[I] 		= DirichletGhost;

				BC->value[I] 		= 0.0;
				BC->type[I] 		= NeumannGhost;

				C = C+1;
			}
			I++;

		}


	}

	else if (BC->SetupType==Stokes_WindTunnel) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================

		//compute VxL =  BC->backStrainRate*Grid->xmin;
		//compute VxR =  BC->backStrainRate*Grid->xmax;
		//compute VyB = -BC->backStrainRate*Grid->ymin;
		//compute VyT = -BC->backStrainRate*Grid->ymax;

		compute VxInput = BC->refValue;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= true;
		BC->IsFreeSlipTop 	= true;





		C = 0;
		for (i=0; i<Grid->nyVx; i++) { // Vx Left
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = VxInput;
				BC->type[I] = Dirichlet;

				C += Grid->nxVx;
			}
			I++;

		}


		C = Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = 0.0;
				BC->type[I] = Neumann;

				C += Grid->nxVx;
			}

			I++;

		}


		C = Grid->nVxTot + 0;
		for (i=0; i<Grid->nxVy; i++) { // Vy Bottom
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = 0.0;
				BC->type[I] = Dirichlet;
				C += 1;
			}
			I++;

		}


		C = Grid->nVxTot + Grid->nxVy*(Grid->nyVy-1);

		for (i=0; i<Grid->nxVy; i++) { // Vy Top
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = 0.0;
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
				BC->value[I]         = 0.0;
				BC->type[I] 		 = Dirichlet;
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
				BC->type[I] 		 = NeumannGhost;
				C = C+1;
			}
			I++;

		}

		C = Grid->nxVx*(Grid->nyVx-1)+1;
		for (i=0;i<Grid->nxVx-2;i++){ // Vx Top
			if (assigning) {
				BC->list[I]         = C;
				BC->value[I]        = 0.0;
				BC->type[I] 		= NeumannGhost;
				C = C+1;
			}
			I++;

		}

		// Inner boundary condtions
		compute radius = (Grid->ymax-Grid->ymin)/11.0;
		compute cx = 0.0;
		compute cy = 0.0;
		int ix, iy;
		compute x, y;
		// Vx inner BC
		for (iy = 0; iy < Grid->nyVx; ++iy) {
			for (ix = 0; ix < Grid->nxVx; ++ix) {
				x = ix*Grid->dx+ Grid->xmin;
				y = iy*Grid->dy + Grid->ymin - 0.5*Grid->dy;
				//printf("DXS[%i] = %.2e\n",ix,Grid->DXS[ix]);
				if ( sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) < radius ) {
					//printf("x = %.2e, cx = %.2e, y = %.2e, cy = %.2e, ix = %i, iy = %i, y2 = %.2e\n", x, cx, y, cy, ix, iy, iy*Grid->dy + Grid->ymin);

					if (assigning) {
					BC->list[I] = ix+iy*Grid->nxVx;
					BC->value[I] = 0.0;
					BC->type[I] = Dirichlet;
					}
					I++;


				}
			}
		}
		// Vy inner BC
		for (iy = 0; iy < Grid->nyVy; ++iy) {
			for (ix = 0; ix < Grid->nxVy; ++ix) {
				x = ix*Grid->dx+ Grid->xmin - 0.5*Grid->dx;
				y = iy*Grid->dy + Grid->ymin;
				if ( sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) < radius ) {

					if (assigning) {
					BC->list[I] = ix+iy*Grid->nxVy + Grid->nVxTot;
					BC->value[I] = 0.0;
					BC->type[I] = Dirichlet;
					}
					I++;


				}
			}
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





	if (!Grid->isPeriodic) {
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


		if (BC->SetupType==Stokes_WindTunnel) {
			// Extra BC for pressure
			C = Grid->nVxTot + Grid->nVyTot + Grid->nxEC-2 + Grid->nxEC;
			//C = Grid->nVxTot + Grid->nVyTot + 1 + Grid->nxEC;
			for (i=0;i<Grid->nyEC-2;i++){ // Prigth
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					BC->type[I] 		= Dirichlet;
				}
				I++;
				C = C+Grid->nxEC;
			}

		} else if (BC->SetupType==Stokes_Sandbox) {
			printf("koko\n");
			// Extra BC for pressure
			int ix, iy;
			for (iy=BC->iyTopRow;iy<Grid->nyEC-1;iy++) {
				for (ix=1;ix<Grid->nxEC-1;ix++) {
					if (assigning) {
						BC->list[I]         = Grid->nVxTot + Grid->nVyTot + ix + iy*Grid->nxEC;
						BC->value[I]        = 0.0;
						BC->type[I] 		= Dirichlet;
					}
					I++;
				}
			}
			printf("soko\n");
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
					BC->value[I]        = 0.0;//-1.0*Physics->rho_g[i]*Physics->gFac[1];
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
					//BC->value[I]        = 1.0*Physics->rho_f_g*Physics->gFac[1];//1.0*Physics->rho_g[i]*Physics->gFac[1];
					BC->value[I]        = 1.0*Physics->rho[i]*Physics->g[1];
					//BC->value[I]        = 1.0*Physics->rho_f_g*Physics->gFac[1] + 0.5*(1.0*Physics->rho_g[i]*Physics->gFac[1]-1.0*Physics->rho_f_g*Physics->gFac[1]);
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
					BC->type[I] 		= NeumannGhost;
				}
				I++;
				C = C+1;
			}
		} else if (iP == 0 ) {
			C = Grid->nVxTot + Grid->nVyTot + (Grid->nyEC-1)*Grid->nxEC + NumberMod;
			for (i=0;i<Grid->nxEC;i++){ // PTop
				if (assigning) {
					BC->list[I]         = C;
					BC->value[I]        = 0.0;
					/*
					if (Physics->y_oceanSurface < 0.0 + 1e-8) {
						// 0.0 is the default value
						BC->value[I]        = 0.0;
					} else {
						BC->value[I] =  Physics->rho_f_g*Physics->gFac[1]*(Grid->ymax-(Physics->y_oceanSurface+Grid->ymin));
					}
					*/
					BC->type[I] 		= Dirichlet;

					/*
					BC->value[I]        = 1.0*Physics->rho_g[i]*Physics->gFac[1];
					BC->type[I] 		= NeumannGhost;
					*/
				}
				I++;
				C = C+1;
			}
		}






		if (iP==1) { // Pc, i.e. Dummy, but actually important for interp
			if (!Grid->isPeriodic) {
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
		} else if (iP == 0 ) {
			if (!Grid->isPeriodic) {
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

						/*
						compute* Plitho_column;

						if (BC->SetupType == Stokes_Sandbox) {
							Plitho_column = (compute*) malloc(Grid->nyEC*sizeof(compute));
							int iy, iCell;
							int ix = Grid->nxEC-1;
							for (iy = 0; iy < Grid->nyEC; ++iy) {
								iCell = ix + iy*Grid->nxEC;
								compute y = Grid->ymin + iy*Grid->dy;
								Plitho_column[iy] = Physics->rho_g[iCell] * fabs(Physics->gFac[1]) * (Grid->ymax-y);
							}
						}
						*/

						C = Grid->nVxTot + Grid->nVyTot + Grid->nxEC-1 + Grid->nxEC + NumberMod;
						for (i=0;i<Grid->nyEC-2;i++){ // PRight
							if (assigning) {
								BC->list[I]         = C;
								BC->value[I]        = 0.0;
								BC->type[I] 		= NeumannGhost;

								/*
								if (BC->SetupType == Stokes_Sandbox) {

									compute y = Grid->ymin + i*Grid->dy;
									compute Plitho, Phydro;
									compute OvPFac = 0.25; // 0.0 is hydrostatic, 1.0 is lithostatic

									if (y<=BC->Sandbox_TopSeg01) {
										Phydro = Physics->rho_f_g * fabs(Physics->gFac[1]) * (Grid->ymax-y);
										Plitho = Plitho_column[i];
										BC->value[I]        = OvPFac*Plitho + (1.0-OvPFac)*Phydro;
										BC->type[I] 		= DirichletGhost;
									}
								}
								*/

							}
							I++;
							C = C+Grid->nxEC;
						}
						/*
						if (BC->SetupType == Stokes_Sandbox) {
							free(Plitho_column);
						}
						*/
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
					//BC->value[I]        = Physics->rho_g[i + (Grid->nyEC-2)*Grid->nxEC + 1] * fabs(Physics->gFac[1]) * Grid->dy/2  - Physics->Pf[i + (Grid->nyEC-2)*Grid->nxEC + 1];
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

	if (BC->SetupType==Thermal_TT_TB_LRNoFlux) {
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



		if (!Grid->isPeriodic) {


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

		}


	} else if (BC->SetupType==Thermal_TT_TBExternal_LRNoFlux) {
		// =======================================
		// =======================================
		// 				Pure Shear
		// =======================================
		// =======================================


		C = 0; // the first element in the numbering map is a ghost (in the sense of empty, i.e. there are no nodes in the corners)
		for (i=0; i<Grid->nxEC; i++) { // Bottom boundary
			if (assigning) {
				BC->list[I] = C;

				BC->value[I] = BC->TB; // Value at the external boundary
				BC->type [I] = Infinity;
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



		if (!Grid->isPeriodic) {


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

		}


	}

	else {
		printf("Unknown Temp BC.SetupType %i", BC->SetupType);
		exit(0);

	}


	BC->counter = I;
}



