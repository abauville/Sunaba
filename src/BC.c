/*
 * BC.c
 *
 *  Created on: Feb 23, 2016
 *      Author: abauville
 */

#include "stokes.h"
#define SQUARE(x) (x)*(x)

typedef enum { Vx, Vy, P } NodeType;
typedef enum { Row, Col } RowOrCol;
void assignBCToRowOrCol(NodeType NodeType, RowOrCol RowOrCol, int indRowOrCol, compute value, BCType BCType, int shift_start, int shift_end, Grid* Grid, BC* BC, bool assigning) {
	// if indRowOrCol is negative it indexes from the end backward. Like in numpy. e.g. -1 is end; -2 is end-1

	int I = BC->counter;
	int C;
	int increment;
	int increment_otherDirection;
	int i_start, i_end;
	int i, i_max;
	int ix_max, iy_max;
	if (NodeType==Vx) {
		C = 0;
		ix_max = Grid->nxVx;
		iy_max = Grid->nyVx;
	} else if (NodeType==Vy) {
		C = Grid->nVxTot;
		ix_max = Grid->nxVy;
		iy_max = Grid->nyVy;
	} else if (NodeType==P ) {
		C = Grid->nVxTot + Grid->nVyTot;
		ix_max = Grid->nxEC;
		iy_max = Grid->nyEC;
	} else {
		printf("error: unknown NodeType %i. Should be Vx, Vy or P", RowOrCol);
		exit(0);
	}

	if (RowOrCol == Col) {
		increment = ix_max;
		increment_otherDirection = 1;
		i_max = iy_max;
	} else if (RowOrCol == Row) {
		increment = 1;
		increment_otherDirection = ix_max;
		i_max = ix_max;
	} else {
		printf("error: unknown RowOrCol %i. Should be Row or Col", RowOrCol);
		exit(0);
	}

	
	C += shift_start*increment;

	if (indRowOrCol>=0) { 
		C += indRowOrCol * increment_otherDirection;
	} else {
		if (RowOrCol == Col) {
			C += (ix_max + indRowOrCol) * increment_otherDirection;
		} else if (RowOrCol == Row) {
			C += (iy_max + indRowOrCol) * increment_otherDirection;
		}
	}

	int C0 = C;


	for (i=shift_start; i<i_max-shift_end; i++) { // Vx Left
		if (assigning) {
			BC->list[I] = C;

			BC->value[I] = value;
			BC->type[I] = BCType;
			C += increment;
		}
		I+=1;
	}

	BC->counter = I;

	//printf("A: C0 = %i, Cend = %i, inc = %i, i_max = %i\n", C0, C, increment, i_max);
}
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


void BC_initStokes(Model* Model)
{
	BC* BC					= &(Model->BCStokes);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	EqSystem* EqSystem 		= &(Model->EqStokes);

	BC->reCompute_SymbolicFactorization = false;
	BC->iyTopRow = Grid->nyS;
	if (!DARCY) {
		BC->counter = 0;
		int nP, nV;
		BC_updateStokes_Vel(Model, false);
		nV = BC->counter;
		BC_updateStokes_P(Model, false);
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
		BC_updateStokes_Vel(Model, false);
		BC_updateStokesDarcy_P(Model, false);
		BC->n = BC->counter;
		EqSystem->nEq = EqSystem->nEqIni - BC->n;
		EqSystem->nRow = EqSystem->nEq;
	}





	BC->list    = (int*)     malloc( BC->n * sizeof(  int  ));
	BC->value   = (compute*) malloc( BC->n * sizeof(compute));
	BC->type   	= (BCType*)  malloc( BC->n * sizeof(BCType ));

	BC->counter = 0;

	if (!DARCY) {
		BC_updateStokes_Vel(Model, true);
		BC_updateStokes_P(Model, true);
	} else if (DARCY) {
		BC_updateStokes_Vel(Model, true);
		BC_updateStokesDarcy_P(Model, true);
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


void BC_initThermal(Model* Model)
{
	BC* BC					= &(Model->BCThermal);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	EqSystem* EqSystem 		= &(Model->EqThermal);

	BC->counter = 0;
	BC_updateThermal(Model, false);
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
	BC_updateThermal(Model, true);




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





void BC_updateStokes_Vel(Model* Model, bool assigning)
{
	BC* BC					= &(Model->BCStokes);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
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

		assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 0, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop
		
		I = BC->counter;
		
	} else if (BC->SetupType==Stokes_SimpleShear) {
		// =======================================
		// =======================================
		// Horizontal simple shear with lateral periodic BC
		// =======================================
		// =======================================1
		compute VxB =  2.0*1.0*BC->backStrainRate*Grid->ymin;
		compute VxT =  2.0*1.0*BC->backStrainRate*Grid->ymax;
		compute VyB =  0;
		compute VyT =  0;

		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vx, Row, 0, 0.0, DirichletGhost, 0, 0, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, DirichletGhost, 0, 0, Grid, BC, assigning); // VxTop
		
		I = BC->counter;

	} else if (BC->SetupType==Stokes_FixedLeftWall) {
		// =======================================
		// =======================================
		// 				Fixed Left Wall
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

		assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 0, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, DirichletGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning);  // VyRight
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning);  // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning);  // VxTop
		I = BC->counter;

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


		assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vx, Row, 0, VxL, DirichletGhost, 1, 1, Grid, BC, assigning); // VxBottom
		I = BC->counter;


		C = 1*Grid->nxVx-1;
		for (i=0; i<Grid->nyVx; i++) { // Vx Right
			if (assigning) {
				BC->list[I] = C;
				BC->value[I] = VxR;
				BC->type[I] = Dirichlet;

				// OutFlow
				compute y00 = (BC->Sandbox_TopSeg00 - (Grid->ymin + (i) * Grid->dy))/BC->Sandbox_TopSeg00;
				compute y01 = (BC->Sandbox_TopSeg01 - (Grid->ymin + (i) * Grid->dy))/(BC->Sandbox_TopSeg01-BC->Sandbox_TopSeg00);
				
				if (y00>0.0) {
					BC->value[I] = VxL;
				}
				else if (y01>0.0) {
					BC->value[I] = y01*VxL;
				}

				if (y01>0.0) {
					if (i>0) {
						integralOutflowVxdy += BC->value[I]*Grid->dy;
					}
				}
			}
			I++;
			C += Grid->nxVx;
		}
		extraOutFlowVy = -integralOutflowVxdy/((Grid->nxVy-2)*Grid->dx);


		int iy, ix;
		int highestTopo_iy = 0;
		int lowestTopo_iy = Grid->nyS+2;

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
		BC->iyTopRow_tolerance = (int) round(0.15*lowestTopo_iy);
		int iyTopRow_ideal = highestTopo_iy + nyAirMin;
		
		
		
		if (abs(iyTopRow_ideal-BC->iyTopRow)>BC->iyTopRow_tolerance &&  !(BC->iyTopRow == Grid->nyS && iyTopRow_ideal>Grid->nyS) ) {
		//if (abs(iyTopRow_ideal-BC->iyTopRow)>BC->iyTopRow_tolerance) {
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
		
		//printf("iyTopRow = %i, nTopRowsWithBC=%i, highestTopo_iy = %i, nyAirMin = %i\n", BC->iyTopRow, nTopRowsWithBC, highestTopo_iy, nyAirMin);
		
		for (iy=BC->iyTopRow-1;iy<Grid->nyVy;iy++) {
			for (ix=0; ix<Grid->nxVy; ix++) { // Vy Top
				if (assigning) {
					BC->list[I] = Grid->nVxTot + ix+iy*Grid->nxVy;//C;
					compute y = (Grid->Y[iy]-Grid->dy-Grid->ymin)/(Grid->ymax-Grid->ymin);

					BC->value[I] = (1.0-y)*VyB + y*VyT + extraOutFlowVy;
					BC->type[I] = Dirichlet;

				}
				I++;
			}
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



		C = Grid->nVxTot + 2*Grid->nxVy-1;
		for (i=0;i<Grid->nyVy-2-(nTopRowsWithBC-1);i++){ // Vy Right
			if (assigning) {
				BC->list[I]          = C;
				BC->value[I]         = 0.0;
				BC->type[I] 		 = NeumannGhost;

				
				if (BC->Sandbox_NoSlipWall) {
					if (Physics->phase[Grid->nxEC-1 + i * Grid->nxEC] != Physics->phaseAir && Physics->phase[Grid->nxEC-1 + i * Grid->nxEC] != Physics->phaseWater) {
						BC->type[I] 		 = DirichletGhost;
					}
				}

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
					

				}
				I++;
			}
		}
		
		

	} else if (BC->SetupType==Stokes_SandboxWeakBackstop) {
		// =======================================
		// =======================================
		// 			Sandbox weak backstop
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

		assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 1, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		//assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop
		
		I = BC->counter;

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
		BC->IsFreeSlipBot 	= true;
		BC->IsFreeSlipTop 	= true;

		//assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 0, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		//assignBCToRowOrCol(Vy, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop
		
		I = BC->counter;
		
		
		compute y0 = .8;
		compute beta = 30.0 * PI/180.0; // basal angle

		compute Vplate = VxL;
		compute VxPlate = Vplate * cos(beta);
		compute VyPlate = Vplate * sin(beta);
		compute x, y;
		

		assignBCToRowOrCol(Vx, Col, 0, VxPlate, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft

		I = BC->counter;

		compute inputVolume = 0.0;


		inputVolume += VxPlate * Grid->dy * (Grid->nyVx-1);

		// Neumann
		// =======================================
		C = Grid->nVxTot + Grid->nxVy;
		for (i=0;i<Grid->nyVy-2;i++){ // Vy Left
			if (assigning) {
				y = Grid->ymin+(i+1)*Grid->dy;
				BC->list[I]          = C;
				if (y<y0) {
					BC->value[I]         = VyPlate;
					BC->type[I] 		 = DirichletGhost;
				} else {
					BC->value[I]         = 0.0;
					BC->type[I] 		 = DirichletGhost;
				}

			}
			I++;
			C = C+Grid->nxVy;
		}
		

		/*

		// =======================================
		// =======================================
		// 			Sandbox internal BC
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

		

		assignBCToRowOrCol(Vx, Col, 0, VxL, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 0, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop
*/

/*
		assignBCToRowOrCol(Vx, Col,-1, VxR, Dirichlet, 0, 0, Grid, BC, assigning); // VxRight
		
		assignBCToRowOrCol(Vy, Row,-1, VyT, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop

		assignBCToRowOrCol(Vy, Row, 0, VyB, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom

		

		compute y0 = .8;
		compute beta = 30.0 * PI/180.0; // basal angle

		compute Vplate = VxL;
		compute VxPlate = Vplate * cos(beta);
		compute VyPlate = Vplate * sin(beta);

		assignBCToRowOrCol(Vx, Col, 0, VxPlate, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
*/
		//I = BC->counter;

		/*
		int ix, iy;
		compute x, y;
		// Vx loop
		for (ix=1;ix<Grid->nxVx-1;ix++) {
			x = Grid->xmin + ix*Grid->dx;
			for (iy=1;ix<Grid->nyVx-1;ix++) {
				y = Grid->ymin + iy*Grid->dy;
				if (y<=y0+x*tan(-beta)) {
					if (assigning) {
						BC->list[I]         = ix+iy*Grid->nxVx;
						BC->value[I]        = VxPlate;
						BC->type[I] 		= Dirichlet;
						I++;
					}
				}
			}
		}
		// Vy loop
		for (ix=1;ix<Grid->nxVy-1;ix++) {
			x = Grid->xmin + ix*Grid->dx;
			for (iy=1;ix<Grid->nyVy-1;ix++) {
				y = Grid->ymin + iy*Grid->dy;
				if (y<=y0+x*tan(-beta)) {
					if (assigning) {
						BC->list[I]         = Grid->nVxTot + ix+iy*Grid->nxVy;
						BC->value[I]        = VyPlate;
						BC->type[I] 		= Dirichlet;
						I++;
					}
				}
				
			}
		}


		C = 1;
		y = Grid->ymin;
		for (i=0;i<Grid->nxVx-1;i++){ // Vx Bottom
			x = Grid->xmin + i*Grid->dx;
			if (y<=y0+x*tan(-beta)) {
				if (assigning) {
					if (Physics->phase[i + 0*Grid->nxC] == BC->specialPhase) { // if special phase => dragging down
						BC->list[I]  = C;
						BC->value[I] = VxPlate;//(VxL+VxR)/2.0;
						BC->type[I]  = DirichletGhost;
					} else { 												   // else free slip
						BC->list[I]          = C;
						BC->value[I]         = 0.0;
						BC->type[I] = DirichletGhost;
					}
				}
				I++;
			C = C+1;
			}
		}

		


		C = Grid->nVxTot + 1;
		y = Grid->ymin;
		for (i=0;i<Grid->nxVy-1;i++){ // Vx Bottom
			x = Grid->xmin + i*Grid->dx;
			if (y<=y0+x*tan(-beta)) {
				if (assigning) {
					if (Physics->phase[i + 0*Grid->nxC] == BC->specialPhase) { // if special phase => dragging down
						BC->list[I]  = C;
						BC->value[I] = VyPlate;//(VxL+VxR)/2.0;
						BC->type[I]  = DirichletGhost;
					} else { 												   // else free slip
						BC->list[I]  = C;
						BC->value[I] = 0.0;
						BC->type[I]  = DirichletGhost;
					}
				}
				I++;
			C = C+1;
			}
		}

		*/

		


		





	} else if (BC->SetupType==Stokes_CornerFlow) {
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
				y = (Grid->ymin + Grid->dy*i);
				iy = i+1;

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
				y = (Grid->ymin + Grid->dy*i);
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


		compute VxInput = BC->refValue;

		BC->IsFreeSlipLeft	= true;
		BC->IsFreeSlipRight = true;
		BC->IsFreeSlipBot 	= true;
		BC->IsFreeSlipTop 	= true;

	

	
		assignBCToRowOrCol(Vx, Col, 0, VxInput, Dirichlet, 0, 0, Grid, BC, assigning); // VxLeft
		assignBCToRowOrCol(Vx, Col,-1, 0.0, Neumann, 0, 0, Grid, BC, assigning); // VxRight
		assignBCToRowOrCol(Vy, Row, 0, 0.0, Dirichlet, 0, 0, Grid, BC, assigning); // VyBottom
		assignBCToRowOrCol(Vy, Row,-1, 0.0, Dirichlet, 0, 0, Grid, BC, assigning); // VyTop

		assignBCToRowOrCol(Vy, Col, 0, 0.0, Dirichlet, 1, 1, Grid, BC, assigning); // Vyleft
		assignBCToRowOrCol(Vy, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VyRight
		assignBCToRowOrCol(Vx, Row, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxBottom
		assignBCToRowOrCol(Vx, Row,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // VxTop
		
		I = BC->counter;
		


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
				if ( sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) < radius ) {

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



	} else {
		printf("Unknown Stokes BC.SetupType %i", BC->SetupType);
		exit(0);
	}

	BC->counter = I;
}



void BC_updateStokes_P(Model* Model, bool assigning)
{
	BC* BC					= &(Model->BCStokes);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);	
	int C, I, i;

	I = BC->counter;
	// Pressure BC for all setup
	// =======================================
	// in normal stokes there is lagrangian operator on the Pressure, only the gradient
	// and it doesn't have to be build for boundary velocity nodes (since they are BC)
	// therefore these are just dummy values in this case
	// However the following acts as Pf Boundary conditions for the two-phase flow

	assignBCToRowOrCol(P, Row, 0, 0.0, NeumannGhost, 0, 0, Grid, BC, assigning); // PBottom
	assignBCToRowOrCol(P, Row,-1, 0.0, NeumannGhost, 0, 0, Grid, BC, assigning); // PTop
	I = BC->counter;
	
	if (!Grid->isPeriodic) {
		assignBCToRowOrCol(P, Col, 0, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // PLeft
		assignBCToRowOrCol(P, Col,-1, 0.0, NeumannGhost, 1, 1, Grid, BC, assigning); // PRight
		I = BC->counter;
	}

	if (BC->SetupType==Stokes_WindTunnel) {
		assignBCToRowOrCol(P, Col,-2, 0.0, Dirichlet, 1, 1, Grid, BC, assigning); // PRight
		I = BC->counter;
		
	} else if (BC->SetupType==Stokes_Sandbox) {
		// Extra BC for pressure
		int ix, iy;
		for (iy=BC->iyTopRow-1;iy<Grid->nyEC-1;iy++) {
			for (ix=1;ix<Grid->nxEC-1;ix++) {
				if (assigning) {
					BC->list[I]         = Grid->nVxTot + Grid->nVyTot + ix + iy*Grid->nxEC;
					BC->value[I]        = 0.0;
					BC->type[I] 		= Dirichlet;
				}
				I++;
			}
		}
	}
	BC->counter = I;


}


#if (DARCY)
void BC_updateStokesDarcy_P(Model* Model, bool assigning) {
	BC* BC					= &(Model->BCStokes);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
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

void BC_updateThermal(Model* Model, bool assigning)
{
	BC* BC					= &(Model->BCThermal);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
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


