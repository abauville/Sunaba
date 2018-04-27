/*
 * Interp.c
 *
 *  Created on: Jul 26, 2017
 *      Author: abauville
 */


#include "stokes.h"

#define TEST_SIGMA_INTERP false
#define TEST_SIGMA_INTERP_FROM_PART_TO_CELL true // if false, eulerian only
#define PART2GRID_SCHEME 0  // 0 local scheme (Taras), each Particle contributes to only one node or cell (domain area: dx*dy)
						   	// 1 wide scheme (Mikito), each Particle contributes to only 4 nodes or cells (domain area: 2*dx * 2*dy)
#define USE_CLOSEST_GRID2PART false // false is linear interpolation, true is closest neighbour
#define USE_SPECIAL_STRESS_INTERP false

inline compute Interp_ECVal_Cell2Particle_Local(compute* A, int ix, int iy, int nxEC, compute locX, compute locY)
{
	// Compute a value on particles from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node the particle is attached to

/*
     	locX=-1     locX=+1
			2 x ------- x 3   locY=+1
			  |         |
			  |    O    |
			  |         |
			1 X ------- x 4   locY=-1

			O: node to which the particles are attached (has index ix, iy)
			x: Cells
			X: Cell with ix,iy index
*/
#if (USE_CLOSEST_GRID2PART)
	if (locX>0.0) {
		ix = ix+1;
	}
	if (locY>0.0) {
		iy = iy+1;
	}
	return A[ix  +(iy  )*nxEC];
#else
	return ( .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxEC]
           + .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxEC]
		   + .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxEC]
		   + .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxEC] );
#endif

}



inline compute Interp_NodeVal_Node2Particle_Local(compute* A, int ix, int iy, int nxS, int nyS, compute locX, compute locY) {
	// Compute a value on particles from a Array of values defined on the shear node grid
	// where ix and iy refer to shear node the particle is attached to
	/*
     	locX=-1     locX=+1
			2 o ------- o 3   locY=+1
			  |         |
			  |         |
			  |         |
			1 O ------- o 4   locY=-1

			O: node to which the particles are attached (has index ix, iy)
			o: Nodes

*/
#if (USE_CLOSEST_GRID2PART)
	return A[ix      +(iy  )    *nxS];
#else
	int signX, signY;
	if (locX<0.0) {
		signX = -1;
	} else {
		signX = 1;
	}
	if (locY<0.0) {
		signY = -1;
	} else {
		signY = 1;
	}
	locX = fabs(locX)-1.0;
	locY = fabs(locY)-1.0;
	return ( .25*(1.0-locX)*(1.0-locY)*A[ix      +(iy  )    *nxS]
		   + .25*(1.0-locX)*(1.0+locY)*A[ix      +(iy+signY)*nxS]
		   + .25*(1.0+locX)*(1.0+locY)*A[ix+signX+(iy+signY)*nxS]
		   + .25*(1.0+locX)*(1.0-locY)*A[ix+signX+(iy  )    *nxS] );
#endif

}

compute Interp_Special_Sxx_Cell2Particle_Local(compute* Sxx, compute* Epxx, int ix, int iy, int nxEC, compute locX, compute locY)
{
	// Compute a value on particles from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node the particle is attached to

	compute* A = Sxx;



	compute particleValue; //value on the given particle;

	compute locXN = locX;
	compute locYN = locY;
	int ixN = ix;
	int iyN = iy;

	
	

	compute defVal = 1.0;

	int ixCSW = ixN;
	int iyCSW = iyN;
	
	int ixCSE = ixN+1;
	int iyCSE = iyN;
	
	int ixCNW = ixN;
	int iyCNW = iyN+1;
	
	int ixCNE = ixN+1;
	int iyCNE = iyN+1;
	
	int Counter = 0;
	if (A[ixCSW+iyCSW*nxEC] > defVal)
		Counter += 1;
	if (A[ixCSE+iyCSE*nxEC] > defVal)
		Counter += 1;
	if (A[ixCNW+iyCNW*nxEC] > defVal)
		Counter += 1;
	if (A[ixCNE+iyCNE*nxEC] > defVal)
		Counter += 1;
		

	int ixC, iyC, signX, signY;


    
	if ( (Counter == 1 && Epxx[ixCSE+iyCSE*nxEC]==defVal && Epxx[ixCNE+iyCNE*nxEC]<defVal) || (Counter == 1 && Epxx[ixCNW+iyCNW*nxEC]==defVal && Epxx[ixCSW+iyCSW*nxEC]<defVal) || (Counter == 2 && Epxx[ixCNW+iyCNW*nxEC]<defVal && Epxx[ixCSE+iyCSE*nxEC]<defVal) ) {// # NW-SE diagonal
	//if (Counter == 1 and A[ixCSE,iyCSE]==defVal and A[ixCNE,iyCNE]<defVal) or (Counter == 1 and A[ixCNW,iyCNW]==defVal and A[ixCSW,iyCSW]<defVal) or (A[ixCNW,iyCNW]<defVal and A[ixCSE,iyCSE]<defVal): # NW-SE diagonal
//                if (Counter == 2 and A[ixCNW,iyCNW]<defVal and A[ixCSE,iyCSE]<defVal): # NW-SE diagonal
	//if (fabs(Epxx[ixCNW+iyCNW*nxEC])<defVal && fabs(Epxx[ixCSE+iyCSE*nxEC])<defVal) { // NW-SE diagonal
		if (locYN<-locXN) {
			ixC = ixN;
			iyC = iyN;
			signX = 1;
			signY = 1;
			locX = (locXN+1.0)/2.0;
			locY = (locYN+1.0)/2.0;
		} else {
			ixC = ixN+1;
			iyC = iyN+1;
			signX = -1;
			signY = -1;
			locX = -(locXN-1.0)/2.0;
            locY = -(locYN-1.0)/2.0;
		}

		particleValue = ( (1.0-locX-locY)*A[ixC      +iyC*nxEC      ]
						   + locX        *A[ixC+signX   +iyC*nxEC      ]
						   + locY        *A[ixC         +(iyC+signY)*nxEC] );
	} else if ( (Counter == 1 && Epxx[ixCNE+iyCNE*nxEC]==defVal && Epxx[ixCSE+iyCSE*nxEC]<defVal) || (Counter == 1 && Epxx[ixCSW+iyCSW*nxEC]==defVal && A[ixCNW+iyCNW*nxEC]<defVal) || (Counter == 2 && Epxx[ixCSW+iyCSW*nxEC]<defVal && Epxx[ixCNE+iyCNE*nxEC]<defVal) ) { //# SW-NE diagonal
//                elif (Counter == 1 and A[ixCNE,iyCNE]==defVal and A[ixCSE,iyCSE]<defVal) or (Counter == 1 and A[ixCSW,iyCSW]==defVal and A[ixCNW,iyCNW]<defVal) or (A[ixCSW,iyCSW]<defVal and A[ixCNE,iyCNE]<defVal): # SW-NE diagonal
//                elif (Counter == 2 and A[ixCSW,iyCSW]<defVal and A[ixCNE,iyCNE]<defVal): # SW-NE diagonal
	//} else if (fabs(Epxx[ixCSW+iyCSW*nxEC])<defVal && fabs(Epxx[ixCNE+iyCNE*nxEC])<defVal) { // SW-NE diagonal
		if (locYN<locXN) {
			ixC = ixN+1;
			iyC = iyN;
			signX = -1;
			signY =  1;
			locX = -(locXN-1.0)/2.0;
			locY =  (locYN+1.0)/2.0;
		} else {
			ixC = ixN;
			iyC = iyN+1;
			signX =  1;
			signY = -1;
			locX =  (locXN+1.0)/2.0;
			locY = -(locYN-1.0)/2.0;
		}

		particleValue = ( (1.0-locX-locY)*A[ixC      +iyC*nxEC      ]
						+ locX           *A[ixC+signX+iyC*nxEC      ]
						+ locY           *A[ixC      +(iyC+signY)*nxEC] );     
	
	} else {
		ixC = ixN;
		iyC = iyN;
		locX = locXN;
		locY = locYN;
		particleValue =  ( .25*(1.0-locX)*(1.0-locY)*A[ixC  + iyC   *nxEC]
					 	 + .25*(1.0-locX)*(1.0+locY)*A[ixC  +(iyC+1)*nxEC]
						 + .25*(1.0+locX)*(1.0+locY)*A[ixC+1+(iyC+1)*nxEC]
						 + .25*(1.0+locX)*(1.0-locY)*A[ixC+1+iyC    *nxEC]);
	
	}
	
	return particleValue;
	
	
	

}

compute Interp_Special_Sxy_Node2Particle_Local(compute* Sxy, compute* Epxy, int ix, int iy, int nxS, int nyS, compute locX, compute locY) {
	// Compute a value on particles from a Array of values defined on the shear node grid
	// where ix and iy refer to shear node the particle is attached to

	compute* A = Sxy;

	compute particleValue;
	compute defVal = 1.0;
	compute locXN = locX;
	compute locYN = locY;

	int ixN = ix;
	int iyN = iy;

	int signX, signY;
	int ixSW, ixSE, ixNE, ixNW;
	int iySW, iySE, iyNE, iyNW;

	if (locXN>0.0) {
		signX =  1;
		ixSW = ixN;
		ixSE = ixN+1;
		ixNE = ixN+1;
		ixNW = ixN;
	} else {
		signX = -1;
		ixSW = ixN-1;
		ixSE = ixN;
		ixNE = ixN;
		ixNW = ixN-1;
	}
	
	if (locYN>0.0) {
		signY =  1;
		iySW = iyN;
		iySE = iyN;
		iyNE = iyN+1;
		iyNW = iyN+1;
	} else {
		signY = -1;
		iySW = iyN-1;
		iySE = iyN-1;
		iyNE = iyN;
		iyNW = iyN;
	}
	

	if (fabs(Epxy[ixSW+iySW*nxS])<defVal && fabs(Epxy[ixNE+iyNE*nxS])<defVal) { // SW-NE diagonal
		locX = fabs((fabs(locX)-1.0)/2.0);
		locY = fabs((fabs(locY)-1.0)/2.0);  
		
		if (locXN>0.0) {
 			locX = locXN/2.0;
		} else {
			locX = locXN/2.0+1.0;
		}
		if (locYN>0.0) {
			locY = locYN/2.0;
		} else {
			locY = locYN/2.0+1.0;
		}
			
		if (locY<locX) {
			locX = fabs(locX-1.0);
			particleValue = ( (1.0-locX-locY)*A[ixSE    +iySE*nxS   ]
							+ locX           *A[ixSW    +iySW*nxS   ]
							+ locY           *A[ixNE    +iyNE*nxS   ] );   
		} else {
			locY = fabs(locY-1.0);
			particleValue = ( (1.0-locX-locY)*A[ixNW     +iyNW*nxS   ]
							+ locX           *A[ixNE     +iyNE*nxS   ]
							+ locY           *A[ixSW     +iySW*nxS   ] );
		}


	} else if (fabs(Epxy[ixNW+iyNW*nxS])<defVal && fabs(Epxy[ixSE+iySE*nxS])<defVal) { // NW-SE diagonal
		if (locXN>0.0) {
			locX = locXN/2.0;
		} else {
			locX = locXN/2.0+1.0;
		}
		if (locYN>0.0) {
			locY = locYN/2.0;
		} else {
			locY = locYN/2.0+1.0;
		}
			
		if (locY<=-locX+1.0) {
			particleValue = ( (1.0-locX-locY)*A[ixSW    +iySW*nxS   ]
							+ locX           *A[ixSE    +iySE*nxS   ]
							+ locY           *A[ixNW    +iyNW*nxS   ] );
		} else {
			locX = fabs(locX-1.0);
			locY = fabs(locY-1.0);
			particleValue = ( (1.0-locX-locY)*A[ixNE    +iyNE*nxS   ]
							+ locX           *A[ixNW    +iyNW*nxS   ]
							+ locY           *A[ixSE    +iySE*nxS   ] );   
		}	
	} else {
		locX = fabs(locXN)-1.0;
		locY = fabs(locYN)-1.0;

		particleValue =  ( .25*(1.0-locX)*(1.0-locY)*A[ixN      + iyN       *nxS]
						 + .25*(1.0-locX)*(1.0+locY)*A[ixN      +(iyN+signY)*nxS]
						 + .25*(1.0+locX)*(1.0+locY)*A[ixN+signX+(iyN+signY)*nxS]
						 + .25*(1.0+locX)*(1.0-locY)*A[ixN+signX+iyN        *nxS]); 
	}

	return particleValue;
	

	

}


inline compute Interp_Product_ECVal_Cell2Particle_Local(compute* A, compute* B, int ix, int iy, int nxEC, compute locX, compute locY)
{
	// Compute the product of the NodeVal, A and B, and interpolates it to the particles

	/*
			locX=-1     locX=+1
				2 x ------- x 3   locY=+1
				|         |
				|    O    |
				|         |
				1 X ------- x 4   locY=-1

				O: node to which the particles are attached (has index ix, iy)
				x: Cells
				X: Cell with ix,iy index
	*/

	return (  .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxEC]*B[ix  +(iy  )*nxEC]
			+ .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxEC]*B[ix  +(iy+1)*nxEC]
			+ .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxEC]*B[ix+1+(iy+1)*nxEC]
			+ .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxEC]*B[ix+1+(iy  )*nxEC] );
}

inline compute Interp_Product_NodeVal_Node2Particle_Local(compute* A, compute* B, int ix, int iy, int nxS, int nyS, compute locX, compute locY) 
{
	// Compute the product of the NodeVal, A and B, and interpolates it to the particles
	// where ix and iy refer to shear node the particle is attached to
	/*
			locX=-1     locX=+1
			2 o ------- o 3   locY=+1
				|         |
				|         |
				|         |
			1 O ------- o 4   locY=-1

			O: node to which the particles are attached (has index ix, iy)
			o: Nodes

	*/
	int signX, signY;
	if (locX<0.0) {
		signX = -1;
	} else {
		signX = 1;
	}
	if (locY<0.0) {
		signY = -1;
	} else {
		signY = 1;
	}
	locX = fabs(locX)-1.0;
	locY = fabs(locY)-1.0;
	return (  .25*(1.0-locX)*(1.0-locY)*A[ix      +(iy  )    *nxS]*B[ix      +(iy  )    *nxS]
			+ .25*(1.0-locX)*(1.0+locY)*A[ix      +(iy+signY)*nxS]*B[ix      +(iy+signY)*nxS]
			+ .25*(1.0+locX)*(1.0+locY)*A[ix+signX+(iy+signY)*nxS]*B[ix+signX+(iy+signY)*nxS]
			+ .25*(1.0+locX)*(1.0-locY)*A[ix+signX+(iy  )    *nxS]*B[ix+signX+(iy  )    *nxS] );

}

inline compute Interp_Product_ECNodeVal_Cell2Particle_Local(compute* A, compute* B, int ix, int iy, int nxEC, int nxS, compute locX, compute locY)
{
	// Compute the product of the ECVal, A, and NodeVal, B, and interpolates it to the particles
	/*
			  locX=-1     locX=+1
				  | 		|
		|    N    |	   N    |    N    |
				  |			|
				2 x ------- x 3   locY=+1
				  |         |
		|    N    |    O    |    N    |
				  |         |
				1 X ------- x 4   locY=-1
				  |			|
		|    N    |	   N	|    N    |
				  | 		|
				O: node to which the particles are attached (has index ix, iy)
				N: surrounding nodes
				x: Cells
				X: Cell with ix,iy index

		First Node values are interpolated at the cell centers at which the product is computed.
		Then, the product is interpolated from the cell centers to the particle
	*/

	compute B_at_Cells[4];
	
	B_at_Cells[0] = Interp_NodeVal_Node2Cell_Local(B, ix  , iy  , nxS);
	B_at_Cells[1] = Interp_NodeVal_Node2Cell_Local(B, ix  , iy+1, nxS);
	B_at_Cells[2] = Interp_NodeVal_Node2Cell_Local(B, ix+1, iy+1, nxS);
	B_at_Cells[3] = Interp_NodeVal_Node2Cell_Local(B, ix+1, iy  , nxS);
	
	return (  .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxEC]*B_at_Cells[0]
			+ .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxEC]*B_at_Cells[1]
			+ .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxEC]*B_at_Cells[2]
			+ .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxEC]*B_at_Cells[3] );
}





inline compute Interp_VxVal_VxNode2Particle_Local(compute* A, int ix, int iy, int nxVx, compute locX, compute locY)
{
	// Compute a value on particles from a Array of values defined on the Vx nodes
	// where ix and iy refer to shear node the particle is attached to

	if (locX>0.0) {
		locX = locX-1.0;
	}
	else {
		locX = locX+1.0;
		ix-=1;
	}

	return  ( .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxVx]
			+ .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxVx]
			+ .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxVx]
			+ .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxVx] );

}


inline compute Interp_VyVal_VyNode2Particle_Local(compute* A, int ix, int iy, int nxVy, compute locX, compute locY)
{
	// Compute a value on particles from a Array of values defined on the Vy nodes
	// where ix and iy refer to shear node the particle is attached to

	if (locY>0.0) {
		locY = locY-1.0;
	}
	else {
		locY = locY+1.0;
		iy -= 1;
	}

	return  ( .25*(1.0-locX)*(1.0-locY)*A[ix  +(iy  )*nxVy]
			+ .25*(1.0-locX)*(1.0+locY)*A[ix  +(iy+1)*nxVy]
			+ .25*(1.0+locX)*(1.0+locY)*A[ix+1+(iy+1)*nxVy]
			+ .25*(1.0+locX)*(1.0-locY)*A[ix+1+(iy  )*nxVy] );


}










inline compute Interp_ECVal_Cell2Node_Local(compute* A, int ix, int iy, int nxEC)
{
	// Compute a value on the shear grid from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy+1)*nxEC] + A[ix+1+(iy+1)*nxEC] + A[ix  +(iy  )*nxEC] + A[ix+1+(iy  )*nxEC])/4.0;
}

inline compute Interp_NodeVal_Node2Cell_Local(compute* A, int ix, int iy, int nxS)
{
	// Compute a value on an embedded cell center from the A Array of values defined on the shear grid
	// where ix and iy refer to cell centers
	return(A[ix  +(iy-1)*nxS] + A[ix-1+(iy-1)*nxS] + A[ix  +(iy  )*nxS] + A[ix-1+(iy  )*nxS])/4.0;
}

inline compute Interp_Product_ECVal_Cell2Node_Local(compute* A, compute* B, int ix, int iy, int nxEC)
{
	// Compute a value on the shear grid from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy+1)*nxEC]*B[ix  +(iy+1)*nxEC]
		 + A[ix+1+(iy+1)*nxEC]*B[ix+1+(iy+1)*nxEC] 
		 + A[ix  +(iy  )*nxEC]*B[ix  +(iy  )*nxEC]
		 + A[ix+1+(iy  )*nxEC]*B[ix+1+(iy  )*nxEC])/4.0;
}

inline compute Interp_Product_NodeVal_Node2Cell_Local(compute* A, compute* B, int ix, int iy, int nxS)
{
	// Compute a value on an embedded cell center from the A Array of values defined on the shear grid
	// where ix and iy refer to cell centers
	return(A[ix  +(iy-1)*nxS]*B[ix  +(iy-1)*nxS] 
		 + A[ix-1+(iy-1)*nxS]*B[ix-1+(iy-1)*nxS] 
		 + A[ix  +(iy  )*nxS]*B[ix  +(iy  )*nxS] 
		 + A[ix-1+(iy  )*nxS]*B[ix-1+(iy  )*nxS])/4.0;
}



void Interp_All_Particles2Grid_Global(Model* Model)
{


	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);





	// Declarations
	// =========================
	int iCell, iNode;
	//int nNeighbours = 4;
	coord locX, locY;


	// Reinitialize Physics array
	bool* changedHead = (bool*) malloc(Grid->nECTot * sizeof(bool));


#pragma omp parallel for private(iCell) OMP_SCHEDULE
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
		Physics->sigma_xx_0 [iCell] = 0.0;
#endif
		Physics->sumOfWeightsCells [iCell] = 0.0;
		Physics->G[iCell] = 0.0;
#if (HEAT)

		Physics->T[iCell] = 0.0;
#endif
#if (DARCY)
		Physics->DeltaP0[iCell] 		= 0.0;
		Physics->phi0[iCell] 		= 0.0;
#endif
#if (STORE_PLASTIC_STRAIN)
		Physics->strain[iCell] 		= 0.0;
#endif

		changedHead[iCell] = false;
	}


#pragma omp parallel for private(iNode) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
		Physics->sigma_xy_0 [iNode] = 0.0;
		Physics->GShear [iNode] = 0.0;
#endif
		Physics->sumOfWeightsNodes [iNode] = 0.0;
	}


	Physics_PhaseList_reinit(Model);






	compute weight;
	int phase;
	int nxEC = Grid->nxEC;
	compute xMod[4], yMod[4];
	int ix,  iy;



	xMod[0] = -1.0; yMod[0] = -1.0;
	xMod[1] =  1.0; yMod[1] = -1.0;
	xMod[2] = -1.0; yMod[2] =  1.0;
	xMod[3] =  1.0; yMod[3] =  1.0;



	// Index of neighbouring cells, with respect to the node ix, iy
	int i;
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right


	SingleParticle* thisParticle = NULL;


	int iColor; // indexing of the color group for nodes. Nodes of the same color don't collide with each other. i.e. similar to matrix coloring
	int ixStart[4] = {0,0,1,1};
	int iyStart[4] = {0,1,0,1};
	SinglePhase* thisPhaseInfo;

	for (iColor = 0; iColor < 4; ++iColor) {
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iCell, weight, thisPhaseInfo) OMP_SCHEDULE
		for (iy = iyStart[iColor]; iy < Grid->nyS; iy+=2) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ix+=2) { // I don't get why though
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				// Loop through the particles in the shifted cell
				// ======================================
				while (thisParticle!=NULL) {
					locX = Particles_getLocX(ix, thisParticle->x,Grid);
					locY = Particles_getLocY(iy, thisParticle->y,Grid);
					/*
					if (fabs(locX)>1.0+1e-2 || fabs(locY)>1.0+1e-2 ) {
						printf("Error locXY, locX = %.4f, locY = %.4f\n", locX, locY);
						exit(0);p
					}
					*/
					


					phase = thisParticle->phase;
#if (PART2GRID_SCHEME == 1)
					for (i=0; i<4; i++) {
						iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
						weight = fabs((locX + xMod[i])   *   (locY + yMod[i]));
#elif (PART2GRID_SCHEME == 0)
						int signX, signY;
						if (locX<0.0) {
							signX = -1;
						} else {
							signX = 1;
						}
						if (locY<0.0) {
							signY = -1;
						} else {
							signY = 1;
						}
						if 		 	(signX>=0 && signY>=0) { // upper right
							i = 3;
						} else if 	(signX<0 && signY>=0) { // upper left
							// the particle is in the SE quadrant, the cell center 1 is NW (wrt to the node ix,iy)
							i = 2;
						} else if 	(signX>=0 && signY<0) { // lower right
							i = 1;
						} else if 	(signX<0 && signY<0) { // lower left
							i = 0;
						} else {
							printf("error in Interp_ECVal_Cell2Particle_Local. No case was triggered\n.");
							exit(0);
						}
						iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
						weight = fabs(locX)*fabs(locY);

#endif
						


						// Get the phase and weight of phase contribution for each cell
						thisPhaseInfo = Physics->phaseListHead[iCell];

						while (thisPhaseInfo->phase != phase) {
							if (thisPhaseInfo->next == NULL) {

								if (!changedHead[iCell]) {
									thisPhaseInfo->phase = phase;
									changedHead[iCell] = true;
									break;
								} else {

									Physics_Phase_addSingle(&Physics->phaseListHead[iCell],phase);
									thisPhaseInfo = Physics->phaseListHead[iCell];
									break;
								}


							} else {
								thisPhaseInfo = thisPhaseInfo->next;
							}
						}
						thisPhaseInfo->weight += weight;


						// For properties that are stored on the markers, sum contributions
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
						Physics->sigma_xx_0		[iCell] += thisParticle->sigma_xx_0 * weight;
						Physics->G				[iCell] += weight / MatProps->G[phase];
#endif
#if (HEAT)
						Physics->T				[iCell] += thisParticle->T * weight;
#endif
#if (DARCY)
						Physics->DeltaP0		[iCell] += thisParticle->DeltaP0 * weight;
						Physics->phi0			[iCell] += thisParticle->phi * weight;
#endif
#if (STORE_PLASTIC_STRAIN)
						Physics->strain			[iCell] += thisParticle->strain * weight;
#endif
						Physics->sumOfWeightsCells	[iCell] += weight;
#if (PART2GRID_SCHEME == 1)
					}
#endif
					thisParticle = thisParticle->next;
				}
			}
		}
	}

#if (PART2GRID_SCHEME == 0)
	// For this scheme, outer cells receive no contribution from particles
	// Not so important because calculation is not made on them
	// But to avoid division by 0, I here copy the values from the neighbours anyway.
	// Also this allows to check for empty cells.
	
	for (iy=1;iy<Grid->nyEC-1;iy++) {
		for (ix=1;ix<Grid->nxEC-1;ix++) {
			iCell = ix + iy*Grid->nxEC;
			if (Physics->sumOfWeightsCells	[iCell] == 0) {
				// If no contributions was given to this cell (i.e. empty cell), then use a higher order interpolation scheme (2x2 cells wide instead of 1x1)
				int iNodeCounter = 0;
				for (iNodeCounter=0;iNodeCounter<4;iNodeCounter++) {

					iNode = ix  + (iy  )*Grid->nxS;
					thisParticle = Particles->linkHead[iNode];
					// Loop through the particles in the shifted cell
					// ======================================
					while (thisParticle!=NULL) {
						locX = Particles_getLocX(ix, thisParticle->x,Grid);
						locY = Particles_getLocY(iy, thisParticle->y,Grid);


						phase = thisParticle->phase;
						for (i=0; i<4; i++) {
							int thisCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
							if (thisCell==iCell) {
								weight = fabs((locX + xMod[i])   *   (locY + yMod[i]));
								// Get the phase and weight of phase contribution for each cell
								thisPhaseInfo = Physics->phaseListHead[iCell];

								while (thisPhaseInfo->phase != phase) {
									if (thisPhaseInfo->next == NULL) {

										if (!changedHead[iCell]) {
											thisPhaseInfo->phase = phase;
											changedHead[iCell] = true;
											break;
										} else {

											Physics_Phase_addSingle(&Physics->phaseListHead[iCell],phase);
											thisPhaseInfo = Physics->phaseListHead[iCell];
											break;
										}


									} else {
										thisPhaseInfo = thisPhaseInfo->next;
									}
								}
								thisPhaseInfo->weight += weight;


								// For properties that are stored on the markers, sum contributions
		#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
								Physics->sigma_xx_0		[iCell] += thisParticle->sigma_xx_0 * weight;
								Physics->G				[iCell] += weight / MatProps->G[phase];
		#endif
		#if (HEAT)
								Physics->T				[iCell] += thisParticle->T * weight;
		#endif
		#if (DARCY)
								Physics->DeltaP0		[iCell] += thisParticle->DeltaP0 * weight;
								Physics->phi0			[iCell] += thisParticle->phi * weight;
		#endif
		#if (STORE_PLASTIC_STRAIN)
								Physics->strain			[iCell] += thisParticle->strain * weight;
		#endif
								Physics->sumOfWeightsCells	[iCell] += weight;
							} // if thisCell=iCell
						} // for neighbour cells
						thisParticle = thisParticle->next;
					} // while Particles
				} // iNodeCounter
			} //if (Physics->sumOfWeightsCells	[iCell] == 0) 
		} // ixCell
	} // iyCell
	
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sumOfWeightsCells, Grid);

	


#endif


	// Copy contribution from one side to the other in case of periodic BC
	if(Grid->isPeriodic) {
		int iCellS, iCellD, j;
#pragma omp parallel for private(iy, j, iCellS, iCellD) OMP_SCHEDULE
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			for (j = 0; j<2; ++j) {
				iCellS = j + iy*Grid->nxEC;
				iCellD = Grid->nxEC-2+j + iy*Grid->nxEC;

				Physics->sigma_xx_0		[iCellD] += Physics->sigma_xx_0		[iCellS];
				Physics->sigma_xx_0		[iCellS]  = Physics->sigma_xx_0		[iCellD];

				Physics->G		[iCellD] += Physics->G		[iCellS];
				Physics->G		[iCellS]  = Physics->G		[iCellD];
#if (HEAT)
				Physics->T				[iCellD] += Physics->T				[iCellS];
				Physics->T				[iCellS]  = Physics->T				[iCellD];
#endif
#if (DARCY)
				Physics->DeltaP0		[iCellD] += Physics->DeltaP0		[iCellS];
				Physics->DeltaP0		[iCellS]  = Physics->DeltaP0		[iCellD];
				Physics->phi0			[iCellD] += Physics->phi0			[iCellS];
				Physics->phi0			[iCellS]  = Physics->phi0			[iCellD];
#endif
#if (STORE_PLASTIC_STRAIN)
				Physics->strain			[iCellD] += Physics->strain			[iCellS];
				Physics->strain			[iCellS]  = Physics->strain			[iCellD];
#endif
				Physics->sumOfWeightsCells	[iCellD] += Physics->sumOfWeightsCells	[iCellS];
				Physics->sumOfWeightsCells	[iCellS]  = Physics->sumOfWeightsCells	[iCellD];

			}
		}
	}

	// Filling side values
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sumOfWeightsCells, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sigma_xx_0, Grid);
#if (HEAT)
	Physics_CellVal_SideValues_getFromBC(Physics->T, Grid, BCThermal, NumThermal);
#endif
#if (DARCY)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->DeltaP0, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->phi0, Grid);
#endif
#if (STORE_PLASTIC_STRAIN)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->strain, Grid);
#endif

	free(changedHead);
	// Dividing by the sum of weights
//#pragma omp parallel for private(iCell) OMP_SCHEDULE
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		//printf("sumOfWeights[%i] = %.2e\n", iCell, Physics->sumOfWeightsCells	[iCell]);
		
		if (Physics->sumOfWeightsCells	[iCell]==0.0) {
			printf("error in Interp_All_Particles2Grid_Global. Cell #%i received no contribution from particles (i.e. empty cell).\n", iCell);
			exit(0);
		}
		

#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
		Physics->sigma_xx_0	[iCell] /= Physics->sumOfWeightsCells	[iCell];


#endif
		Physics->G			[iCell]  = Physics->sumOfWeightsCells	[iCell] / Physics->G	[iCell];
		
		//Physics->sigma_xx_0	[iCell] = Physics->sumOfWeightsCells	[iCell] / Physics->sigma_xx_0	[iCell];
#if (HEAT)
		Physics->T			[iCell] /= Physics->sumOfWeightsCells	[iCell];
#endif
#if (DARCY)
		Physics->DeltaP0	[iCell] /= Physics->sumOfWeightsCells	[iCell];
		Physics->phi0		[iCell] /= Physics->sumOfWeightsCells	[iCell];
#endif
#if (STORE_PLASTIC_STRAIN)
		Physics->strain		[iCell] /= Physics->sumOfWeightsCells	[iCell];
#endif

	}





	// Filling side values
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->sigma_xx_0, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);
#if (HEAT)
	Physics_CellVal_SideValues_getFromBC(Physics->T, Grid, BCThermal, NumThermal);
#endif
#if (DARCY)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->DeltaP0, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->phi0, Grid);
#endif
#if (STORE_PLASTIC_STRAIN)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->strain, Grid);
#endif



#if (HEAT)
	// Should probably be moved to a specific function
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->k[iCell] = 0.0;
		thisPhaseInfo = Physics->phaseListHead[iCell];
		while (thisPhaseInfo != NULL) {
			Physics->k[iCell] += MatProps->k[thisPhaseInfo->phase] * thisPhaseInfo->weight;
			thisPhaseInfo = thisPhaseInfo->next;
		}
		Physics->k[iCell] /= Physics->sumOfWeightsCells[iCell];
	}
#endif

	// ==================================
	// Interpolate to nodes
	// ==================================
//	int signX;
//	int signY;

	int iNodeNeigh;
	xMod[0] =  1; yMod[0] =  1;
	xMod[1] =  0; yMod[1] =  1;
	xMod[2] =  1; yMod[2] =  0;
	xMod[3] =  0; yMod[3] =  0;


	//int iColor; // indexing of the color group for nodes. Nodes of the same color don't collide with each other. i.e. similar to matrix coloring
	int ixStartS[9] = {0,0,0,1,1,1,2,2,2};
	int iyStartS[9] = {0,1,2,0,1,2,0,1,2};
	//SinglePhase* thisPhaseInfo;


	for (iColor = 0; iColor < 9; ++iColor) {
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY, phase, i, iNodeNeigh, weight) OMP_SCHEDULE
		for (iy = iyStartS[iColor]; iy < Grid->nyS; iy+=3) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStartS[iColor]; ix < Grid->nxS; ix+=3) { // I don't get why though
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];

				// Loop through the particles in the shifted cell
				// ======================================
				while (thisParticle!=NULL) {
					locX = Particles_getLocX(ix, thisParticle->x,Grid);
					locY = Particles_getLocY(iy, thisParticle->y,Grid);
/*
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
*/
					locX = fabs(locX);
					locY = fabs(locY);

					phase = thisParticle->phase;
#if (PART2GRID_SCHEME == 1)
					
					for (i=0; i<4; i++) {
						iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

						if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
							printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
							//printf("IX = %i, IY = %i, locX = %.2e, locY = %.2e, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
							printf("thisParticle->x = %.2e , y = %.2e \n", thisParticle->x, thisParticle->y);
							exit(0);
						}

						

						weight = (fabs(locX) + xMod[i])   *   (fabs(locY) + yMod[i]);
#else
						weight = (1.0 - fabs(locX)) * (1.0 - fabs(locY));
						iNodeNeigh = iNode;
#endif

						
#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)

						Physics->sigma_xy_0 		[iNodeNeigh] += thisParticle->sigma_xy_0 * weight;
#endif
						Physics->GShear 			[iNodeNeigh] += weight / MatProps->G[phase];// * weight;
						Physics->sumOfWeightsNodes	[iNodeNeigh] += weight; // using the same arrays
#if (PART2GRID_SCHEME == 1)
					}
#endif
					thisParticle = thisParticle->next;

				}

			}
		}
	}


	// Adding contribution to the other side for periodic BC
	if(Grid->isPeriodic) {
		int iCellS, iCellD;
#pragma omp parallel for private(iy, iCellS, iCellD,i) OMP_SCHEDULE
		for (iy = 0; iy < Grid->nyS; ++iy) {
			iCellS = 0 + iy*Grid->nxS; // Source
			iCellD = Grid->nxS-1 + iy*Grid->nxS; // destination
			Physics->sigma_xy_0 		[iCellD] += Physics->sigma_xy_0  [iCellS];
			Physics->sigma_xy_0 		[iCellS]  = Physics->sigma_xy_0  [iCellD];
			Physics->GShear		 		[iCellD] += Physics->GShear  [iCellS];
			Physics->GShear		 		[iCellS]  = Physics->GShear  [iCellD];
			Physics->sumOfWeightsNodes	[iCellD] += Physics->sumOfWeightsNodes[iCellS];
			Physics->sumOfWeightsNodes	[iCellS]  = Physics->sumOfWeightsNodes[iCellD];
		}
	}



#if (TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
	// Dividing by the sum of weights
#pragma omp parallel for private(iNode) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {

		Physics->sigma_xy_0 [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic caverage
		Physics->GShear [iNode]  = Physics->sumOfWeightsNodes[iNode] / Physics->GShear[iNode]; // arithmetic caverage
		//Physics->sigma_xy_0 [iNode] = Physics->sumOfWeightsNodes[iNode] / Physics->sigma_xy_0 [iNode]; // harmonic average
	}

#endif


#if (!TEST_SIGMA_INTERP_FROM_PART_TO_CELL)
	//Physics_CellVal_advectEulerian(Physics->sigma_xx_0, Model);
	//Physics_NodeVal_advectEulerian(Physics->sigma_xy_0, Model);
#endif



}
































#if (HEAT)
void Interp_Temperature_Grid2Particles_Global(Model* Model)
{
Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps, BC* BCThermal

Grid* Grid 				= &(Model->Grid);
Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	MatProps* MatProps 		= &(Model->MatProps);
	BC* BCThermal 			= &(Model->BCThermal);





	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

	compute* DT_sub_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* sumOfWeights_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* DT_rem_OnTheCells = (compute*) malloc( Grid->nECTot *sizeof(compute) );

	int i, iCell;

#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i=0;i<4*Grid->nECTot;++i) {
		DT_sub_OnTheCells[i] = 0;
		sumOfWeights_OnTheCells[i] = 0;
	}

	compute TFromCells, DT_sub_OnThisPart, PFromCells;
	compute rhoParticle;
	compute dtDiff;
	compute d = 0.98;

	int phase;

	compute weight;
	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right

	int xModCell[4], yModCell[4];
	xModCell[0] = -1; yModCell[0] = -1;
	xModCell[1] =  1; yModCell[1] = -1;
	xModCell[2] = -1; yModCell[2] =  1;
	xModCell[3] =  1; yModCell[3] =  1;





	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->DT[i] = (Physics->T[i] - Physics->T0[i]) * Physics->dtAdv/Physics->dtT;
	}

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, phase, locX, locY, TFromCells, PFromCells, rhoParticle, dtDiff, DT_sub_OnThisPart, i, iCell, weight) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			phase = thisParticle->phase;
			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);


				
				TFromCells  = 		Interp_Any_Cell2Particle_Local(Physics->T, ix, iy, Grid->nxEC, locX, locY);
				PFromCells  = 		Interp_Any_Cell2Particle_Local(Physics->P, ix, iy, Grid->nxEC, locX, locY);


				rhoParticle = MatProps->rho0[phase];// * (1+MatProps->beta[phase]*PFromCells) * (1-MatProps->alpha[phase]*thisParticle->T);
				if (rhoParticle<0) {
					printf("error: Negative density on Particles in Physisc_interTempFromCellsParticle\n");
					exit(0);
				}

				dtDiff = (Physics->Cp*rhoParticle)/(  MatProps->k[phase]*( 2.0/(Grid->dx*Grid->dx) + 2.0/(Grid->dy*Grid->dy) )  );


				DT_sub_OnThisPart = ( TFromCells - thisParticle->T ) * ( 1.0 - exp(-d * Physics->dtAdv/dtDiff) );

				// redefine locX, locY (used to compute surface based weight, not used as weight directly)
				locX = (thisParticle->x-Grid->xmin)/dx - ix;
				locY = (thisParticle->y-Grid->ymin)/dy - iy;
				// Interp Dsigma_xx_sub from particles to Cells
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
					weight = fabs((locX + xModCell[i]*0.5)   *   (locY + yModCell[i]*0.5));

					DT_sub_OnTheCells[iCell*4+i] += DT_sub_OnThisPart * weight;
					sumOfWeights_OnTheCells	[iCell*4+i] += weight;

				}

				thisParticle->T += DT_sub_OnThisPart;

				thisParticle = thisParticle->next;
			}
		}
	}


	compute DT_sub_OnThisCell;
	compute sum;
	int I;
#pragma omp parallel for private(iCell, I, sum, DT_sub_OnThisCell) OMP_SCHEDULE
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		I = 4*iCell;
		sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];
		if (sum==0) {
			printf("error in Physics_interpFromParticlesToCell: cell #%i received no contribution from particles\n", iCell );
			exit(0);
		}

		DT_sub_OnThisCell = ( DT_sub_OnTheCells[I+0] + DT_sub_OnTheCells[I+1] + DT_sub_OnTheCells[I+2] + DT_sub_OnTheCells[I+3]) / sum ; // harmonic average
		DT_rem_OnTheCells[iCell] = Physics->DT[iCell] - DT_sub_OnThisCell;

	}

	Physics_CellVal_SideValues_copyNeighbours_Global(DT_rem_OnTheCells, Grid);



	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);
				
				thisParticle->T  += Interp_ECVal_Cell2Particle_Local(DT_rem_OnTheCells, ix, iy, Grid->nxEC, locX, locY);

				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->T = BCThermal->TT;
				}

				thisParticle = thisParticle->next;
			}
		}
	}



	free(DT_sub_OnTheCells);
	free(sumOfWeights_OnTheCells);
	free(DT_rem_OnTheCells);

}
#endif
















#if (DARCY)
void Interp_Phi_Grid2Particles_Global(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;

#if (DARCY)
		int i;
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DDeltaP[i] *=  Physics->dtAdv/Physics->dt;
			Physics->Dphi[i] *=  Physics->dtAdv/Physics->dt;
			//printf("Physics->dtAdv/Physics->d = %.2e\n",Physics->dtAdv/Physics->dt);
		}
#endif

	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);



				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {

					thisParticle->DeltaP0 = 0.0;


				} else {

				
				thisParticle->DeltaP0 += Interp_ECVal_Cell2Particle_Local(Physics->DDeltaP, ix, iy, Grid->nxEC, locX, locY);
				}
				thisParticle->phi += Interp_ECVal_Cell2Particle_Local(Physics->Dphi, ix, iy, Grid->nxEC, locX, locY);



				thisParticle = thisParticle->next;
			}
		}
	}


}
#endif


#if (STORE_PLASTIC_STRAIN)
void Interp_Strain_Grid2Particles_Global(Model* Model)
{

Grid* Grid 				= &(Model->Grid);
Particles* Particles 	= &(Model->Particles);
Physics* Physics 		= &(Model->Physics);
	

	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;


	int iCell;
	compute SII;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			if (Physics->khi[iCell]<1e30 && Physics->phase[iCell] != Physics->phaseAir && Physics->phase[iCell] != Physics->phaseWater) {
				SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);// //(Physics, Grid, ix, iy, &SII);
				Physics->Dstrain[iCell] = SII/(2.0*Physics->khi[iCell])*Physics->dtAdv; // Recovering the incremental plastic strain
			} else {
				Physics->Dstrain[iCell] = 0.0;
			}
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Dstrain, Grid);


	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {
				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);

				thisParticle->strain += Interp_ECVal_Cell2Particle_Local(Physics->Dstrain, ix, iy, Grid->nxEC, locX, locY);
				thisParticle = thisParticle->next;
			}
		}
	}


}
#endif




#if (TEST_SIGMA_INTERP)

void Interp_Stresses_Grid2Particles_Global(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);


	compute locX, locY;
	int ix, iy;

	INIT_PARTICLE
	compute Dsigma_xx_0_Grid;
	compute Dsigma_xy_0_Grid;

	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);


#if (USE_SPECIAL_STRESS_INTERP)
					Dsigma_xx_0_Grid  = Interp_Special_Sxx_Cell2Particle_Local(Physics->Dsigma_xx_0, Physics->Lambda , ix, iy, Grid->nxEC, locX, locY);
					Dsigma_xy_0_Grid = Interp_Special_Sxy_Node2Particle_Local(Physics->Dsigma_xy_0, Physics->LambdaShear, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
#else
					Dsigma_xx_0_Grid = Interp_ECVal_Cell2Particle_Local(Physics->Dsigma_xx_0, ix, iy, Grid->nxEC, locX, locY);
					Dsigma_xy_0_Grid = Interp_NodeVal_Node2Particle_Local(Physics->Dsigma_xy_0, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
#endif			

					
					if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
						thisParticle->sigma_xx_0 = 0.0;
						thisParticle->sigma_xy_0 = 0.0;
					} else {
						thisParticle->sigma_xx_0 += Dsigma_xx_0_Grid;
						thisParticle->sigma_xy_0 += Dsigma_xy_0_Grid;
					}
				thisParticle = thisParticle->next;
			}
		}
	}


}





#else
void Interp_Stresses_Grid2Particles_Global(Model* Model)
{
	// General
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);



	INIT_PARTICLE

	int ix, iy;

	compute locX, locY;


	int signX, signY;


	// See Taras' book pp. 186-187
	compute sigma_xx_0_fromCells;
	compute sigma_xy_0_fromNodes;

	compute d_ve = Numerics->stressSubGridDiffFac;
	compute dtm = Physics->dtAdv;
	compute dtMaxwell;

	compute eta, G;

	compute* Dsigma_xy_sub_OnTheNodes = (compute*) malloc( 4*Grid->nSTot  *sizeof(compute) );
	compute* sumOfWeights_OnTheNodes  = (compute*) malloc( 4*Grid->nSTot  *sizeof(compute) );
	compute* Dsigma_xx_sub_OnTheCells = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );
	compute* sumOfWeights_OnTheCells  = (compute*) malloc( 4*Grid->nECTot *sizeof(compute) );


	compute* Dsigma_xy_rem_OnTheNodes = (compute*) malloc( Grid->nSTot  *sizeof(compute) );
	compute* Dsigma_xx_rem_OnTheCells = (compute*) malloc( Grid->nECTot *sizeof(compute) );

	int i;
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i=0;i<4*Grid->nSTot;++i) {
		Dsigma_xy_sub_OnTheNodes[i] = 0.0;
		sumOfWeights_OnTheNodes[i] = 0.0;
	}
#pragma omp parallel for private(i) OMP_SCHEDULE
	for (i=0;i<4*Grid->nECTot;++i) {
		Dsigma_xx_sub_OnTheCells[i] = 0.0;
		sumOfWeights_OnTheCells[i] = 0.0;
	}

	int iNodeNeigh;
	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right


	int iCell;
	compute xModNode[4], yModNode[4], xModCell[4], yModCell[4];
	
	compute weight;
	xModNode[0] =  1.0; yModNode[0] =  1.0;
	xModNode[1] =  0.0; yModNode[1] =  1.0;
	xModNode[2] =  1.0; yModNode[2] =  0.0;
	xModNode[3] =  0.0; yModNode[3] =  0.0;


	xModCell[0] = -1.0; yModCell[0] = -1.0;
	xModCell[1] =  1.0; yModCell[1] = -1.0;
	xModCell[2] = -1.0; yModCell[2] =  1.0;
	xModCell[3] =  1.0; yModCell[3] =  1.0;

	compute Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart;

	compute khi, eta_vp;

	//printf("dtM/dt = %.2e, ( 1.0 - exp(-d_ve * dtm/dtMaxwell) =%.2e\n", Numerics->subgridStressDiffTimeScale/Physics->dtAdv, ( 1.0 - exp(-d_ve * Physics->dtAdv/Numerics->subgridStressDiffTimeScale) ) );

	// ===============================================================================================
	// Stress rotation stuff


	// Compute the Alpha array
	// add a condi	ztion with signX signY to avoid recomputing alpha if not necessary
	compute* alphaArray = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* Rotxy = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* dVxdyGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* dVydxGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute alpha;
#pragma omp parallel for private(iy, ix, iNode) OMP_SCHEDULE
	for (iy=0; iy<Grid->nyS; iy++) {
		for (ix=0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			
#if (USE_UPPER_CONVECTED)
			alphaArray[iNode]  =  0.0;

			compute dVxdy, dVydx;

			//
			
			/*
			if (iy>0 && iy<Grid->nyS-1 && ix>0 && ix<Grid->nxS-1) {
				dVxdy = 2.0/3.0 * ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
				dVxdy += 1.0/3.0 * 0.5*((VxCell[ix   + (iy+1)*Grid->nxEC] - VxCell[ix   +(iy  )*Grid->nxEC])/Grid->DYEC[iy]
									+(VxCell[ix+1 + (iy+1)*Grid->nxEC] - VxCell[ix+1 +(iy  )*Grid->nxEC])/Grid->DYEC[iy]);
				//
				dVydx = 2.0/3.0 * ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
				dVydx += 1.0/3.0 *  0.5*((VyCell[ix+1 + (iy  )*Grid->nxEC] - VyCell[ix   +(iy  )*Grid->nxEC])/Grid->DXEC[ix]
										+(VyCell[ix+1 + (iy+1)*Grid->nxEC] - VyCell[ix   +(iy+1)*Grid->nxEC])/Grid->DXEC[ix]);
			} else {
				dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
				dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			}
			*/
			

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;

			dVxdyGrid[iNode] =  dVxdy;
			dVydxGrid[iNode] =  dVydx;

			Rotxy[iNode] = 0.5*(dVxdy-dVydx);
#else
/*
			alphaArray[iNode]  = 2.0/3.0 * .5*Physics->dtAdv*((Physics->Vy[ix+1 + (iy  )*Grid->nxVy] - Physics->Vy[ix   +(iy  )*Grid->nxVy])/Grid->DXEC[ix]
						- (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   +(iy  )*Grid->nxVx])/Grid->DYEC[iy]);

			alphaArray[iNode]  += 1.0/3.0 *  .5*Physics->dtAdv*( 0.5*((VyCell[ix+1 + (iy  )*Grid->nxEC] - VyCell[ix   +(iy  )*Grid->nxEC])/Grid->DXEC[ix]
									+(VyCell[ix+1 + (iy+1)*Grid->nxEC] - VyCell[ix   +(iy+1)*Grid->nxEC])/Grid->DXEC[ix])
								- 0.5*((VxCell[ix   + (iy+1)*Grid->nxEC] - VxCell[ix   +(iy  )*Grid->nxEC])/Grid->DYEC[iy]
									+(VxCell[ix+1 + (iy+1)*Grid->nxEC] - VxCell[ix+1 +(iy  )*Grid->nxEC])/Grid->DYEC[iy]));
*/
			alphaArray[iNode]  = 3.0/3.0 * .5*Physics->dtAdv*((Physics->Vy[ix+1 + (iy  )*Grid->nxVy] - Physics->Vy[ix   +(iy  )*Grid->nxVy])/Grid->DXEC[ix]
														- (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   +(iy  )*Grid->nxVx])/Grid->DYEC[iy]);
	
			
			//alphaArray[iNode]  =  0.0;

			

			
			
			

#endif
		}
	}
	
	compute* Exx = (compute*) malloc(Grid->nECTot * sizeof(compute));
	
	compute dVxdx, dVydy;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						- Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						- Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Exx[iCell] = 0.5*(dVxdx-dVydy);


		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx, Grid);






	// Stress rotation stuff
	// ===============================================================================================
















	// compute Dsigma_xx_0_sub on the particles and interpolate to the grid
//#pragma omp parallel for private(iy, ix, i, iNode, thisParticle, locX, locY, signX, signY, sigma_xx_0_fromCells, sigma_xy_0_fromNodes, eta, khi, G, eta_vp, dtMaxwell, Dsigma_xx_sub_OnThisPart, Dsigma_xy_sub_OnThisPart, iNodeNeigh, weight, iCell) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			compute EII;
			dtMaxwell = Numerics->subgridStressDiffTimeScale;
			
			compute eta = Physics->etaShear[iNode];
			compute khi = Physics->khiShear[iNode];
			compute G = Physics->GShear[iNode];
			Physics_StrainRateInvariant_getLocalNode(Model, ix, iy, &EII);
			compute VP_EP = (1.0/(1.0/(eta) + 1.0/khi)) / (1.0/(1.0/G + Physics->dt/khi));
			compute fAngle = MatProps->frictionAngle[thisParticle->phase];
			compute coh = MatProps->frictionAngle[thisParticle->phase];
			compute P = Physics->P[ix+iy*Grid->nxEC];
			P = fmin(P,Physics->P[ix+1+iy*Grid->nxEC]);
			P = fmin(P,Physics->P[ix+1+(iy+1)*Grid->nxEC]);
			P = fmin(P,Physics->P[ix+(iy+1)*Grid->nxEC]);


			int Count = 0;
			
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);

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
				
				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {

					// Compute Dsigma sub grid
					Dsigma_xx_sub_OnThisPart =   0.0;
					Dsigma_xy_sub_OnThisPart =   0.0;

					// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
					thisParticle->sigma_xx_0 = 0.0;
					thisParticle->sigma_xy_0 = 0.0;

				}
				else {
					sigma_xx_0_fromCells  = Interp_ECVal_Cell2Particle_Local(Physics->sigma_xx_0, ix, iy, Grid->nxEC, locX, locY);
					sigma_xy_0_fromNodes = Interp_NodeVal_Node2Particle_Local(Physics->sigma_xy_0, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
			
					
					
					if (P<0.0) {
						P = 0.0;
					}
					compute Ty = P * sin(fAngle) + coh * cos(fAngle);
					compute refTime_noPlast;
					if ((2.0*eta*EII - Ty )>0.0) {
						refTime_noPlast = eta/G* log(2.0*eta*EII / (2.0*eta*EII - Ty ));
					} else {
						refTime_noPlast = 1e100;
					}

					dtMaxwell = fmin(VP_EP,refTime_noPlast);
					
					//dtMaxwell = refTime_noPlast;
					/*
					if (ix == Grid->nxS-3) {
						if (Count==0) {
							//printf("dt/dtM = %.2e, ( 1.0 - exp(-d_ve * dtm/dtMaxwell) =%.2e, dt = %.2e, VP_EP = %.2e, refTime_noPlast = %.2e\n", Physics->dtAdv/dtMaxwell, ( 1.0 - exp(-d_ve * Physics->dtAdv/dtMaxwell) ) , Physics->dt, VP_EP, refTime_noPlast);
							//printf("dt/dtM = %.2e, ( 1.0 - exp(-d_ve * dtm/dtMaxwell) =%.2e, dt = %.2e, VP_EP =%.2e, refTime = %.2e\n", Physics->dtAdv/dtMaxwell, ( 1.0 - exp(-d_ve * Physics->dtAdv/dtMaxwell) ) , Physics->dt, VP_EP, refTime_noPlast);
							//printf("dt/dtM = %.2e, ( 1.0 - exp(-d_ve * dtm/dtMaxwell) =%.2e, dt = %.2e, refTime = %.2e\n", Physics->dtAdv/dtMaxwell, ( 1.0 - exp(-d_ve * Physics->dtAdv/dtMaxwell) ) , Physics->dt, refTime_noPlast);
							printf("refTime = %.2e\n", refTime_noPlast);
							Count++;
						}
					}
					*/
					
								
					// Compute Dsigma sub grid
					Dsigma_xx_sub_OnThisPart = ( sigma_xx_0_fromCells - thisParticle->sigma_xx_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
					Dsigma_xy_sub_OnThisPart = ( sigma_xy_0_fromNodes - thisParticle->sigma_xy_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
					/*
					if ( ( 1.0 - exp(-d_ve * dtm/dtMaxwell)<0.0) || ( 1.0 - exp(-d_ve * dtm/dtMaxwell)>1.0) || isnan(Dsigma_xx_sub_OnThisPart) ) {
						printf("Problem with Fac: Fac = %.2e, eta_vp = %.2e, eta = %.2e, khi = %.2e, dtm =%.2e, dtMaxwell = %.2e, eta/Physics->G[iCell]=%.2e, log(2.0*eta*EII / (2.0*eta*EII - Ty )) = %.2e, Ty =%.2e, EII = %.2e, (2.0*eta*EII - Ty ) = %.2e\n", ( 1.0 - exp(-d_ve * dtm/dtMaxwell) ) , eta_vp, eta, khi, dtm, dtMaxwell, eta/G, log(2.0*eta*EII / (2.0*eta*EII - Ty )), Ty, EII, (2.0*eta*EII - Ty ));
						exit(0);
					}
					*/
					

#if (USE_UPPER_CONVECTED)
				compute Sxx0 = thisParticle->sigma_xx_0; 
				compute Sxy0 = thisParticle->sigma_xy_0;
				compute Z = Interp_ECVal_Cell2Particle_Local(Physics->Z, ix, iy, Grid->nxEC, locX, locY);
				compute G = Interp_ECVal_Cell2Particle_Local(Physics->G, ix, iy, Grid->nxEC, locX, locY);

				compute ExxPart = Interp_ECVal_Cell2Particle_Local(Exx, ix, iy, Grid->nxEC, locX, locY);
				compute dVxdyPart = Interp_NodeVal_Node2Particle_Local(dVxdyGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				compute dVydxPart = Interp_NodeVal_Node2Particle_Local(dVydxGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				thisParticle->sigma_xx_0 +=  - Z/G*( - 2.0*Sxx0*ExxPart   - 2.0*Sxy0*dVxdyPart);
				thisParticle->sigma_xy_0 +=  - Z/G*( - 1.0*Sxx0*dVydxPart + 1.0*Sxx0*dVxdyPart);
				
#endif

					// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
					thisParticle->sigma_xx_0 += Dsigma_xx_sub_OnThisPart;
					thisParticle->sigma_xy_0 += Dsigma_xy_sub_OnThisPart;
				}





				// Interp Dsigma_xy_sub from particles to Nodes
#if (PART2GRID_SCHEME == 1)
				for (i=0; i<4; i++) {

					iNodeNeigh = ix+IxN[i]*signX  +  (iy+IyN[i]*signY)*Grid->nxS;

					if (ix+IxN[i]*signX>Grid->nxS || ix+IxN[i]*signX<0 || (iy+IyN[i]*signY)>Grid->nyS || (iy+IyN[i]*signY)<0) {
						printf("error in interpFromParticlesToCells: trying to access a non existing node\n");
						//printf("IX = %i, IY = %i, locX = %.2e, locY = %.2e, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i, Counter = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX, Counter);
						printf("thisParticle->x = %.2e , y = %.2e \n", thisParticle->x, thisParticle->y);
						exit(0);
					}


					//locX = fabs(locX);
					//locY = fabs(locY);

					weight = (fabs(locX) + xModNode[i])   *   (fabs(locY) + yModNode[i]);
					//weight = fabs((locX + xMod[i])   *   (locY + yMod[i]));

					Dsigma_xy_sub_OnTheNodes[iNodeNeigh*4+i] += Dsigma_xy_sub_OnThisPart * weight;
					sumOfWeights_OnTheNodes [iNodeNeigh*4+i] += weight; // using the same arrays
				}
				
				// Interp Dsigma_xx_sub from particles to Cells
				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
					weight = fabs((locX + xModCell[i]*0.5)   *   (locY + yModCell[i]*0.5));

					Dsigma_xx_sub_OnTheCells[iCell*4+i] += Dsigma_xx_sub_OnThisPart * weight;
					sumOfWeights_OnTheCells	[iCell*4+i] += weight;

				}
				

#else
				weight = (1.0 - fabs(locX)) * (1.0 - fabs(locY));
				Dsigma_xy_sub_OnTheNodes[iNode*4] += Dsigma_xy_sub_OnThisPart * weight;
				sumOfWeights_OnTheNodes [iNode*4] += weight; // using the same arrays

				// Interp Dsigma_xx_sub from particles to Cells
				int signX, signY;
				if (locX<0.0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0.0) {
					signY = -1;
				} else {
					signY = 1;
				}
				if 		 	(signX>=0 && signY>=0) { // upper right
					i = 3;
				} else if 	(signX<0 && signY>=0) { // upper left
					// the particle is in the SE quadrant, the cell center 1 is NW (wrt to the node ix,iy)
					i = 2;
				} else if 	(signX>=0 && signY<0) { // lower right
					i = 1;
				} else if 	(signX<0 && signY<0) { // lower left
					i = 0;
				} else {
					printf("error in Interp_ECVal_Cell2Particle_Local. No case was triggered\n.");
					exit(0);
				}
				//printf("i = %i, locX = %.2f, locY = %.2f\n",i,locX,locY);
				iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
				weight = fabs(locX)*fabs(locY);
				Dsigma_xx_sub_OnTheCells[iCell*4+i] += Dsigma_xx_sub_OnThisPart * weight;
				sumOfWeights_OnTheCells	[iCell*4+i] += weight;
#endif


				thisParticle = thisParticle->next;

			}
		}
	}



	int I;
	compute sum;
	compute Dsigma_xy_sub_OnThisNode;
#pragma omp parallel for private(iNode, I, sum, Dsigma_xy_sub_OnThisNode) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		I = 4*iNode;
		sum = sumOfWeights_OnTheNodes[I+0] + sumOfWeights_OnTheNodes[I+1] + sumOfWeights_OnTheNodes[I+2] + sumOfWeights_OnTheNodes[I+3];
		if (sum==0) {
			printf("error in Physics_interpStressesFromCellsToParticles: node #%i received no contribution from particles\n", iNode );
			exit(0);
		}

		Dsigma_xy_sub_OnThisNode = ( Dsigma_xy_sub_OnTheNodes[I+0] +  Dsigma_xy_sub_OnTheNodes[I+1] +  Dsigma_xy_sub_OnTheNodes[I+2] +  Dsigma_xy_sub_OnTheNodes[I+3]) / sum ; // harmonic average

		Dsigma_xy_rem_OnTheNodes[iNode] = Physics->Dsigma_xy_0[iNode] - Dsigma_xy_sub_OnThisNode;
		if (isnan(Dsigma_xy_rem_OnTheNodes[iNode])) {
			printf("Dsigma_xy_rem_OnTheNodes[iNode] is nan, Physics->Dsigma_xy_0[iNode] = %.2e, Dsigma_xx_sub_OnThisCell = %.2e\n",Physics->Dsigma_xy_0[iNode],Dsigma_xy_sub_OnThisNode);
		}
	}



#if (PART2GRID_SCHEME == 0)
	// For this scheme, outer cells receive no contribution from particles
	// Not so important because calculation is not made on them
	// But to avoid division by 0, I here copy the values from the neighbours anyway.
	// Also this allows to check for empty cells.
	int phase;
	SinglePhase* thisPhaseInfo;
	int nxEC = Grid->nxEC;
	for (iy=1;iy<Grid->nyEC-1;iy++) {
		for (ix=1;ix<Grid->nxEC-1;ix++) {
			iCell = ix + iy*Grid->nxEC;
			I = 4*iCell;
			sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];
			
			if (sum == 0.0) {
				printf("diffSum = %.2e\n", fabs(sum-Physics->sumOfWeightsCells[iCell]));
				//exit(0);
				printf("Trying to save something!\n");
				// If no contributions was given to this cell (i.e. empty cell), then use a higher order interpolation scheme (2x2 cells wide instead of 1x1)
				int iNodeCounter = 0;
				for (iNodeCounter=0;iNodeCounter<4;iNodeCounter++) {

					iNode = ix  + (iy  )*Grid->nxS;
					thisParticle = Particles->linkHead[iNode];
					// Loop through the particles in the shifted cell
					// ======================================
					while (thisParticle!=NULL) {
						locX = Particles_getLocX(ix, thisParticle->x,Grid);
						locY = Particles_getLocY(iy, thisParticle->y,Grid);


						sigma_xx_0_fromCells  = Interp_ECVal_Cell2Particle_Local(Physics->sigma_xx_0, ix, iy, Grid->nxEC, locX, locY);
						
						eta  				  = Interp_ECVal_Cell2Particle_Local(Physics->eta, ix, iy, Grid->nxEC, locX, locY);
						khi  				  = Interp_ECVal_Cell2Particle_Local(Physics->khi, ix, iy, Grid->nxEC, locX, locY);
						eta_vp = 1.0 / (1.0/eta + 1.0/khi);
						//eta_vp = fmax(eta_vp,Numerics->etaMin);
						G = MatProps->G[thisParticle->phase];
						dtMaxwell = eta_vp/G;




						// Compute Dsigma sub grid
						Dsigma_xx_sub_OnThisPart = ( sigma_xx_0_fromCells - thisParticle->sigma_xx_0 ) * ( 1.0 - exp(-d_ve * dtm/dtMaxwell) );
						thisParticle->sigma_xx_0 += Dsigma_xx_sub_OnThisPart;
					
/*					
#if (USE_UPPER_CONVECTED)
				compute Sxx0 = thisParticle->sigma_xx_0; 
				compute Sxy0 = thisParticle->sigma_xy_0;
				compute Z = Interp_ECVal_Cell2Particle_Local(Physics->Z, ix, iy, Grid->nxEC, locX, locY);
				compute G = Interp_ECVal_Cell2Particle_Local(Physics->G, ix, iy, Grid->nxEC, locX, locY);

				compute ExxPart = Interp_ECVal_Cell2Particle_Local(Exx, ix, iy, Grid->nxEC, locX, locY);
				compute dVxdyPart = Interp_NodeVal_Node2Particle_Local(dVxdyGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				compute dVydxPart = Interp_NodeVal_Node2Particle_Local(dVydxGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				thisParticle->sigma_xx_0 +=  - Z/G*( - 2.0*Sxx0*ExxPart   - 2.0*Sxy0*dVxdyPart);
				thisParticle->sigma_xy_0 +=  - Z/G*( - 1.0*Sxx0*dVydxPart + 1.0*Sxx0*dVxdyPart);
				
#endif
*/

						for (i=0; i<4; i++) {
							int thisCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
							if (thisCell==iCell) {
								weight = fabs((locX + xModCell[i])   *   (locY + yModCell[i]));
								// Get the phase and weight of phase contribution for each cell
								
								Dsigma_xx_sub_OnTheCells[iCell*4+i] += Dsigma_xx_sub_OnThisPart * weight;
								sumOfWeights_OnTheCells	[iCell*4+i] += weight;

							} // if thisCell=iCell
						} // for neighbour cells
						thisParticle = thisParticle->next;
					} // while Particles
				} // iNodeCounter
			} //if (Physics->sumOfWeightsCells	[iCell] == 0) 
		} // ixCell
	} // iyCell
	
#endif



	compute Dsigma_xx_sub_OnThisCell;
#pragma omp parallel for private(iy, ix, iCell, I, sum, Dsigma_xx_sub_OnThisCell) OMP_SCHEDULE
	for (iy=1;iy<Grid->nyEC-1;iy++) {
		for (ix=1;ix<Grid->nxEC-1;ix++) {
			iCell = ix +iy*Grid->nxEC;
			I = 4*iCell;
			sum = sumOfWeights_OnTheCells[I+0] + sumOfWeights_OnTheCells[I+1] + sumOfWeights_OnTheCells[I+2] + sumOfWeights_OnTheCells[I+3];


			Dsigma_xx_sub_OnThisCell = ( Dsigma_xx_sub_OnTheCells[I+0] + Dsigma_xx_sub_OnTheCells[I+1] + Dsigma_xx_sub_OnTheCells[I+2] + Dsigma_xx_sub_OnTheCells[I+3]) / sum ; // harmonic average
			Dsigma_xx_rem_OnTheCells[iCell] = Physics->Dsigma_xx_0[iCell] - Dsigma_xx_sub_OnThisCell;
			if (isnan(Dsigma_xx_rem_OnTheCells[iCell])) {
				printf("Dsigma_xx_rem_OnTheCells[iCell] is nan, Physics->Dsigma_xx_0[iCell] = %.2e, Dsigma_xx_sub_OnThisCell = %.2e, sum = %.2e\n",Physics->Dsigma_xx_0[iCell],Dsigma_xx_sub_OnThisCell, sum);
			}
		}
	}

	// Copy values to sides
	Physics_CellVal_SideValues_copyNeighbours_Global(Dsigma_xx_rem_OnTheCells, Grid);







	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, signX, signY) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);


				// Stress rotation for upper convected (uses old stresses)




				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->sigma_xx_0 = 0.0;
					thisParticle->sigma_xy_0 = 0.0;
				} else {
					thisParticle->sigma_xx_0  += Interp_ECVal_Cell2Particle_Local(Dsigma_xx_rem_OnTheCells, ix, iy, Grid->nxEC, locX, locY);
					thisParticle->sigma_xy_0  += Interp_NodeVal_Node2Particle_Local(Dsigma_xy_rem_OnTheNodes, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				}
				if (isnan(thisParticle->sigma_xx_0)) {
					printf("Sxx on particle is nan\n");
				}
				if (isnan(thisParticle->sigma_xy_0)) {
					printf("Sxy on particle is nan\n");
				}

#if (!USE_UPPER_CONVECTED)
				// Rotation of stresses without assuming a small angle
				compute alpha = Interp_NodeVal_Node2Particle_Local(alphaArray, ix, iy, Grid->nxS, Grid->nyS, locX, locY);				
				compute sigma_xx_temp = thisParticle->sigma_xx_0*cos(alpha)*cos(alpha) - thisParticle->sigma_xx_0*sin(alpha)*sin(alpha)  -  thisParticle->sigma_xy_0*sin(2.0*alpha);
				thisParticle->sigma_xy_0 = thisParticle->sigma_xy_0*cos(2.0*alpha)  +  thisParticle->sigma_xx_0*sin(2.0*alpha);
				thisParticle->sigma_xx_0 = sigma_xx_temp;
#endif



				thisParticle = thisParticle->next;
			}
		}
	}


	free(Dsigma_xy_sub_OnTheNodes);
	free(Dsigma_xx_sub_OnTheCells);

	free(sumOfWeights_OnTheNodes);
	free(sumOfWeights_OnTheCells);


	free(Dsigma_xy_rem_OnTheNodes);
	free(Dsigma_xx_rem_OnTheCells);


	free(alphaArray);



	free(Exx);
	free(Rotxy);
	free(dVxdyGrid);
	free(dVydxGrid);


}



#endif
