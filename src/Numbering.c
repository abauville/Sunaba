/*
 * numbering.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Numbering_allocateMemory(Numbering* Numbering, EqSystem* EqSystem, Grid* Grid) {

		Numbering->map  		= (int*) 		malloc(EqSystem->nEqIni 	* sizeof(int)); // Numbering map
		Numbering->IX   = (int*) malloc( EqSystem->nEq     * sizeof(int)); // contains ix indices for all equations without dirichlet (the numbering starts  from 0 for Vx, Vy and P equations)
		Numbering->IY   = (int*) malloc( EqSystem->nEq     * sizeof(int));
}

void Numbering_freeMemory(Numbering* Numbering)
{

	free( Numbering->map );
	free( Numbering->IX );
	free( Numbering->IY );

}




void Numbering_init(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering, Physics* Physics)
{

	//==========================================================================
	//
	//                  INIT Numbering->map, EqSystem->I and nNonzeros
	//
	//==========================================================================
	// Numbering->map stores the number of all equations on the grid in a continuous manner.
	// Dirichlet equations are not numbered and set as -1 in Numbering->map

	int I = 0;
	int InoDir = 0;
	int sum = 0;

	bool jumping; // to jump nodes, for periodic BC

	int iSubEqSystem;

	int nx, ny;


	// For the local stencil
	StencilType thisStencil;
	int Jloc[11];
	compute Vloc[11];
	compute bloc;
	int order[11] = {0,1,2,3,4,5,6,7,8,9,10};
	int nLoc = 0;
	int shift = 0;
	int i;
	int SetupType = BC->SetupType;

	int IC = 0;






	int ix, iy;


	EqSystem->I[0] = 0;

	// Init Numbering->IX, Numbering->IY with -1
	for (i=0;i<EqSystem->nEq;i++){
		Numbering->IX[i] = -1;
		Numbering->IY[i] = -1;
	}
	// 1. fill with ones, except first value = 0
	for (i=0;i<EqSystem->nEqIni;i++){
		Numbering->map[i] = 1;
	}

	// 2. replace value for Dirichlet and Neumann equations by 0
	for (i=0; i<BC->n; i++) { // Velocity Dirichlet
		Numbering->map[ BC->list[i] ] = 0;
	}






	// compute the number of non zeros using the 0 and 1 filled Numbering->map
	// + Fill EqSystem->I
	// + Fill Numbering->IX and Numbering->IY vectors





	EqSystem->nnz = 0;

	Numbering->subEqSystem0[0] = 0;
	for (iSubEqSystem=0;iSubEqSystem<Numbering->nSubEqSystem;iSubEqSystem++) {
		thisStencil = Numbering->Stencil[iSubEqSystem];


		switch (thisStencil) {
		case Stencil_Stokes_Momentum_x:
			nx = Grid->nxVx;
			ny = Grid->nyVx;
			break;

		case Stencil_Stokes_Momentum_y:
			nx = Grid->nxVy;
			ny = Grid->nyVy;
			break;

		case Stencil_Stokes_Darcy_Continuity:
		case Stencil_Stokes_Continuity:
			nx = Grid->nxC;
			ny = Grid->nyC;
			break;

		case Stencil_Poisson:
		case Stencil_Heat:
			nx = Grid->nxC+2;
			ny = Grid->nyC+2;
			Numbering->map[0] = 0;
			Numbering->map[nx-1] = 0;
			Numbering->map[nx*(ny-1)] = 0;
			Numbering->map[nx*(ny)-1] = 0;
			break;
		default:
			printf("error: unknwon Stencil %i", thisStencil);
			exit(0);
			break;
		}


		for (iy=0; iy<ny; iy++)
		{
			for (ix=0; ix<nx; ix++)
			{

				if (Numbering->map[I] != 0) // Free equation, i.e. neither a Dirichlet equation nor Neumann
				{



					// Check if this node should be jumped
					jumping = false;
					switch (thisStencil) {
					case Stencil_Stokes_Momentum_x:
						if ((BC->SetupType==SimpleShearPeriodic  && ix==nx-1) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}
						break;
					case Stencil_Stokes_Momentum_y:
						if ( (BC->SetupType==SimpleShearPeriodic && ix>=nx-2) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}
						break;
					case Stencil_Stokes_Continuity:
						break;
					case Stencil_Heat:
						/*
						// jump corners
						if ( (ix==0) && (iy=0) ) // To jump the rightmost nodes for periodic bc
							jumping = true;
						if ( (ix==nx+2) && (iy=0) ) // To jump the rightmost nodes for periodic bc
							jumping = true;
						if ( (ix==0) && (iy=ny+2) ) // To jump the rightmost nodes for periodic bc
							jumping = true;
						if ( (ix==nx+2) && (iy=ny+2) ) // To jump the rightmost nodes for periodic bc
							jumping = true;
						 */
						if ( (BC->SetupType==SimpleShearPeriodic && ix>=nx-2) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}

						break;
					default:
						printf("error: unknwon Stencil %i", thisStencil);
						exit(0);
						break;
					}




					// Fill I, IX and IY
					if (!jumping) {
						// Get the nnz
						//	Numbering_getLocalNNZ(ix, iy, Numbering, Grid, BC, true, thisStencil, &sum, Physics);
						LocalStencil_Call(thisStencil, order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);

						sum = 0;
						for (i = shift; i < nLoc; ++i) {
							sum += Numbering->map[ Jloc[i] ];
						}
						if (sum==0) {
							sum = 1; // For compatibility with Pardiso
						}
						EqSystem->nnz += sum;
						EqSystem->I[InoDir+1] = EqSystem->nnz;
						Numbering->IX[InoDir] = ix;
						Numbering->IY[InoDir] = iy;
						InoDir++;
					}



				}

				I++;
			}
		}
		Numbering->subEqSystem0[iSubEqSystem+1] = InoDir;


	}



	printf("InoDir = %i nEq = %i, nRow = %i\n",InoDir, EqSystem->nEq, EqSystem->nRow);

	// 3. Numbering->map = cumsum(Numbering->map);
	// Start numbering at 0, but jump the first dirichlet nodes
	i = 0;
	while (Numbering->map[i]==0){
		i++;
	}
	Numbering->map[i] = 0;




	if (BC->SetupType==SimpleShearPeriodic) // Number the Equations on the Right boundary with the number from the left one
	{

		int Ileft;
		i = 0;
		int iPrevious = 0;
		bool replaceNode;
		bool replaceSecondNode;

		for (iSubEqSystem=0;iSubEqSystem<Numbering->nSubEqSystem;iSubEqSystem++) {
			thisStencil = Numbering->Stencil[iSubEqSystem];
			switch (thisStencil) {
			case Stencil_Stokes_Momentum_x:
				nx = Grid->nxVx-1;
				ny = Grid->nyVx;
				replaceNode = true;
				replaceSecondNode = false;
				break;
			case Stencil_Stokes_Momentum_y:
				nx = Grid->nxVy-2;
				ny = Grid->nyVy;
				replaceNode = true;
				replaceSecondNode = true;
				break;
			case Stencil_Stokes_Continuity:
				nx = Grid->nxC;
				ny = Grid->nyC;
				replaceNode = false;
				replaceSecondNode = false;
				break;
			case Stencil_Heat:
				nx = Grid->nxC+2-2;
				ny = Grid->nyC+2;
				replaceNode = true;
				replaceSecondNode = true;
				break;
			default:
				printf("error: unknwon Stencil %i", thisStencil);
				exit(0);
				break;
			}


			for (iy=0; iy<ny; iy++) {
				Ileft = i;
				for (ix=0; ix<nx; ix++) {
					Numbering->map[i] += Numbering->map[iPrevious];
					i++;
					iPrevious = i-1;
				}
				if (replaceNode) {
					Numbering->map[i] = Numbering->map[Ileft];
					i++;
				}
				if (replaceSecondNode) {
					Numbering->map[i] = Numbering->map[Ileft+1];
					i++;
				}
			}
		}
	}

	else
	{
		for (i=1;i<EqSystem->nEqIni;i++){
			Numbering->map[i] += Numbering->map[i-1];
		}
	}


	// 4. replace value for Dirichlet equations by -1
	I = -1;
	for (i=0; i<BC->n; i++) { // Velocity Dirichlet
		Numbering->map[ BC->list[i] ] = I;
		I--;
	}




	if (DEBUG) {
		printf("===== Numbering->map =====\n");  printListi  (Numbering->map       ,EqSystem->nEqIni);
		printf("===== EqSystem->I =====\n");  printListi  (EqSystem->I       ,EqSystem->nEq+1);
		printf("===== IX =====\n");  printListi  (Numbering->IX       ,EqSystem->nEq);
		printf("===== IY =====\n");  printListi  (Numbering->IY       ,EqSystem->nEq);
	}




	if (DEBUG) {
		int C = 0;
		for (iSubEqSystem=0;iSubEqSystem<Numbering->nSubEqSystem;iSubEqSystem++) {
			thisStencil = Numbering->Stencil[iSubEqSystem];
			switch (thisStencil) {
			case Stencil_Stokes_Momentum_x:
				nx = Grid->nxVx;
				ny = Grid->nyVx;
				break;
			case Stencil_Stokes_Momentum_y:
				nx = Grid->nxVy;
				ny = Grid->nyVy;
				break;
			case Stencil_Stokes_Continuity:
				nx = Grid->nxC;
				ny = Grid->nyC;
				break;
			case Stencil_Heat:
				nx = Grid->nxC+2;
				ny = Grid->nyC+2;
				break;
			default:
				printf("error: unknwon Stencil %i", thisStencil);
				exit(0);
				break;
			}

			printf("=== numMap Stencil %i\n", thisStencil);
			for (iy = 0; iy < ny; ++iy) {
				for (ix = 0; ix < nx; ++ix) {
					printf("%i  ",Numbering->map[C]);
					C++;
				}
				printf("\n");
			}
		}



	}




}


void Numbering_getLocalNNZ(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, bool useNumMap, StencilType StencilType, int* sum, Physics* Physics)
{

	int Jloc[11];
	compute Vloc[11];
	compute bloc;
	int order[11] = {0,1,2,3,4,5,6,7,8,9,10};
	int nLoc = 0;
	int shift = 0;
	int i;
	int SetupType = BC->SetupType;

	int IC = 0;


	if (StencilType == Stencil_Stokes_Momentum_x) {
		LocalStencil_Stokes_Momentum_x(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);
	}

	else if (StencilType==Stencil_Stokes_Momentum_y) {
		LocalStencil_Stokes_Momentum_y(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);
	}

	else if (StencilType == Stencil_Stokes_Continuity) {
		LocalStencil_Stokes_Continuity(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);
	}

	else if (StencilType == Stencil_Heat) {
		LocalStencil_Heat(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);
	}

	else {
		printf("error: unknwon stencil %i", StencilType);
		exit(0);
	}


	*sum = 0;
	for (i = shift; i < nLoc; ++i) {
		*sum += Numbering->map[ Jloc[i] ];
	}
	if (*sum==0) {
		*sum = 1; // For compatibility with Pardiso
	}

}





