/*
 * numbering.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Numbering_Memory_allocate(Numbering* Numbering, EqSystem* EqSystem, Grid* Grid) {

		Numbering->map  	= (int*) 		malloc(EqSystem->nEqIni 	* sizeof(int)); // Numbering map
		Numbering->IX   = (int*) malloc( EqSystem->nEq     * sizeof(int)); // contains ix indices for all equations without dirichlet (the numbering starts  from 0 for Vx, Vy and P equations)
		Numbering->IY   = (int*) malloc( EqSystem->nEq     * sizeof(int));
}

void Numbering_Memory_free(Numbering* Numbering)
{

	free( Numbering->map );
	free( Numbering->IX );
	free( Numbering->IY );

}




void Numbering_init(Model* Model, EqSystemType EqSystemType)
{

	//==========================================================================
	//
	//                  INIT Numbering->map, EqSystem->I and nNonzeros
	//
	//==========================================================================
	// Numbering->map stores the number of all equations on the grid in a continuous manner.
	// Dirichlet equations are not numbered and set as -1 in Numbering->map
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);

	BC* BC;
	EqSystem* EqSystem;
	Numbering* Numbering;

	if (EqSystemType == EqSystemType_Stokes) {
		BC			= &(Model->BCStokes);
		EqSystem 	= &(Model->EqStokes);
		Numbering 	= &(Model->NumStokes);
	} else if (EqSystemType == EqSystemType_Thermal) {
		BC			= &(Model->BCThermal);
		EqSystem 	= &(Model->EqThermal);
		Numbering 	= &(Model->NumThermal);
	} else {
		printf("Error: Unknwon EqSystemType %i", EqSystemType);
		exit(0);
	}

	int I = 0;
	int InoDir = 0;
	int sum = 0;

	bool jumping; // to jump nodes, for periodic BC

	int iSubEqSystem;

	int nx, ny;


	// For the local stencil
	StencilType thisStencil;
	int Jloc[13];
	compute Vloc[13];
	compute bloc;
	int order[13] = {0,1,2,3,4,5,6,7,8,9,10,11,12};
	int nLoc = 0;
	int shift = 0;
	int i;

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
	Numbering->subEqSystem0Dir[0] = 0;
	for (iSubEqSystem=0;iSubEqSystem<Numbering->nSubEqSystem;iSubEqSystem++) {
		thisStencil = Numbering->Stencil[iSubEqSystem];


		switch (thisStencil) {
		case Stencil_Stokes_Darcy_Momentum_x:
		case Stencil_Stokes_Momentum_x:
			nx = Grid->nxVx;
			ny = Grid->nyVx;
			break;

		case Stencil_Stokes_Darcy_Momentum_y:
		case Stencil_Stokes_Momentum_y:
			nx = Grid->nxVy;
			ny = Grid->nyVy;
			break;

		case Stencil_Stokes_Darcy_Darcy:
		case Stencil_Stokes_Darcy_Continuity:
		case Stencil_Stokes_Continuity:
		case Stencil_Poisson:
		case Stencil_Heat:
			nx = Grid->nxEC;
			ny = Grid->nyEC;

			// Not sure why the following is important
			Numbering->map[I+0] = 0;
			Numbering->map[I+nx-1] = 0;
			Numbering->map[I+nx*(ny-1)] = 0;
			Numbering->map[I+nx*(ny)-1] = 0;
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
					int ixShift = 0;


					// Check if this node should be jumped
					jumping = false;
					switch (thisStencil) {
					case Stencil_Stokes_Darcy_Momentum_x:
					case Stencil_Stokes_Momentum_x:
						if ((Grid->isPeriodic  && ix==nx-1) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}
						break;
					case Stencil_Stokes_Darcy_Momentum_y:
					case Stencil_Stokes_Momentum_y:
						if (Grid->isPeriodic){
							if (ix>=nx-2)  { // To jump the rightmost nodes for periodic bc
								jumping = true;
							} else if (ix==0) {
								ixShift = nx-2;
							}
						}
						break;
					case Stencil_Stokes_Darcy_Continuity:
					case Stencil_Stokes_Darcy_Darcy:
					case Stencil_Stokes_Continuity:
					case Stencil_Poisson:
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
						if (Grid->isPeriodic){
							if (ix>=nx-2)  { // To jump the rightmost nodes for periodic bc
								jumping = true;
							} else if (ix==0) {
								ixShift = nx-2;
							}
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
						LocalStencil_Call(Model, thisStencil, order, Jloc, Vloc, &bloc, ix+ixShift, iy, &shift, &nLoc, &IC);

						sum = 0;
						for (i = shift; i < nLoc; ++i) {
							sum += Numbering->map[ Jloc[i] ];
						}
						if (sum==0) {
							sum = 1; // For compatibility with Pardiso
						}
						EqSystem->nnz += sum;
						EqSystem->I[InoDir+1] = EqSystem->nnz;
						Numbering->IX[InoDir] = ix+ixShift;
						Numbering->IY[InoDir] = iy;
						InoDir++;
					}



				}

				I++;
			}
		}
		Numbering->subEqSystem0[iSubEqSystem+1] = InoDir;
		Numbering->subEqSystem0Dir[iSubEqSystem+1] = I;

	}



	printf("InoDir = %i nEq = %i, nRow = %i\n",InoDir, EqSystem->nEq, EqSystem->nRow);



	// 3. Numbering->map = cumsum(Numbering->map);
	// Start numbering at 0, but jump the first dirichlet nodes
	i = 0;
	while (Numbering->map[i]==0){
		i++;
	}
	Numbering->map[i] = 0;






	if (Grid->isPeriodic) // Number the Equations on the Right boundary with the number from the left one
	{

		int Ileft;
		i = 0;
		int iPrevious = 0;
		bool replaceNode;
		bool replaceSecondNode;

		for (iSubEqSystem=0;iSubEqSystem<Numbering->nSubEqSystem;iSubEqSystem++) {
			thisStencil = Numbering->Stencil[iSubEqSystem];
			switch (thisStencil) {
			case Stencil_Stokes_Darcy_Momentum_x:
			case Stencil_Stokes_Momentum_x:
				nx = Grid->nxVx-1;
				ny = Grid->nyVx;
				replaceNode = true;
				replaceSecondNode = false;
				break;
			case Stencil_Stokes_Darcy_Momentum_y:
			case Stencil_Stokes_Momentum_y:
				nx = Grid->nxVy-2;
				ny = Grid->nyVy;
				replaceNode = true;
				replaceSecondNode = true;
				break;
			case Stencil_Stokes_Darcy_Continuity:
			case Stencil_Stokes_Darcy_Darcy:
			case Stencil_Stokes_Continuity:
			case Stencil_Poisson:
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
			case Stencil_Stokes_Darcy_Momentum_x:
			case Stencil_Stokes_Momentum_x:
				nx = Grid->nxVx;
				ny = Grid->nyVx;
				break;
			case Stencil_Stokes_Darcy_Momentum_y:
			case Stencil_Stokes_Momentum_y:
				nx = Grid->nxVy;
				ny = Grid->nyVy;
				break;
			case Stencil_Stokes_Darcy_Continuity:
			case Stencil_Stokes_Darcy_Darcy:
			case Stencil_Stokes_Continuity:
			case Stencil_Poisson:
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



