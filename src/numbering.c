/*
 * numbering.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"


void Numbering_init(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering)
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
	StencilType thisStencil;


	EqSystem->I 	= (int*) malloc((EqSystem->nEq+1)  * sizeof(int));
	Numbering->IX   = (int*) malloc( EqSystem->nEq     * sizeof(int)); // contains ix indices for all equations without dirichlet (the numbering starts  from 0 for Vx, Vy and P equations)
	Numbering->IY   = (int*) malloc( EqSystem->nEq     * sizeof(int));

	int ix, iy;

	int i;
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
		case Vx:
			nx = Grid->nxVx;
			ny = Grid->nyVx;
			break;
		case Vy:
			nx = Grid->nxVy;
			ny = Grid->nyVy;
			break;
		case P:
			nx = Grid->nxC;
			ny = Grid->nyC;
			break;
		case T:
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
					// Get the nnz
					Numbering_getLocalNNZ(ix, iy, Numbering, Grid, BC, true, thisStencil, &sum);


					// Check if this node should be jumped
					jumping = false;
					switch (thisStencil) {
					case Vx:
						if ((BC->SetupType==SimpleShearPeriodic  && ix==nx-1) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}
						break;
					case Vy:
						if ( (BC->SetupType==SimpleShearPeriodic && ix>=nx-2) ) { // To jump the rightmost nodes for periodic bc
							jumping = true;
						}
						break;
					case P:
						if (UPPER_TRI)
							sum = 1; // particular condition for Pardiso, otherwise in UPPER_TRI it would be 0
						break;
					case T:
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
			case Vx:
				nx = Grid->nxVx-1;
				ny = Grid->nyVx;
				replaceNode = true;
				replaceSecondNode = false;
				break;
			case Vy:
				nx = Grid->nxVy-2;
				ny = Grid->nyVy;
				replaceNode = true;
				replaceSecondNode = true;
				break;
			case P:
				nx = Grid->nxC;
				ny = Grid->nyC;
				replaceNode = false;
				replaceSecondNode = false;
				break;
			case T:
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
			case Vx:
				nx = Grid->nxVx;
				ny = Grid->nyVx;
				break;
			case Vy:
				nx = Grid->nxVy;
				ny = Grid->nyVy;
				break;
			case P:
				nx = Grid->nxC;
				ny = Grid->nyC;
				break;
			case T:
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


void Numbering_getLocalNNZ(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, bool useNumMap, StencilType StencilType, int* sum)
{

	if (StencilType == Vx) {

		LocalNumberingVx LocVx;

		int nxVx 	= Grid->nxVx;
		int nxVy 	= Grid->nxVy;
		int nVxTot 	= Grid->nVxTot;
		int nVyTot 	= Grid->nVyTot;
		int nxC 	= Grid->nxC;

		LocVx.VxC     = ix      + iy*nxVx                 			;
		LocVx.VxN     = ix      + iy*nxVx     + nxVx      			;
		LocVx.VxS     = ix      + iy*nxVx     - nxVx      			;
		LocVx.VxE     = ix      + iy*nxVx     + 1         			;
		LocVx.VxW     = ix      + iy*nxVx     - 1        			;
		LocVx.VyNE    = ix+1    + iy*nxVy     + nVxTot   			;
		LocVx.VyNW    = ix+0    + iy*nxVy     + nVxTot   			;
		LocVx.VySE    = ix+1    + (iy-1)*nxVy + nVxTot   			;
		LocVx.VySW    = ix+0    + (iy-1)*nxVy + nVxTot   			;

		LocVx.PE      = ix      + (iy-1)*nxC  + nVxTot + nVyTot	;
		LocVx.PW      = ix-1    + (iy-1)*nxC  + nVxTot + nVyTot	;

		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				LocVx.VxW += nxVx-1;
				LocVx.PW  += nxC;
			}
		}



		if (useNumMap) {
			LocVx.VxC     = Numbering->map[ LocVx.VxC ];
			LocVx.VxN     = Numbering->map[ LocVx.VxN ];
			LocVx.VxS     = Numbering->map[ LocVx.VxS ];
			LocVx.VxE     = Numbering->map[ LocVx.VxE ];
			LocVx.VxW     = Numbering->map[ LocVx.VxW ];
			LocVx.VyNE    = Numbering->map[ LocVx.VyNE];
			LocVx.VyNW    = Numbering->map[ LocVx.VyNW];
			LocVx.VySE    = Numbering->map[ LocVx.VySE];
			LocVx.VySW    = Numbering->map[ LocVx.VySW];

			LocVx.PE      = Numbering->map[ LocVx.PE];
			LocVx.PW      = Numbering->map[ LocVx.PW];


			if (UPPER_TRI) {
				LocVx.VxS = 0;
				LocVx.VxW = 0;

				if (BC->SetupType==SimpleShearPeriodic) {
					if (ix==0) {
						LocVx.VxW = 1;
					}
					if (ix==nxVx-2) {
						LocVx.VxE = 0;
					}
				}
			}
		}

		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix==nxVx-1) {
				LocVx.VxC     = 0	;
				LocVx.VxN     = 0 			;
				LocVx.VxS     = 0		;
				LocVx.VxE     = 0  			;
				LocVx.VxW     = 0   			;
				LocVx.VyNE    = 0 			;
				LocVx.VyNW    = 0			;
				LocVx.VySE    = 0		;
				LocVx.VySW    = 0			;

				LocVx.PE      = 0	;
				LocVx.PW      = 0	;
			}
		}


		*sum = LocVx.VxS + LocVx.VxW + LocVx.VxC + LocVx.VxE + LocVx.VxN
				+ LocVx.VySW + LocVx.VySE + LocVx.VyNW + LocVx.VyNE
				+ LocVx.PW + LocVx.PE; // VxS etc... are 1 if the equation they refer to is free, 0 otherwise

	}

	else if (StencilType==Vy) {
		LocalNumberingVy LocVy;
		int nxVx 	= Grid->nxVx;
		int nxVy 	= Grid->nxVy;
		int nVxTot 	= Grid->nVxTot;
		int nVyTot 	= Grid->nVyTot;
		int nxC 	= Grid->nxC;


		LocVy.VyC     = ix     + iy*nxVy + nVxTot          	;
		LocVy.VyN     = ix     + iy*nxVy + nVxTot   + nxVy     ;
		LocVy.VyS     = ix     + iy*nxVy + nVxTot   - nxVy     ;
		LocVy.VyE     = ix     + iy*nxVy + nVxTot   + 1        ;
		LocVy.VyW     = ix     + iy*nxVy + nVxTot   - 1        ;
		LocVy.VxNE    = ix     + (iy+1)*nxVx               	;
		LocVy.VxNW    = ix     + (iy+1)*nxVx - 1          		;
		LocVy.VxSE    = ix     + (iy  )*nxVx               	;
		LocVy.VxSW    = ix     + (iy  )*nxVx - 1           	;

		LocVy.PN      = ix-1   + (iy  )*nxC + nVxTot + nVyTot 	;
		LocVy.PS      = ix-1   + (iy-1)*nxC + nVxTot + nVyTot 	;


		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				LocVy.VxSW += nxVx-1;
				LocVy.VxNW += nxVx-1;
				LocVy.VyW  += nxVy-2;
				LocVy.PS   += nxC;
				LocVy.PN   += nxC;
			}
		}



		if (useNumMap) {
			LocVy.VyC     = Numbering->map[ LocVy.VyC      	];
			LocVy.VyN     = Numbering->map[ LocVy.VyN  ];
			LocVy.VyS     = Numbering->map[ LocVy.VyS     ];
			LocVy.VyE     = Numbering->map[ LocVy.VyE       ];
			LocVy.VyW     = Numbering->map[ LocVy.VyW     ];
			LocVy.VxNE    = Numbering->map[ LocVy.VxNE       	];
			LocVy.VxNW    = Numbering->map[ LocVy.VxNW          		];
			LocVy.VxSE    = Numbering->map[ LocVy.VxSE     	];
			LocVy.VxSW    = Numbering->map[ LocVy.VxSW        	];

			LocVy.PN      = Numbering->map[ LocVy.PN 	];
			LocVy.PS      = Numbering->map[ LocVy.PS	];

			if (UPPER_TRI) {
				LocVy.VyS = 0;
				LocVy.VyW = 0;
				LocVy.VxNE = 0;
				LocVy.VxNW = 0;
				LocVy.VxSE = 0;
				LocVy.VxSW = 0;

				if (BC->SetupType==SimpleShearPeriodic) {
					if (ix==0) {
						LocVy.VyW = 1;
					}
					else if (ix==nxVy-3) {
						LocVy.VyE = 0;
					}
				}
			}
		}



		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix>=nxVy-2) {
				LocVy.VyC     = 0;
				LocVy.VyN     = 0;
				LocVy.VyS     = 0;
				LocVy.VyE     = 0;
				LocVy.VyW     = 0;
				LocVy.VxNE    = 0;
				LocVy.VxNW    = 0;
				LocVy.VxSE    = 0;
				LocVy.VxSW    = 0;

				LocVy.PN      = 0;
				LocVy.PS      = 0;
			}
		}

		*sum =  LocVy.VxSW + LocVy.VxSE + LocVy.VxNW + LocVy.VxNE
				+ LocVy.VyS + LocVy.VyW + LocVy.VyC + LocVy.VyE + LocVy.VyN
				+ LocVy.PS + LocVy.PN; // VxSW etc... are 1 if the equation they refer to is free, 0 otherwise

	}

	else if (StencilType == P) {
		LocalNumberingP LocP;

		int nxVx 	= Grid->nxVx;
		int nxVy 	= Grid->nxVy;
		int nVxTot 	= Grid->nVxTot;


		LocP.VxE     = ix+1 + (iy+1)*nxVx               ;
		LocP.VxW     = ix   + (iy+1)*nxVx               ;
		LocP.VyN     = ix+1 + (iy+1)*(nxVy)  + nVxTot   ;
		LocP.VyS     = ix+1 + iy*(nxVy)      + nVxTot	;

		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				LocP.VxW += nxVx-1;
			}
		}

		if (useNumMap) {
			LocP.VxE     = Numbering->map[ ix+1 + (iy+1)*nxVx               ];
			LocP.VxW     = Numbering->map[ ix   + (iy+1)*nxVx               ];
			LocP.VyN     = Numbering->map[ ix+1 + (iy+1)*(nxVy)  + nVxTot   ];
			LocP.VyS     = Numbering->map[ ix+1 + iy*(nxVy)      + nVxTot	];
			if (UPPER_TRI) {
				LocP.VxE     = 0;
				LocP.VxW     = 0;
				LocP.VyN     = 0;
				LocP.VyS     = 0;
			}
		}
		*sum = LocP.VxW + LocP.VxE + LocP.VyS + LocP.VyN;

	}


	else if (StencilType == T) {
		int TS, TW, TC, TE, TN;

		TS =  ix 	+ (iy-1)*(Grid->nxC+2);
		TW = (ix-1) +  iy   *(Grid->nxC+2);
		TC =  ix 	+  iy   *(Grid->nxC+2);
		TE = (ix+1) +  iy   *(Grid->nxC+2);
		TN =  ix 	+ (iy+1)*(Grid->nxC+2);

		if (BC->SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				TW  += (Grid->nxC+2)-2  ; // VyW
			}
		}
		if (useNumMap) {
			TS     = Numbering->map[  ix 	+ (iy-1)*(Grid->nxC+2)         ];
			TW     = Numbering->map[  (ix-1) +  iy   *(Grid->nxC+2)      ];
			TC     = Numbering->map[ ix 	+  iy   *(Grid->nxC+2)  ];
			TE     = Numbering->map[ (ix+1) +  iy   *(Grid->nxC+2)	];
			TN     = Numbering->map[ ix 	+ (iy+1)*(Grid->nxC+2)	];


			if (UPPER_TRI) {
				TS = 0;
				TW = 0;


				if (BC->SetupType==SimpleShearPeriodic) {
					if (ix==0) {
						TW = 1;
					}
					else if (ix==(Grid->nxC+2)-3) {
						TE = 0;
					}
					if (ix>=(Grid->nxC+2)-2) {
						TS     = 0;
						TW    = 0;
						TC    = 0;
						TE     = 0;
						TN    = 0;

					}
				}
			}
		}

		*sum = TW + TC + TE + TS + TN;

	}

	else {
		printf("error: unknwon stencil %i", StencilType);
		exit(0);
	}

}





