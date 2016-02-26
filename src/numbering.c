/*
 * numbering.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"


void Numbering_initMapAndSparseTripletIJ(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering)
{

	//==========================================================================
	//
	//                  INIT Numbering->map, EqSystem->I and nNonzeros
	//
	//==========================================================================
	// Numbering->map stores the number of all equations on the grid in a continuous manner.
	// Dirichlet equations are not numbered and set as -1 in Numbering->map
	EqSystem->I 	= (int*) malloc((EqSystem->nEq+1)  * sizeof(int));
    Numbering->IX   = (int*) malloc( EqSystem->nEq     * sizeof(int)); // contains ix indices for all equations without dirichlet (the numbering starts  from 0 for Vx, Vy and P equations)
    Numbering->IY   = (int*) malloc( EqSystem->nEq     * sizeof(int));

	LocalNumberingVx LocVx;
	LocalNumberingVy LocVy;
	LocalNumberingP LocP ;
	int ix, iy;

	int i;
	EqSystem->I[0] = 0;
	if (UPPER_TRI) {
		EqSystem->nRow = EqSystem->nEq + BC->nPDir - Grid->nCTot;
	}
	else {
		EqSystem->nRow = EqSystem->nEq;
	}
	// Init Numbering->IX, Numbering->IY with -1
	for (i=0;i<EqSystem->nEq;i++){
		Numbering->IX[i] = -1;
		Numbering->IY[i] = -1;
	}
	// 1. fill with ones, except first value = 0
	for (i=0;i<EqSystem->nEqIni;i++){
		Numbering->map[i] = 1;
	}
	// 2. replace value for Dirichlet equations by 0
	for (i=0; i<BC->nDir; i++) { // Velocity Dirichlet
		Numbering->map[ BC->listDir[i] ] = 0;
	}


	// compute the number of non zeros using the 0 and 1 filled Numbering->map
	// + Fill EqSystem->I
	// + Fill Numbering->IX and Numbering->IY vectors

	// Vx
	EqSystem->nnz = 0;
	int I = 0;
	int InoDir = 0;
	int iNeu;
	//int VxEq0 = InoDir;

	for (iy=0; iy<Grid->nyVx; iy++)
	{
		for (ix=0; ix<Grid->nxVx; ix++)
		{

			if (Numbering->map[I] != 0) // If not a Dirichlet equation
			{
				if (BC->isNeu[I]) // If Neumann equation
				{
					if (UPPER_TRI) {
						iNeu = 0;
						while (BC->listNeu[iNeu] != I) {
							iNeu++;
							if (iNeu==BC->nNeu) {
								fprintf(stderr,"Couldn't find the Neumann equations");
								exit(0);
							}
						}
						if (BC->listNeu[iNeu]>BC->listNeuNeigh[iNeu]) {
							EqSystem->nnz += 1;
						}
						else {
							EqSystem->nnz += 2;
						}
					}
					else {
						EqSystem->nnz += 2;
					}
				}
				else { // if free equation (i.e. not Neumann, not Dir)
					Numbering_getLocalVx(ix, iy, Numbering, Grid, &LocVx, true);

					if (UPPER_TRI) {
						EqSystem->nnz += LocVx.VxC + LocVx.VxE + LocVx.VxN
								   + LocVx.VySW + LocVx.VySE + LocVx.VyNW + LocVx.VyNE
								   + LocVx.PW + LocVx.PE; // VxS etc... are 1 if the equation they refer to is free, 0 otherwise
					}
					else {
						EqSystem->nnz += LocVx.VxS + LocVx.VxW + LocVx.VxC + LocVx.VxE + LocVx.VxN
								   + LocVx.VySW + LocVx.VySE + LocVx.VyNW + LocVx.VyNE
								   + LocVx.PW + LocVx.PE; // VxS etc... are 1 if the equation they refer to is free, 0 otherwise
					}
				}

				EqSystem->I[InoDir+1] = EqSystem->nnz;
				Numbering->IX[InoDir] = ix;
				Numbering->IY[InoDir] = iy;
				InoDir++;
			}

			I++;
		}
	}
	// Vy
	printf("C\n");
	EqSystem->VyEq0 = InoDir;

	for (iy=0; iy<Grid->nyVy; iy++)
	{
		for (ix=0; ix<Grid->nxVy; ix++)
		{
			if (Numbering->map[I] != 0) // If not a Dirichlet equation
			{
				if (BC->isNeu[I]) // If Neumann equation
				{
					if (UPPER_TRI)
					{
						iNeu = 0;
						while (BC->listNeu[iNeu] != I)
						{
							iNeu++;
							if (iNeu==BC->nNeu) {
								fprintf(stderr,"Couldn't find the Neumann equations");
								exit(0);
							}
						}

						if (BC->listNeu[iNeu]>BC->listNeuNeigh[iNeu]) {
							EqSystem->nnz += 1;
						}
						else {
							EqSystem->nnz += 2;
						}
					}
					else {
						EqSystem->nnz += 2;
					}

				}
				else { // if free equation (i.e. not Neumann, not Dir)
					Numbering_getLocalVy(ix, iy, Numbering, Grid, &LocVy, true);
					if (UPPER_TRI) {
						EqSystem->nnz +=  LocVy.VyC + LocVy.VyE + LocVy.VyN + LocVy.PS + LocVy.PN; // VxSW etc... are 1 if the equation they refer to is free, 0 otherwise
					}
					else {
						EqSystem->nnz +=  LocVy.VxSW + LocVy.VxSE + LocVy.VxNW + LocVy.VxNE
								    + LocVy.VyS + LocVy.VyW + LocVy.VyC + LocVy.VyE + LocVy.VyN
									+ LocVy.PS + LocVy.PN; // VxSW etc... are 1 if the equation they refer to is free, 0 otherwise
					}

				}
				EqSystem->I[InoDir+1] = EqSystem->nnz;
				Numbering->IX[InoDir] = ix;
				Numbering->IY[InoDir] = iy;
				InoDir++;
			}
			I++;
		}
	}



	EqSystem->PEq0 = InoDir;


	for (iy=0; iy<Grid->nyC; iy++)
	{
		for (ix=0; ix<Grid->nxC; ix++)
		{


			if (Numbering->map[I] != 0) // If not a Dirichlet equation
			{

				if (UPPER_TRI) {
					EqSystem->nnz+=1;
					EqSystem->I[InoDir+1] = EqSystem->nnz;
				}
				else {

					Numbering_getLocalP(ix, iy, Numbering, Grid, &LocP, true);
					EqSystem->nnz += LocP.VxW + LocP.VxE + LocP.VyS + LocP.VyN;

					EqSystem->I[InoDir+1] = EqSystem->nnz;


				}

				Numbering->IX[InoDir] = ix;
				Numbering->IY[InoDir] = iy;
				InoDir++;
			}
			I++;

		}
	}


	// 3. Numbering->map = cumsum(Numbering->map);
	i = 0;
	while (Numbering->map[i]==0){
		i++;
	}
	Numbering->map[i] = 0;
	for (i=1;i<EqSystem->nEqIni;i++){
		Numbering->map[i] += Numbering->map[i-1];
	}

	// 4. replace value for Dirichlet equations by -1
	for (i=0; i<BC->nDir; i++) { // Velocity Dirichlet
		Numbering->map[ BC->listDir[i] ] = -1;
	}


	// Apply Numbering->map to BC->listNeu
	for (i=0; i<BC->nNeu; i++) {
		BC->listNeu[i] = Numbering->map[ BC->listNeu[i] ];
		BC->listNeuNeigh[i] = Numbering->map[ BC->listNeuNeigh[i] ];

	}


	if (DEBUG) {
		printf("===== Numbering->map =====\n");  printListi  (Numbering->map       ,EqSystem->nEqIni);
		printf("===== EqSystem->I =====\n");  printListi  (EqSystem->I       ,EqSystem->nRow+1);
	}

	//Check BC indices
	// ==============
	if (DEBUG) {
		printf("===== BC INDICES =====\n");
		printf("BC->listNeu        : ");    printListi(BC->listNeu       ,BC->nNeu);
		printf("BC->listNeuNeigh  : ");    printListi(BC->listNeuNeigh ,BC->nNeu);
		printf("BC->coeffNeu       : ");    printListd(BC->coeffNeu      ,BC->nNeu);
		printf("BC->coeffNeuNeigh : ");    printListd(BC->coeffNeuNeigh,BC->nNeu);
		printf("BC->valueNeu       : ");    printListd(BC->valueNeu      ,BC->nNeu);
		printf("\n");

		//printf("BC->listDir            : ");    printListi(BC->listDir           ,nDir)
		//printf("BC->valueDir           : ");    printListd(BC->valueDir          ,nDir)
		//printf("\n");
	}


	if (DEBUG) {

		int ix, iy;
		int C = 0;
		printf("=== numMap Vx\n");
		for (iy = 0; iy < Grid->nyVx; ++iy) {
			for (ix = 0; ix < Grid->nxVx; ++ix) {
				printf("%i  ",Numbering->map[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== numMap Vy\n");
		for (iy = 0; iy < Grid->nyVy; ++iy) {
			for (ix = 0; ix < Grid->nxVy; ++ix) {
				printf("%i  ",Numbering->map[C]);
				C++;
			}
			printf("\n");
		}

		printf("=== numMap P\n");
		for (iy = 0; iy < Grid->nyC; ++iy) {
			for (ix = 0; ix < Grid->nxC; ++ix) {
				printf("%i  ",Numbering->map[C]);
				C++;
			}
			printf("\n");
		}

	}




}


void Numbering_getLocalVx(int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingVx* LocVx, bool useNumMap)
{


	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;
	int nVyTot 	= Grid->nVyTot;
	int nxC 	= Grid->nxC;
	int nxS 	= Grid->nxS;
	if (useNumMap) {
		LocVx->VxC     = Numbering->map[ ix      + iy*nxVx                 			];
		LocVx->VxN     = Numbering->map[ ix      + iy*nxVx     + nxVx           	];
		LocVx->VxS     = Numbering->map[ ix      + iy*nxVx     - nxVx           	];
		LocVx->VxE     = Numbering->map[ ix      + iy*nxVx     + 1             	 	];
		LocVx->VxW     = Numbering->map[ ix      + iy*nxVx     - 1             		];
		LocVx->VyNE    = Numbering->map[ ix+1    + iy*nxVy     + nVxTot  			];
		LocVx->VyNW    = Numbering->map[ ix+0    + iy*nxVy     + nVxTot    			];
		LocVx->VySE    = Numbering->map[ ix+1    + (iy-1)*nxVy + nVxTot    			];
		LocVx->VySW    = Numbering->map[ ix+0    + (iy-1)*nxVy + nVxTot    			];
		LocVx->NormalE = Numbering->map[ ix      + (iy-1)*nxC              			];
		LocVx->NormalW = Numbering->map[ ix-1    + (iy-1)*nxC               		];
		LocVx->ShearN  = Numbering->map[ ix      + iy*nxS                   		];
		LocVx->ShearS  = Numbering->map[ ix      + (iy-1)*nxS               		];
		LocVx->PE      = Numbering->map[ ix      + (iy-1)*nxC  + nVxTot + nVyTot  	];
		LocVx->PW      = Numbering->map[ ix-1    + (iy-1)*nxC  + nVxTot + nVyTot	];
	}
	else {
		LocVx->VxC     = ix      + iy*nxVx                 			;
		LocVx->VxN     = ix      + iy*nxVx     + nxVx      			;
		LocVx->VxS     = ix      + iy*nxVx     - nxVx      			;
		LocVx->VxE     = ix      + iy*nxVx     + 1         			;
		LocVx->VxW     = ix      + iy*nxVx     - 1        			;
		LocVx->VyNE    = ix+1    + iy*nxVy     + nVxTot   			;
		LocVx->VyNW    = ix+0    + iy*nxVy     + nVxTot   			;
		LocVx->VySE    = ix+1    + (iy-1)*nxVy + nVxTot   			;
		LocVx->VySW    = ix+0    + (iy-1)*nxVy + nVxTot   			;
		LocVx->NormalE =  ix     + (iy-1)*nxC             			;
		LocVx->NormalW = ix-1    + (iy-1)*nxC              			;
		LocVx->ShearN  = ix      + iy*nxS                  			;
		LocVx->ShearS  = ix      + (iy-1)*nxS              			;
		LocVx->PE      = ix      + (iy-1)*nxC  + nVxTot + nVyTot	;
		LocVx->PW      = ix-1    + (iy-1)*nxC  + nVxTot + nVyTot	;
	}


}


void Numbering_getLocalVy(int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingVy* LocVy, bool useNumMap)
{
	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;
	int nVyTot 	= Grid->nVyTot;
	int nxC 	= Grid->nxC;
	int nxS 	= Grid->nxS;
	if (useNumMap) {
		LocVy->VyC     = Numbering->map[ ix     + iy*nxVy + nVxTot          	];
		LocVy->VyN     = Numbering->map[ ix     + iy*nxVy + nVxTot   + nxVy     ];
		LocVy->VyS     = Numbering->map[ ix     + iy*nxVy + nVxTot   - nxVy     ];
		LocVy->VyE     = Numbering->map[ ix     + iy*nxVy + nVxTot   + 1        ];
		LocVy->VyW     = Numbering->map[ ix     + iy*nxVy + nVxTot   - 1        ];
		LocVy->VxNE    = Numbering->map[ ix     + (iy+1)*nxVx               	];
		LocVy->VxNW    = Numbering->map[ ix     + (iy+1)*nxVx - 1          		];
		LocVy->VxSE    = Numbering->map[ ix     + (iy  )*nxVx               	];
		LocVy->VxSW    = Numbering->map[ ix     + (iy  )*nxVx - 1           	];
		LocVy->NormalN = Numbering->map[ ix-1   + (iy  )*nxC                	];
		LocVy->NormalS = Numbering->map[ ix-1   + (iy-1)*nxC                	];
		LocVy->ShearE  = Numbering->map[ ix     + iy*nxS                    	];
		LocVy->ShearW  = Numbering->map[ ix-1   + iy*nxS                    	];
		LocVy->PN      = Numbering->map[ ix-1   + (iy  )*nxC + nVxTot + nVyTot 	];
		LocVy->PS      = Numbering->map[ ix-1   + (iy-1)*nxC + nVxTot + nVyTot 	];
	}
	else {
		LocVy->VyC     = ix     + iy*nxVy + nVxTot          	;
		LocVy->VyN     = ix     + iy*nxVy + nVxTot   + nxVy     ;
		LocVy->VyS     = ix     + iy*nxVy + nVxTot   - nxVy     ;
		LocVy->VyE     = ix     + iy*nxVy + nVxTot   + 1        ;
		LocVy->VyW     = ix     + iy*nxVy + nVxTot   - 1        ;
		LocVy->VxNE    = ix     + (iy+1)*nxVx               	;
		LocVy->VxNW    = ix     + (iy+1)*nxVx - 1          		;
		LocVy->VxSE    = ix     + (iy  )*nxVx               	;
		LocVy->VxSW    = ix     + (iy  )*nxVx - 1           	;
		LocVy->NormalN = ix-1   + (iy  )*nxC                	;
		LocVy->NormalS = ix-1   + (iy-1)*nxC                	;
		LocVy->ShearE  = ix     + iy*nxS                    	;
		LocVy->ShearW  = ix-1   + iy*nxS                    	;
		LocVy->PN      = ix-1   + (iy  )*nxC + nVxTot + nVyTot 	;
		LocVy->PS      = ix-1   + (iy-1)*nxC + nVxTot + nVyTot 	;
	}

}



void Numbering_getLocalP(int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingP* LocP, bool useNumMap)
{
	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;

	if (useNumMap) {
		LocP->VxE     = Numbering->map[ ix+1 + (iy+1)*nxVx               ];
		LocP->VxW     = Numbering->map[ ix   + (iy+1)*nxVx               ];
		LocP->VyN     = Numbering->map[ ix+1 + (iy+1)*(nxVy)  + nVxTot   ];
		LocP->VyS     = Numbering->map[ ix+1 + iy*(nxVy)      + nVxTot	];
	}
	else {
		LocP->VxE     = ix+1 + (iy+1)*nxVx               ;
		LocP->VxW     = ix   + (iy+1)*nxVx               ;
		LocP->VyN     = ix+1 + (iy+1)*(nxVy)  + nVxTot   ;
		LocP->VyS     = ix+1 + iy*(nxVy)      + nVxTot	;
	}

}











