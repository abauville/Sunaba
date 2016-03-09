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
	if (EqSystem->penaltyMethod) {
		EqSystem->nEq += BC->nPDir - Grid->nCTot;
	}

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
		if (!EqSystem->penaltyMethod) {
			EqSystem->nRow = EqSystem->nEq + BC->nPDir - Grid->nCTot;
		}
		else {
			EqSystem->nRow = EqSystem->nEq;
		}

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
					Numbering_getLocalVx(ix, iy, Numbering, Grid, BC, &LocVx, true, EqSystem->penaltyMethod);


					EqSystem->nnz += LocVx.VxS + LocVx.VxW + LocVx.VxC + LocVx.VxE + LocVx.VxN
							+ LocVx.VySW + LocVx.VySE + LocVx.VyNW + LocVx.VyNE
							+ LocVx.PW + LocVx.PE; // VxS etc... are 1 if the equation they refer to is free, 0 otherwise

				}

				if ( !(BC->SetupType==1 && ix==Grid->nxVx-1) ) { // To jump the rightmost nodes for periodic bc
					EqSystem->I[InoDir+1] = EqSystem->nnz;
					Numbering->IX[InoDir] = ix;
					Numbering->IY[InoDir] = iy;
					InoDir++;
				}

			}

			I++;
		}
	}









	// Vy
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
					Numbering_getLocalVy(ix, iy, Numbering, Grid, BC, &LocVy, true, EqSystem->penaltyMethod);

					EqSystem->nnz +=  LocVy.VxSW + LocVy.VxSE + LocVy.VxNW + LocVy.VxNE
							+ LocVy.VyS + LocVy.VyW + LocVy.VyC + LocVy.VyE + LocVy.VyN
							+ LocVy.PS + LocVy.PN; // VxSW etc... are 1 if the equation they refer to is free, 0 otherwise


				}
				if ( !(BC->SetupType==1 && ix>=Grid->nxVy-2) ) { // To jump the rightmost nodes for periodic bc
					EqSystem->I[InoDir+1] = EqSystem->nnz;
					Numbering->IX[InoDir] = ix;
					Numbering->IY[InoDir] = iy;
					InoDir++;
				}
			}
			I++;
		}
	}



	EqSystem->PEq0 = InoDir;

	if (!EqSystem->penaltyMethod) {


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

						Numbering_getLocalP(ix, iy, Numbering, Grid, BC, &LocP, true, EqSystem->penaltyMethod);
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
	}

	printf("InoDir = %i nEq = %i, nRow = %i\n",InoDir, EqSystem->nEq, EqSystem->nRow);

	// 3. Numbering->map = cumsum(Numbering->map);

	// Start numbering at 0, but jump the first dirichlet nodes
	i = 0;
	while (Numbering->map[i]==0){
		i++;
	}
	Numbering->map[i] = 0;


	if (BC->SetupType==1) // Number the Equations on the Right boundary with the number from the left one
	{
		int Ileft;
		i = 0;
		int iPrevious = 0;
		// Vx
		for (iy=0; iy<Grid->nyVx; iy++) {
			Ileft = i;
			for (ix=0; ix<Grid->nxVx-1; ix++) {
				Numbering->map[i] += Numbering->map[iPrevious];
				i++;
				iPrevious = i-1;
			}
			Numbering->map[i] = Numbering->map[Ileft];
			i++;
			//iPrevious = i-2;
		}

		// Vy
		for (iy=0; iy<Grid->nyVy; iy++) {
			Ileft = i;
			for (ix=0; ix<Grid->nxVy-2; ix++) {
				Numbering->map[i] += Numbering->map[iPrevious];
				i++;
				iPrevious = i-1;
			}
			Numbering->map[i] = Numbering->map[Ileft];
			i++;
			Numbering->map[i] = Numbering->map[Ileft+1];
			i++;
			//iPrevious = i-2;
		}

		for (iy=0; iy<Grid->nyC; iy++) {
			for (ix=0; ix<Grid->nxC; ix++) {
				Numbering->map[i] += Numbering->map[iPrevious];
				i++;
				iPrevious = i-1;
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
		printf("===== EqSystem->I =====\n");  printListi  (EqSystem->I       ,EqSystem->nEq+1);
		printf("===== IX =====\n");  printListi  (Numbering->IX       ,EqSystem->nEq);
		printf("===== IY =====\n");  printListi  (Numbering->IY       ,EqSystem->nEq);
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


void Numbering_getLocalVx(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingVx* LocVx, bool useNumMap, bool penaltyMethod)
{


	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;
	int nVyTot 	= Grid->nVyTot;
	int nxC 	= Grid->nxC;
	int nxS 	= Grid->nxS;

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

	if (BC->SetupType==1) {
		if (ix==0) {
			LocVx->VxW += nxVx-1;
			LocVx->PW  += nxC;
		}
	}



	if (useNumMap) {
		LocVx->VxC     = Numbering->map[ LocVx->VxC ];
		LocVx->VxN     = Numbering->map[ LocVx->VxN ];
		LocVx->VxS     = Numbering->map[ LocVx->VxS ];
		LocVx->VxE     = Numbering->map[ LocVx->VxE ];
		LocVx->VxW     = Numbering->map[ LocVx->VxW ];
		LocVx->VyNE    = Numbering->map[ LocVx->VyNE];
		LocVx->VyNW    = Numbering->map[ LocVx->VyNW];
		LocVx->VySE    = Numbering->map[ LocVx->VySE];
		LocVx->VySW    = Numbering->map[ LocVx->VySW];
		LocVx->NormalE = Numbering->map[ LocVx->NormalE];
		LocVx->NormalW = Numbering->map[ LocVx->NormalW];
		LocVx->ShearN  = Numbering->map[ LocVx->ShearN];
		LocVx->ShearS  = Numbering->map[ LocVx->ShearS];
		LocVx->PE      = Numbering->map[ LocVx->PE];
		LocVx->PW      = Numbering->map[ LocVx->PW];


		if (UPPER_TRI) {
			LocVx->VxS = 0;
			LocVx->VxW = 0;

			if (BC->SetupType==1) {
				if (ix==0) {
					LocVx->VxW = 1;
				}
				if (ix==nxVx-2) {
					LocVx->VxE = 0;
				}
			}

		}

		if (penaltyMethod) {
			LocVx->PE = 0;
			LocVx->PW = 0;
		}




	}

	if (BC->SetupType==1) {
		if (ix==nxVx-1) {
			LocVx->VxC     = 0	;
			LocVx->VxN     = 0 			;
			LocVx->VxS     =0		;
			LocVx->VxE     = 0  			;
			LocVx->VxW     = 0   			;
			LocVx->VyNE    = 0 			;
			LocVx->VyNW    = 0			;
			LocVx->VySE    = 0		;
			LocVx->VySW    = 0			;
			LocVx->NormalE =  0         			;
			LocVx->NormalW = 0       			;
			LocVx->ShearN  = 0          			;
			LocVx->ShearS  = 0            			;
			LocVx->PE      = 0	;
			LocVx->PW      = 0	;
		}
	}


}


void Numbering_getLocalVy(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingVy* LocVy, bool useNumMap, bool penaltyMethod)
{
	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;
	int nVyTot 	= Grid->nVyTot;
	int nxC 	= Grid->nxC;
	int nxS 	= Grid->nxS;


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


	if (BC->SetupType==1) {
		if (ix==0) {
			LocVy->VxSW += nxVx-1;
			LocVy->VxNW += nxVx-1;
			LocVy->VyW  += nxVy-2;
			LocVy->PS   += nxC;
			LocVy->PN   += nxC;
		}
	}



	if (useNumMap) {
		LocVy->VyC     = Numbering->map[ LocVy->VyC      	];
		LocVy->VyN     = Numbering->map[ LocVy->VyN  ];
		LocVy->VyS     = Numbering->map[ LocVy->VyS     ];
		LocVy->VyE     = Numbering->map[ LocVy->VyE       ];
		LocVy->VyW     = Numbering->map[ LocVy->VyW     ];
		LocVy->VxNE    = Numbering->map[ LocVy->VxNE       	];
		LocVy->VxNW    = Numbering->map[ LocVy->VxNW          		];
		LocVy->VxSE    = Numbering->map[ LocVy->VxSE     	];
		LocVy->VxSW    = Numbering->map[ LocVy->VxSW        	];
		LocVy->NormalN = Numbering->map[ LocVy->NormalN       	];
		LocVy->NormalS = Numbering->map[ LocVy->NormalS             	];
		LocVy->ShearE  = Numbering->map[ LocVy->ShearE              	];
		LocVy->ShearW  = Numbering->map[ LocVy->ShearW         	];
		LocVy->PN      = Numbering->map[ LocVy->PN 	];
		LocVy->PS      = Numbering->map[ LocVy->PS	];

		if (UPPER_TRI) {
			LocVy->VyS = 0;
			LocVy->VyW = 0;
			LocVy->VxNE = 0;
			LocVy->VxNW = 0;
			LocVy->VxSE = 0;
			LocVy->VxSW = 0;

			if (BC->SetupType==1) {
				if (ix==0) {
					LocVy->VyW = 1;
				}
				else if (ix==nxVy-3) {
					LocVy->VyE = 0;
				}
			}

		}


	}



	if (BC->SetupType==1) {
		if (ix>=nxVy-2) {
			LocVy->VyC     = 0;
			LocVy->VyN     = 0;
			LocVy->VyS     = 0;
			LocVy->VyE     = 0;
			LocVy->VyW     = 0;
			LocVy->VxNE    = 0;
			LocVy->VxNW    = 0;
			LocVy->VxSE    = 0;
			LocVy->VxSW    = 0;
			LocVy->NormalN = 0;
			LocVy->NormalS = 0;
			LocVy->ShearE  = 0;
			LocVy->ShearW  = 0;
			LocVy->PN      = 0;
			LocVy->PS      = 0;
		}
	}

	if (penaltyMethod) {
		LocVy->PN = 0;
		LocVy->PS = 0;
	}

}



void Numbering_getLocalP(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingP* LocP, bool useNumMap, bool penaltyMethod)
{
	int nxVx 	= Grid->nxVx;
	int nxVy 	= Grid->nxVy;
	int nVxTot 	= Grid->nVxTot;


	LocP->VxE     = ix+1 + (iy+1)*nxVx               ;
	LocP->VxW     = ix   + (iy+1)*nxVx               ;
	LocP->VyN     = ix+1 + (iy+1)*(nxVy)  + nVxTot   ;
	LocP->VyS     = ix+1 + iy*(nxVy)      + nVxTot	;

	if (BC->SetupType==1) {
		if (ix==0) {
			LocP->VxW += nxVx-1;
		}
	}

	if (useNumMap) {
		LocP->VxE     = Numbering->map[ ix+1 + (iy+1)*nxVx               ];
		LocP->VxW     = Numbering->map[ ix   + (iy+1)*nxVx               ];
		LocP->VyN     = Numbering->map[ ix+1 + (iy+1)*(nxVy)  + nVxTot   ];
		LocP->VyS     = Numbering->map[ ix+1 + iy*(nxVy)      + nVxTot	];
		if (UPPER_TRI) {
			LocP->VxE     = 0;
			LocP->VxW     = 0;
			LocP->VyN     = 0;
			LocP->VyS     = 0;
		}
	}


}











