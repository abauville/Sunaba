/*
 * eqSystem.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"

void EqSystem_allocateI (EqSystem* EqSystem)
{
	EqSystem->I 			= (int*) malloc((EqSystem->nEq+1)  * sizeof(int));
}

void EqSystem_allocateMemory(EqSystem* EqSystem)
{
	EqSystem->J = (int*)     malloc(EqSystem->nnz * sizeof(int));
	EqSystem->V = (compute*) malloc(EqSystem->nnz * sizeof(compute));
	EqSystem->b = (compute*) malloc( EqSystem->nEq * sizeof(compute));
	EqSystem->x = (compute*) malloc( EqSystem->nEq * sizeof(compute));
}

void EqSystem_freeMemory(EqSystem* EqSystem) {
	free(EqSystem->I);
	free(EqSystem->J);
	free(EqSystem->V);
	free(EqSystem->b);
	free(EqSystem->x);
}

void EqSystem_assemble(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering)
{


	//==========================================================================
	//
	//                          BUILD SPARSE TRIPLET
	//
	//==========================================================================
	INIT_TIMER
	TIC
	// Init
	// =======================

	int i;
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->J[i] = -1;
	}
	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->b[i] = 0.0;
		EqSystem->x[i]   = 0.0;
	}

	if (DEBUG) {
		/*
		printf("===== Isparse before =====\n");
		printListi(   EqSystem->I,EqSystem->nEq+1);
		printf("nEq = %i",EqSystem->nEq);
		printf("===== EqSystem->J before filling =====\n");
		printListi(   EqSystem->J,EqSystem->nnz);
		printf("===== Numbering->IX =====\n");
		printListi(Numbering->IX,EqSystem->nEq);

		printf("===== Numbering->IY =====\n");
		printListi(Numbering->IY,EqSystem->nEq);
		 */
	}




	// Fill J, V and EqSystem->b for Free and Dirichlet nodes
	// ===============================================
	if (DEBUG) {
		printf("Start Filling loop\n");
		printf("nEq = %i, EqSystem->nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	}
	int iEq, ix, iy, I;
	int Type = 0;
	int INumMap;

	for (iEq=0; iEq<EqSystem->nEq; iEq++) {

		I = EqSystem->I[iEq];
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];
		if (iEq<EqSystem->VyEq0) {
			INumMap = ix+iy*Grid->nxVx;
		}
		else if (iEq<EqSystem->PEq0) {
			Type = 1;
			INumMap = ix+iy*Grid->nxVy + Grid->nVxTot;
		}
		else {
			Type = 2;
			INumMap = ix+iy*Grid->nxC + Grid->nVxTot + Grid->nVyTot;

		}

		if (BC->isNeu[INumMap]==false) { // If Free equation (i.e. not Neumann equation)
			fill_J_V_local(Type, BC->SetupType, ix, iy, I, iEq, EqSystem->J, EqSystem->V, EqSystem->b, Grid->nxC, Grid->nyC, Grid->dx, Grid->dy, Numbering->map, Physics->etaShear, Physics->eta, BC->listDir, BC->valueDir, BC->nDir);
		}

	}
	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
			EqSystem->J[EqSystem->I[i]] = i;
			EqSystem->V[EqSystem->I[i]] = 0.0;
		}
	}

	if (DEBUG) {
		// List J per row
		/*
		int j;
		printf("===== J per row before Neu =====\n");
		for (i=0; i<EqSystem->nEq; i++) {
			I = EqSystem->I[i];
			printf("row #%*i :",3,i);
			for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
				printf("%*i",4,EqSystem->J[I+j]);
			}
			printf("\n");
		}
		printf("\n");
		 */
	}

	// Fill J, V, b for Neumann nodes
	// ================================
	// Velocity Neumann
	int J=0;
	for (i=0; i<BC->nNeu; i++) {
		I = BC->listNeu[i];
		J = EqSystem->I[I];

		if (BC->listNeu[i]<BC->listNeuNeigh[i]) {
			EqSystem->J[J+0] = BC->listNeu[i];
			EqSystem->J[J+1] = BC->listNeuNeigh[i];

			EqSystem->V[J+0] = BC->coeffNeu[i];
			EqSystem->V[J+1] = BC->coeffNeuNeigh[i];

			EqSystem->b[I  ] = BC->valueNeu[i];
		}
		else {
			if (UPPER_TRI) {
				EqSystem->J[J+0] = BC->listNeu[i];
				EqSystem->V[J+0] = BC->coeffNeu[i];
				EqSystem->b[I  ] = BC->valueNeu[i];
			}
			else {
				EqSystem->J[J+0] = BC->listNeuNeigh[i];
				EqSystem->J[J+1] = BC->listNeu[i];
				EqSystem->V[J+0] = BC->coeffNeuNeigh[i];
				EqSystem->V[J+1] = BC->coeffNeu[i];
				EqSystem->b[I  ] = BC->valueNeu[i];
			}
		}
	}

	TOC
	printf("Building the system of equation took: %.2f s\n", toc);

}

void EqSystem_check(EqSystem* EqSystem)
{
	int i, j, I;
	// Check

	printf(" ===== Isparse =====\n");
	for (i=0;i<EqSystem->nRow+1;i++) {
		printf("%i  ", EqSystem->I[i]);
	}
	printf(" \n");


	printf(" ===== Jsparse =====\n");
	for (i=0;i<EqSystem->nnz;i++) {
		printf("%i  ", EqSystem->J[i]);
	}
	printf(" \n");

	printf(" ===== Vsparse =====\n");
	for (i=0;i<EqSystem->nnz;i++) {
		printf("%.3f  ", EqSystem->V[i]);
	}
	printf(" \n");




	printf("===== SPY =====\n");
	printf("   ");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("%*i ",2,i);
	}
	printf("\n");
	int padding;
	for (i=0; i<EqSystem->nRow; i++) {
		I = EqSystem->I[i];
		printf("%*i ",2,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			if (j==0){
				padding = EqSystem->J[I+j]+1;
			}
			else{
				padding = EqSystem->J[I+j]-EqSystem->J[I+j-1];
			}
			printf("%*s",3*padding,"x ");
		}
		printf("\n");
	}


	// List J per row
	printf("===== J per row =====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*i",4,EqSystem->J[I+j]);
		}
		printf("\n");
	}
	printf("\n");


	printf("===== V per row=====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*.2f",6,EqSystem->V[I+j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("===== RHS =====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("RHS[%i] = %.2f\n", i, EqSystem->b[i]);
	}
	printf("\n");
}





void fill_J_V_local(int Type, int BCtype, int ix, int iy,int I, int iEq, int *Jsparse, compute* Vsparse, compute* RHS, int nxC, int nyC, compute dx, compute dy, int *NumMap, compute* EtaShear, compute* EtaNormal, int* BC_List_Dir, compute* BC_Value_Dir, int nDir)
{
	// Computes the column indices and corresponding value for a local equation defined by ix, iy
	// Type=       0: Vx      1: Vy      2: p
	// (ix,iy) indices for the Vx grid (0,0) corresponds to the first Vx equation or first Vy or first P
	// locJ contains the column indices for the current row
	// locV contains the corresponding values
	// NumMap is the numbering Map (continuous numbering without Dirichlet and -1 when Dirichlet)
	// EtaShear and EtaNormal are lists of viscosities

	// Init variables
	// ===============================
	int nLoc, nxVx, nyVx, nxVy, nyVy, nVxTot, nVyTot, nxN, nxS;
	int iDir, i, J;
	int shift;

	int NormalE, NormalW, NormalS, NormalN;
	int ShearE, ShearW, ShearS, ShearN;

	compute EtaN, EtaS, EtaW, EtaE;


	// Define size variables
	// ===============================


	nxVx = nxC+1; // number of Vx nodes in x
	nyVx = nyC+2; // number of Vx nodes in y
	nxVy = nxC+2; // number of Vy nodes in x
	nyVy = nyC+1; // number of Vy nodes in y

	nVxTot = nxVx*nyVx;
	nVyTot = nxVy*nyVy;

	nxN = nxC;

	nxS = nxC+1;

	int Jloc[11];
	compute Vloc[11];
	if (Type==0)
	{
		// =========================================================================
		// =========================================================================
		//
		//                                   Vx
		//
		// =========================================================================
		// =========================================================================


		nLoc = 11; // Maximum number of non zeros for Stokes on the staggered grid
		if (UPPER_TRI) {
			shift = 2;
		}
		else {
			shift = 0;
		}


		// =====================================================================
		//                               locJ
		// =====================================================================

		// Fill Jloc: list of all J indices (including Dirichlet)
		// ================================================================
		Jloc[ 0] =   ix      + iy*nxVx     - nxVx          ; // VxS
		Jloc[ 1] =   ix      + iy*nxVx     - 1             ; // VxW
		Jloc[ 2] =   ix      + iy*nxVx                     ; // VxC
		Jloc[ 3] =   ix      + iy*nxVx     + 1             ; // VxE
		Jloc[ 4] =   ix      + iy*nxVx     + nxVx          ; // VxN
		Jloc[ 5] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySw
		Jloc[ 6] =   ix+1    + (iy-1)*nxVy + nVxTot        ; // VySE
		Jloc[ 7] =   ix+0    + iy*nxVy     + nVxTot        ; // VyNW
		Jloc[ 8] =   ix+1    + iy*nxVy     + nVxTot        ; // VyNE
		Jloc[ 9] =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot ; // PW
		Jloc[10] =   ix      + (iy-1)*nxN  + nVxTot+nVyTot ; // PE


		// =====================================================================
		//                               locV
		// =====================================================================
		// Get Viscosities
		// ================
		NormalE = ix      + (iy-1)*nxN;
		NormalW = ix-1    + (iy-1)*nxN;
		ShearN  = ix      + iy*nxS    ;
		ShearS  = ix      + (iy-1)*nxS;

		EtaN    = EtaShear[ShearN];
		EtaS    = EtaShear[ShearS];
		EtaE    = EtaNormal[NormalE];
		EtaW    = EtaNormal[NormalW];

		// Fill Vloc: list of coefficients
		// ================================
		Vloc[ 0] =  EtaS/dy/dy;
		Vloc[ 1] =  2.0 * EtaW/dx/dx;
		Vloc[ 2] = -2.0 * EtaE/dx/dx   -2.0 * EtaW/dx/dx   -1.0 * EtaN/dy/dy   -1.0 * EtaS/dy/dy;
		Vloc[ 3] =  2.0 * EtaE/dx/dx;
		Vloc[ 4] =  EtaN/dy/dy;
		Vloc[ 5] =  EtaS/dx/dy;
		Vloc[ 6] = -EtaS/dx/dy;
		Vloc[ 7] = -EtaN/dx/dy;
		Vloc[ 8] =  EtaN/dx/dy;
		Vloc[ 9] =  1.0/dx;
		Vloc[10] = -1.0/dx;


		// Adjustment for the case where boundary conditions are periodic horizontally
		if (BCtype==1) {
			if (ix==0) {
				if (UPPER_TRI) {
					shift = 1;
				}
				Jloc[ 0] =   ix      + iy*nxVx     - nxVx          ; // VxS
				Jloc[ 1] =   ix      + iy*nxVx                     ; // VxC
				Jloc[ 2] =   ix      + iy*nxVx     + 1             ; // VxE
				Jloc[ 3] =   ix      + iy*nxVx     - 1      +nxVx-1; // VxW **
				Jloc[ 4] =   ix      + iy*nxVx     + nxVx          ; // VxN
				Jloc[ 5] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySw
				Jloc[ 6] =   ix+1    + (iy-1)*nxVy + nVxTot        ; // VySE
				Jloc[ 7] =   ix+0    + iy*nxVy     + nVxTot        ; // VyNW
				Jloc[ 8] =   ix+1    + iy*nxVy     + nVxTot        ; // VyNE
				Jloc[ 9] =   ix      + (iy-1)*nxN  + nVxTot+nVyTot ; // PE
				Jloc[10] =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot  + nxN; // PW **

				// =====================================================================
				//                               locV
				// =====================================================================
				// Get Viscosities
				// ================
				NormalW = ix-1    + (iy-1)*nxN + nxN;
				EtaW    = EtaNormal[NormalW];

				// Fill Vloc: list of coefficients
				// ================================
				Vloc[ 0] =  EtaS/dy/dy;
				Vloc[ 1] = -2.0 * EtaE/dx/dx   -2.0 * EtaW/dx/dx   -1.0 * EtaN/dy/dy   -1.0 * EtaS/dy/dy;
				Vloc[ 2] =  2.0 * EtaE/dx/dx;
				Vloc[ 3] =  2.0 * EtaW/dx/dx;
				Vloc[ 4] =  EtaN/dy/dy;
				Vloc[ 5] =  EtaS/dx/dy;
				Vloc[ 6] = -EtaS/dx/dy;
				Vloc[ 7] = -EtaN/dx/dy;
				Vloc[ 8] =  EtaN/dx/dy;
				Vloc[ 9] = -1.0/dx;
				Vloc[10] =  1.0/dx;



			}


			if (ix==nxVx-2) {
				if (UPPER_TRI) {
					shift = 3;
				}
				Jloc[ 0] =   ix      + iy*nxVx     - nxVx          ; // VxS
				Jloc[ 1] =   ix      + iy*nxVx     + 1             ; // VxE
				Jloc[ 2] =   ix      + iy*nxVx     - 1             ; // VxW
				Jloc[ 3] =   ix      + iy*nxVx                     ; // VxC

				Jloc[ 4] =   ix      + iy*nxVx     + nxVx          ; // VxN
				Jloc[ 5] =   ix+1    + (iy-1)*nxVy + nVxTot        ; // VySE **
				Jloc[ 6] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySw
				Jloc[ 7] =   ix+1    + iy*nxVy     + nVxTot        ; // VyNE **
				Jloc[ 8] =   ix+0    + iy*nxVy     + nVxTot        ; // VyNW

				Jloc[ 9] =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot ; // PW
				Jloc[10] =   ix      + (iy-1)*nxN  + nVxTot+nVyTot ; // PE



				// Fill Vloc: list of coefficients
				// ================================
				Vloc[ 0] =  EtaS/dy/dy;
				Vloc[ 1] =  2.0 * EtaE/dx/dx;
				Vloc[ 2] =  2.0 * EtaW/dx/dx;
				Vloc[ 3] = -2.0 * EtaE/dx/dx   -2.0 * EtaW/dx/dx   -1.0 * EtaN/dy/dy   -1.0 * EtaS/dy/dy;

				Vloc[ 4] =  EtaN/dy/dy;
				Vloc[ 5] = -EtaS/dx/dy;
				Vloc[ 6] =  EtaS/dx/dy;
				Vloc[ 7] =  EtaN/dx/dy;
				Vloc[ 8] = -EtaN/dx/dy;

				Vloc[ 9] =  1.0/dx;
				Vloc[10] = -1.0/dx;


			}






		}


	}
	else if (Type==1)
	{
		// =========================================================================
		// =========================================================================
		//
		//                                   Vy
		//
		// =========================================================================
		// =========================================================================

		nLoc = 11; // Maximum number of non zeros for Stokes on the staggered grid
		if (UPPER_TRI) {
			shift = 6;
		}
		else {
			shift = 0;
		}
		// =====================================================================
		//                               locJ
		// =====================================================================

		// Fill Jloc: list of all J indices (including Dirichlet)
		// ================================================================
		Jloc[ 0] =   ix      + (iy  )*nxVx - 1                   ; // VxSW
		Jloc[ 1] =   ix      + (iy  )*nxVx                       ; // VxSE
		Jloc[ 2] =   ix      + (iy+1)*nxVx - 1                   ; // VxNW
		Jloc[ 3] =   ix      + (iy+1)*nxVx                       ; // VxNE
		Jloc[ 4] =   ix      + iy*nxVy     + nVxTot    - nxVy    ; // VyS
		Jloc[ 5] =   ix      + iy*nxVy     + nVxTot   - 1        ; // VyW
		Jloc[ 6] =   ix      + iy*nxVy     + nVxTot              ; // VyC
		Jloc[ 7] =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE
		Jloc[ 8] =   ix      + iy*nxVy     + nVxTot    + nxVy    ; // VyN
		Jloc[ 9] =   ix-1    + (iy-1)*nxN + nVxTot+nVyTot        ; // PS
		Jloc[10] =   ix-1    + (iy  )*nxN + nVxTot+nVyTot        ; // PN


		// =====================================================================
		//                               locV
		// =====================================================================
		// Get Viscosities
		// ================
		NormalN = ix-1    + (iy  )*nxN;
		NormalS = ix-1    + (iy-1)*nxN;
		ShearE  = ix      + iy*nxS    ;
		ShearW  = ix-1    + iy*nxS    ;


		EtaN    = EtaNormal[NormalN];
		EtaS    = EtaNormal[NormalS];
		EtaE    = EtaShear[ShearE];
		EtaW    = EtaShear[ShearW];


		// Fill Vloc: list of coefficients
		// ================================
		Vloc[ 0] =  EtaW/dy/dx;
		Vloc[ 1] = -EtaE/dy/dx;
		Vloc[ 2] = -EtaW/dy/dx;
		Vloc[ 3] =  EtaE/dy/dx;
		Vloc[ 4] =  2.0 * EtaS/dy/dy;
		Vloc[ 5] =  EtaW/dx/dx;
		Vloc[ 6] = -2.0 * EtaN/dy/dy   -2.0 * EtaS/dy/dy   -1.0 * EtaE/dx/dx   -1.0 * EtaW/dx/dx;
		Vloc[ 7] =  EtaE/dx/dx;
		Vloc[ 8] =  2.0 * EtaN/dy/dy;
		Vloc[ 9] =  1.0/dy;
		Vloc[10] = -1.0/dy;




		if (BCtype==1) {
			if (ix==0) {
				if (UPPER_TRI) {
					shift = 5;
				}

				Jloc[ 0] =   ix      + (iy  )*nxVx                       ; // VxSE
				Jloc[ 1] =   ix      + (iy  )*nxVx - 1      + nxVx-1     ; // VxSW **
				Jloc[ 2] =   ix      + (iy+1)*nxVx                       ; // VxNE
				Jloc[ 3] =   ix      + (iy+1)*nxVx - 1      + nxVx-1     ; // VxNW **
				Jloc[ 4] =   ix      + iy*nxVy     + nVxTot    - nxVy    ; // VyS
				Jloc[ 5] =   ix      + iy*nxVy     + nVxTot              ; // VyC
				Jloc[ 6] =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE
				Jloc[ 7] =   ix      + iy*nxVy     + nVxTot   - 1      +nxVy-2  ; // VyW **
				Jloc[ 8] =   ix      + iy*nxVy     + nVxTot    + nxVy    ; // VyN
				Jloc[ 9] =   ix-1    + (iy-1)*nxN + nVxTot+nVyTot    + nxN    ; // PS **
				Jloc[10] =   ix-1    + (iy  )*nxN + nVxTot+nVyTot    + nxN    ; // PN **


				// Get Viscosities
				// ================
				NormalN = ix-1    + (iy  )*nxN + nxN; //**
				NormalS = ix-1    + (iy-1)*nxN + nxN; //**
				ShearW  = ix-1    + iy*nxS   + nxS-1 ;


				EtaN    = EtaNormal[NormalN];
				EtaS    = EtaNormal[NormalS];
				EtaW    = EtaShear[ShearW];


				// Fill Vloc: list of coefficients
				// ================================

				Vloc[ 0] = -EtaE/dy/dx;
				Vloc[ 1] =  EtaW/dy/dx;
				Vloc[ 2] =  EtaE/dy/dx;
				Vloc[ 3] = -EtaW/dy/dx;
				Vloc[ 4] =  2.0 * EtaS/dy/dy;
				Vloc[ 5] = -2.0 * EtaN/dy/dy   -2.0 * EtaS/dy/dy   -1.0 * EtaE/dx/dx   -1.0 * EtaW/dx/dx;
				Vloc[ 6] =  EtaE/dx/dx;
				Vloc[ 7] =  EtaW/dx/dx;
				Vloc[ 8] =  2.0 * EtaN/dy/dy;
				Vloc[ 9] =  1.0/dy;
				Vloc[10] = -1.0/dy;




			}
			if (ix==nxVy-3) {
				if (UPPER_TRI) {
					shift = 7;
				}

				// Fill Jloc: list of all J indices (including Dirichlet)
				// ================================================================

				Jloc[ 5] =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE
				Jloc[ 6] =   ix      + iy*nxVy     + nVxTot   - 1        ; // VyW
				Jloc[ 7] =   ix      + iy*nxVy     + nVxTot              ; // VyC


				// =====================================================================
				//                               locV
				// =====================================================================

				Vloc[ 5] =  EtaE/dx/dx;
				Vloc[ 6] =  EtaW/dx/dx;
				Vloc[ 7] = -2.0 * EtaN/dy/dy   -2.0 * EtaS/dy/dy   -1.0 * EtaE/dx/dx   -1.0 * EtaW/dx/dx;





			}
		}




	}
	else if (Type==2)
	{
		// =========================================================================
		// =========================================================================
		//
		//                                   P
		//
		// =========================================================================
		// =========================================================================
		nLoc = 4; // Maximum number of non zeros for Stokes on the staggered grid
		if (UPPER_TRI) {
			shift = 4;
		}
		else {
			shift = 0;
		}
		// =====================================================================
		//                               locJ
		// =====================================================================

		// Fill Jloc: list of all J indices (including Dirichlet)
		// ================================================================
		Jloc[0]     =   ix   + (iy+1)*nxVx              ; // VxW
		Jloc[1]     =   ix+1 + (iy+1)*nxVx              ; // VxE
		Jloc[2]     =   ix+1 + iy*(nxVy)      + nVxTot  ; // VyS
		Jloc[3]     =   ix+1 + (iy+1)*(nxVy)  + nVxTot  ; // VyN



		// =====================================================================
		//                               locV
		// =====================================================================
		// Fill Vloc: list of coefficients
		// ================================
		Vloc[0] = -1.0/dx;
		Vloc[1] =  1.0/dx;
		Vloc[2] = -1.0/dy;
		Vloc[3] =  1.0/dy;


		if (BCtype==1) {
			if (ix==nxN-1) {
				if (UPPER_TRI) {
					shift = 4;
				}


				Jloc[0]     =   ix+1 + (iy+1)*nxVx              ; // VxE
				Jloc[1]     =   ix   + (iy+1)*nxVx              ; // VxW
				Jloc[2]     =   ix+1 + iy*(nxVy)      + nVxTot  ; // VyS
				Jloc[3]     =   ix+1 + (iy+1)*(nxVy)  + nVxTot  ; // VyN

				Vloc[0] =  1.0/dx;
				Vloc[1] = -1.0/dx;
				Vloc[2] = -1.0/dy;
				Vloc[3] =  1.0/dy;

			}
		}




	}









	// =========================================================================
	// =========================================================================
	//
	//                        Fill Jsparse, Vsparse
	//
	// =========================================================================
	// =========================================================================

	// Jloc gets the numbering with Dirichlet;
	// When values are inputted in Jsparse they are transformed to NumMap[Jloc[i]]
	// i.e. numbering without dirichlet

	J = 0;
	for (i=0; i<nLoc; i++) {
		if (NumMap[Jloc[i]] !=-1) { // if not Dirichlet
			if (i>=shift) {
				Jsparse[I+J] = NumMap[Jloc[i]];
				Vsparse[I+J] = Vloc[i];
				J++;
			}
		}
		else { // if it is a Dirichlet

			// Loop through the Dirichlet list to find the index inside the list of the equation corresponding to
			// Jloc[i]. The index is then used to retrieve the boundary value in BC_Value_Dir
			iDir = 0;
			while (BC_List_Dir[iDir] !=Jloc[i]) {

				iDir++;
				if (iDir==nDir) {
					fprintf(stderr,"Couldn't find the Dirichlet equations");
					exit(0);
				}
			}


			RHS[iEq] += -Vloc[i] * BC_Value_Dir[iDir];


		}
	}



}

void EqSystem_solve(EqSystem* EqSystem) {
	//int i;
	INIT_TIMER
	TIC
	int err;

	if (UPPER_TRI){
		err = pardisoSolveSymmetric(EqSystem->I ,EqSystem->J, EqSystem->V ,EqSystem->x , EqSystem->b, EqSystem->nEq);
	}
	else {
		err = pardisoSolveAssymmetric(EqSystem->I ,EqSystem->J, EqSystem->V ,EqSystem->x , EqSystem->b, EqSystem->nEq);
	}

	TOC
	if (err==1) {
		printf("The solve failed\n");
	}
	else {
		printf("YEAY!! \n");
	}

	if (DEBUG) {
		/*
		printf("===== SOLUTION =====\n");
		for (i=0; i<EqSystem->nEq; i++) {
			printf("x[%i] = %.1f\n", i, EqSystem->x[i]);
		}
		 */
	}


	printf("Direct solve: %.2f s\n", toc);


	printf("\n");


}





int pardisoSolveAssymmetric(int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n)
{
	/* Matrix data. */
	//int    n = 8;
	/*
    int    ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
    int    ja[20] = { 0,    2,       5, 6,
        1, 2,    4,
        2,             7,
        3,       6,
        1,
        2,       5,    7,
        1,             6,
        2,          6, 7 };
    double  a[20] = { 7.0,      1.0,           2.0, 7.0,
        -4.0, 8.0,      2.0,
        1.0,                     5.0,
        7.0,           9.0,
        -4.0,
        7.0,           3.0,      8.0,
        1.0,                    11.0,
        -3.0,                2.0, 5.0 };
	 */

	int      nnz = ia[n];
	int      mtype = 11;        /* Real unsymmetric matrix */

	/* RHS and solution vectors. */
	//double   b[8], x[8], diag[8];

	double *diag = malloc(n * sizeof(double));

	int      nrhs = 1;          /* Number of right hand sides. */

	/* Internal solver memory pointer pt,                  */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
	/* or void *pt[64] should be OK on both architectures  */
	void    *pt[64];

	/* Pardiso control parameters. */
	int      iparm[64];
	double   dparm[64];
	int      solver;
	int      maxfct, mnum, phase, error, msglvl;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;
	int      i, k;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters and initialize the solvers      */
	/*     internal adress pointers. This is only necessary for the FIRST   */
	/*     call of the PARDISO solver.                                      */
	/* ---------------------------------------------------------------------*/

	error = 0;
	solver = 0; /* use sparse direct solver */
	pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

	if (error != 0)
	{
		if (error == -10 )
			printf("No license file found \n");
		if (error == -11 )
			printf("License is expired \n");
		if (error == -12 )
			printf("Wrong username or hostname \n");
		return 1;
	}
	else
		printf("[PARDISO]: License check was successful ... \n");


	/* Numbers of processors, value of OMP_NUM_THREADS */
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	iparm[2]  = num_procs;

	iparm[10] = 0; /* no scaling  */
	iparm[12] = 0; /* no matching */

	maxfct = 1;         /* Maximum number of numerical factorizations.  */
	mnum   = 1;         /* Which factorization to use. */

	msglvl = 1;         /* Print statistical information  */
	error  = 0;         /* Initialize error flag */


	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}


	/* -------------------------------------------------------------------- */
	/*  .. pardiso_chk_matrix(...)                                          */
	/*     Checks the consistency of the given matrix.                      */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */

	pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
	if (error != 0) {
		printf("\nERROR in consistency of matrix: %d", error);
		exit(1);
	}

	/* -------------------------------------------------------------------- */
	/* ..  pardiso_chkvec(...)                                              */
	/*     Checks the given vectors for infinite and NaN values             */
	/*     Input parameters (see PARDISO user manual for a description):    */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */

	pardiso_chkvec (&n, &nrhs, b, &error);
	if (error != 0) {
		printf("\nERROR  in right hand side: %d", error);
		exit(1);
	}

	/* -------------------------------------------------------------------- */
	/* .. pardiso_printstats(...)                                           */
	/*    prints information on the matrix to STDOUT.                       */
	/*    Use this functionality only for debugging purposes                */
	/* -------------------------------------------------------------------- */

	pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
	if (error != 0) {
		printf("\nERROR right hand side: %d", error);
		exit(1);
	}

	/* -------------------------------------------------------------------- */
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */
	phase = 11;

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */
	phase = 22;

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ...\n ");

	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	phase = 33;

	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	printf("\nSolve completed ... ");
	//  printf("\nThe solution of the system is: ");
	//  for (i = 0; i < n; i++) {
	//      printf("\n x [%d] = % f", i, x[i] );
	//  }
	//  printf ("\n");

	/* -------------------------------------------------------------------- */
	/* ..  Back substitution with tranposed matrix A^t x=b                  */
	/* -------------------------------------------------------------------- */

	phase = 33;

	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
	iparm[11] = 1;       /* Solving with transpose matrix. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (i = 0; i < n; i++) {
		printf("\n x [%d] = % f", i, x[i] );
	}
	printf ("\n");



	/* -------------------------------------------------------------------- */
	/* ... compute diagonal elements of the inverse.                        */
	/* -------------------------------------------------------------------- */

	phase = 33;
	iparm[11] = 0;       /* Solving with nontranspose matrix. */
	/* solve for n right hand sides */
	for (k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++) {
			b[i] = 0;
		}
		/* Set k-th right hand side to one. */
		b[k] = 1;

		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
				&n, a, ia, ja, &idum, &nrhs,
				iparm, &msglvl, b, x, &error,  dparm);

		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			exit(3);
		}

		/* save diagonal element */
		diag[k] = x[k];
	}

	/* -------------------------------------------------------------------- */
	/* ... Inverse factorization.                                           */
	/* -------------------------------------------------------------------- */

	if (solver == 0)
	{
		printf("\nCompute Diagonal Elements of the inverse of A ... \n");
		phase = -22;
		iparm[35]  = 0; /*  overwrite internal factor L */

		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
				&n, a, ia, ja, &idum, &nrhs,
				iparm, &msglvl, b, x, &error,  dparm);

		/* print diagonal elements */
		for (k = 0; k < n; k++)
		{
			int j = ia[k]-1;
			printf ("Diagonal element of A^{-1} = %32.24e =  %32.24e \n", a[j], diag[k]);
		}
	}


	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}

	/* -------------------------------------------------------------------- */
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */
	phase = -1;                 /* Release internal memory. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error,  dparm);

	return 0;

}






int pardisoSolveSymmetric(int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n)
{
	INIT_TIMER

	// Matrix data.
	/*
    int    n = 8;
    int    ia[ 9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
    int    ja[18] = { 0,    2,       5, 6,
        1, 2,    4,
        2,             7,
        3,       6,
        4, 5, 6,
        5,    7,
        6,
        7 };
    double  a[18] = { 7.0,      1.0,           2.0, 7.0,
        -4.0, 8.0,           2.0,
        1.0,                     5.0,
        7.0,           9.0,
        5.0, 1.0, 5.0,
        0.0,      5.0,
        11.0,
        5.0 };
	 */

	int i;
	//int  j, I;



	for (i=0; i<n; i++) {
		x[i] = 0;
	}



	int      nnz = ia[n];
	int      mtype = -2;        /* Real symmetric matrix */

	/* RHS and solution vectors. */
	//double   b[8], x[8];
	int      nrhs = 1;          /* Number of right hand sides. */

	/* Internal solver memory pointer pt,                  */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
	/* or void *pt[64] should be OK on both architectures  */
	void    *pt[64];

	/* Pardiso control parameters. */
	int      iparm[64];
	double   dparm[64];
	int      maxfct, mnum, phase, error, msglvl, solver;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */



	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters.                                */
	/* -------------------------------------------------------------------- */
	if (TIMER) {
		TIC
	}

	error = 1;
	solver=0;// use sparse direct solver
	msglvl = 0;         // Print statistical information
	pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

	if (error != 0)
	{
		if (error == -10 )
			printf("No license file found \n");
		if (error == -11 )
			printf("License is expired \n");
		if (error == -12 )
			printf("Wrong username or hostname \n");
		return 1;
	}
	else {
		if (DEBUG) {
			printf("[PARDISO]: License check was successful ... \n");
		}
	}



	// Numbers of processors, value of OMP_NUM_THREADS
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	iparm[2]  = num_procs;

	maxfct = 1;		// Maximum number of numerical factorizations.
	mnum   = 1;         // Which factorization to use.

	msglvl = 0;         // Print statistical information
	error  = 0;         // Initialize error flag

	if (TIMER) {
		TOC
		printf("Pardiso setup: %.3f s\n", toc);
	}

	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */
	if (TIMER) {
		TIC
	}

	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	if (TIMER) {
		TOC
		printf("0 to 1 based indexing: %.3f s\n", toc);
	}

	/*
    // Set right hand side to i.
    for (i = 0; i < n; i++) {
        b[i] = i;
    }
	 */
	/* -------------------------------------------------------------------- */
	/*  .. pardiso_chk_matrix(...)                                          */
	/*     Checks the consistency of the given matrix.                      */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */
	if  (DEBUG) {
		printf("--  chkmatrix\n");
		pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
		if (error != 0) {
			printf("\nERROR in consistency of matrix: %d", error);
			exit(1);
		}
	}
	/* -------------------------------------------------------------------- */
	/* ..  pardiso_chkvec(...)                                              */
	/*     Checks the given vectors for infinite and NaN values             */
	/*     Input parameters (see PARDISO user manual for a description):    */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */
	if  (DEBUG) {
		printf("--  chkvec\n");
		pardiso_chkvec (&n, &nrhs, b, &error);
		if (error != 0) {
			printf("\nERROR  in right hand side: %d", error);
			exit(1);
		}
	}

	/* -------------------------------------------------------------------- */
	/* .. pardiso_printstats(...)                                           */
	/*    prints information on the matrix to STDOUT.                       */
	/*    Use this functionality only for debugging purposes                */
	/* -------------------------------------------------------------------- */
	if  (DEBUG) {
		printf("--  printstats\n");
		pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
		if (error != 0) {
			printf("\nERROR right hand side: %d", error);
			exit(1);
		}
	}
	/* -------------------------------------------------------------------- */
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */



	if (TIMER) {
		TIC
	}
	phase = 11;

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	if (TIMER) {
		TOC
		printf("\nPhase 11: %.3f s\n", toc);
	}



	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */



	if (TIMER) {
		TIC
	}
	phase = 22;
	iparm[32] = 1; /* compute determinant */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ...\n ");

	if (TIMER) {
		TOC
		printf("Phase 33: %.3f s\n", toc);
	}



	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */




	if (TIMER) {
		TIC
	}
	phase = 33;


	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	printf("\nSolve completed ... ");
	if  (DEBUG) {
		/*
		printf("\nThe solution of the system is: ");

		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		 */
	}
	printf ("\n\n");

	if (TIMER) {
		TOC
		printf("Phase 33: %.3f s\n", toc);
	}



	/* -------------------------------------------------------------------- */
	/* ... Inverse factorization.                                           */
	/* -------------------------------------------------------------------- */

	if (solver == 0)
	{
		if (TIMER) {
			TIC
		}
		printf("\nCompute Diagonal Elements of the inverse of A ... \n");
		phase = -22;
		iparm[35]  = 1; //  no not overwrite internal factor L

		pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
				iparm, &msglvl, b, x, &error,  dparm);
		if (TIMER) {
			TOC
			printf("Phase -22: %.3f s\n", toc);
		}


		/*
        // print diagonal elements
        for (k = 0; k < n; k++)
        {
            int j = ia[k]-1;
            printf ("Diagonal element of A^{-1} = %d %d %32.24e\n", k, ja[j]-1, a[j]);
        }
		 */
	}


	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}

	/* -------------------------------------------------------------------- */
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */
	phase = -1;                 /* Release internal memory. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error,  dparm);

	return 0;









}


