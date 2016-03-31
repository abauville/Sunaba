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

void EqSystem_freeMemory(EqSystem* EqSystem, Solver* Solver)
{
	//Free Pardiso

	int error = 0;
	int idum = 0;
	double ddum = 0;
	int phase = -1;                 // Release internal memory.
	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, &ddum, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);

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
	//INIT_TIMER
	//TIC
	// Init
	// =======================

	int i;
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->J[i] = -1;
	}
	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->b[i] = 0.0;
	}





	// Fill J, V and EqSystem->b for Free and Dirichlet nodes
	// ===============================================
	if (DEBUG) {
		printf("Start Filling loop\n");
		printf("nEq = %i, EqSystem->nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	}
	int iEq, ix, iy, I;
	StencilType Stencil;
	//int INumMap;

	for (iEq=0; iEq<EqSystem->nEq; iEq++) {

		I = EqSystem->I[iEq];
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];

		i = 1;
		while (iEq>=Numbering->subEqSystem0[i]) {
			i++;
		}
		Stencil = Numbering->Stencil[i-1];

		fill_J_V_local(Stencil, ix, iy, I, iEq, EqSystem, Grid, Numbering, Physics, BC);


	}
	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
			EqSystem->J[EqSystem->I[i]] = i;
			EqSystem->V[EqSystem->I[i]] = 0.0;
		}
	}





}

void EqSystem_check(EqSystem* EqSystem)
{
	int i, j, I;
	// Check

	printf(" ===== Isparse =====\n");
	for (i=0;i<EqSystem->nEq+1;i++) {
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

	printf(" ===== b=====\n");
	for (i=0;i<EqSystem->nEq;i++) {
		printf("%.3f  ", EqSystem->b[i]);
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
	//for (i=0; i<EqSystem->nEq; i++) {
	for (i=0; i<EqSystem->nRow; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*i",7,EqSystem->J[I+j]);
		}
		printf("\n");
	}
	printf("\n");


	printf("===== V per row=====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*.2f",10,EqSystem->V[I+j]);
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





void fill_J_V_local(StencilType Stencil, int ix, int iy,int I, int iEq, EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC)
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
	int nLoc, nxVx, nyVx, nxVy, nyVy, nVxTot, nVyTot, nxN, nxS, nxEC;
	int i, J;
	int shift;

	int NormalE, NormalW, NormalS, NormalN;
	int ShearE, ShearW, ShearS, ShearN;

	int SetupType = BC->SetupType;

	compute EtaN, EtaS, EtaW, EtaE;
	compute kW, kE, kN, kS;

	compute dx = Grid->dx;
	compute dy = Grid->dy;
	compute dt = Physics->dt;

	int nxC = Grid->nxC;
	int nyC = Grid->nyC;

	// Define size variables
	// ===============================
	nxVx = nxC+1; // number of Vx nodes in x
	nyVx = nyC+2; // number of Vx nodes in y
	nxVy = nxC+2; // number of Vy nodes in x
	nyVy = nyC+1; // number of Vy nodes in y

	nVxTot = nxVx*nyVx;
	nVyTot = nxVy*nyVy;

	nxN = nxC;;
	nxEC = Grid->nxEC;

	nxS = nxC+1;

	int Jloc[11];
	compute Vloc[11];
	compute bloc;
	int order[11] = {0,1,2,3,4,5,6,7,8,9,10};

	int VxW, VxC, VxE, VxS, VxN;
	int VyW, VyC, VyE, VyS, VyN;
	int VxSW, VxSE, VxNW, VxNE;
	int VySW, VySE, VyNW, VyNE;
	int PS, PN, PW, PE;
	int TC, TS, TN, TW, TE;



	if (Stencil==Vx)
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
		//                            Local numbering
		// =====================================================================

		// Fill Jloc: list of all J indices (including Dirichlet)
		// ================================================================
		VxS  =   ix      + iy*nxVx     - nxVx          ; // VxS
		VxW  =   ix      + iy*nxVx     - 1             ; // VxW
		VxC  =   ix      + iy*nxVx                     ; // VxC
		VxE  =   ix      + iy*nxVx     + 1             ; // VxE
		VxN  =   ix      + iy*nxVx     + nxVx          ; // VxN
		VySW =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySW
		VySE =   ix+1    + (iy-1)*nxVy + nVxTot        ; // VySE
		VyNW =   ix+0    + iy*nxVy     + nVxTot        ; // VyNW
		VyNE =   ix+1    + iy*nxVy     + nVxTot        ; // VyNE
		PW   =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot ; // PW
		PE   =   ix      + (iy-1)*nxN  + nVxTot+nVyTot ; // PE


		NormalE = ix  +1    + (iy-1+1)*nxEC  ;
		NormalW = ix-1+1    + (iy-1+1)*nxEC  ;
		ShearN  = ix      + iy*nxS     ;
		ShearS  = ix      + (iy-1)*nxS ;


		// Special case for periodic BC
		if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				if (UPPER_TRI) {
					shift = 1;
				}
				VxW += nxVx-1;
				PW  += nxN;

				NormalW += nxN;

				order[ 0] =  0; // VxS
				order[ 1] =  3; // VxW
				order[ 2] =  1; // VxC
				order[ 3] =  2; // VxE
				order[ 4] =  4; // VxN
				order[ 5] =  5; // VySW
				order[ 6] =  6; // VySE
				order[ 7] =  7; // VyNW
				order[ 8] =  8; // VyNE
				order[ 9] = 10; // PW
				order[10] =  9; // PE
			}
			if (ix==nxVx-2) {
				if (UPPER_TRI) {
					shift = 3;
				}
				order[ 0] =  0; // VxS
				order[ 1] =  2; // VxW
				order[ 2] =  3; // VxC
				order[ 3] =  1; // VxE
				order[ 4] =  4; // VxN
				order[ 5] =  6; // VySW
				order[ 6] =  5; // VySE
				order[ 7] =  8; // VyNW
				order[ 8] =  7; // VyNE
				order[ 9] =  9; // PW
				order[10] = 10; // PE
			}
		}

		// =====================================================================
		//                           Fill Jloc, Vloc, bloc
		// =====================================================================

		Jloc[order[ 0]] = VxS;
		Jloc[order[ 1]] = VxW;
		Jloc[order[ 2]] = VxC;
		Jloc[order[ 3]] = VxE;
		Jloc[order[ 4]] = VxN;
		Jloc[order[ 5]] = VySW;
		Jloc[order[ 6]] = VySE;
		Jloc[order[ 7]] = VyNW;
		Jloc[order[ 8]] = VyNE;
		Jloc[order[ 9]] = PW;
		Jloc[order[10]] = PE;


		EtaN    = Physics->etaShear[ShearN];
		EtaS    = Physics->etaShear[ShearS];
		EtaE    = Physics->eta[NormalE];
		EtaW    = Physics->eta[NormalW];


		// Fill Vloc: list of coefficients
		// ================================
		Vloc[order[ 0]] =  EtaS/dy/dy;
		Vloc[order[ 1]] =  2.0 * EtaW/dx/dx;
		Vloc[order[ 2]] = -2.0 * EtaE/dx/dx   -2.0 * EtaW/dx/dx   -1.0 * EtaN/dy/dy   -1.0 * EtaS/dy/dy;
		Vloc[order[ 3]] =  2.0 * EtaE/dx/dx;
		Vloc[order[ 4]] =  EtaN/dy/dy;
		Vloc[order[ 5]] =  EtaS/dx/dy;
		Vloc[order[ 6]] = -EtaS/dx/dy;
		Vloc[order[ 7]] = -EtaN/dx/dy;
		Vloc[order[ 8]] =  EtaN/dx/dy;
		Vloc[order[ 9]] =  1.0/dx;
		Vloc[order[10]] = -1.0/dx;

		bloc = - Physics->g[0] * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );


	}



	else if (Stencil==Vy)
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
		//                          Local numbering
		// =====================================================================
		VxSW =   ix      + (iy  )*nxVx - 1                   ; // VxSW
		VxSE =   ix      + (iy  )*nxVx                       ; // VxSE
		VxNW =   ix      + (iy+1)*nxVx - 1                   ; // VxNW
		VxNE =   ix      + (iy+1)*nxVx                       ; // VxNE
		VyS  =   ix      + iy*nxVy     + nVxTot    - nxVy    ; // VyS
		VyW  =   ix      + iy*nxVy     + nVxTot   - 1        ; // VyW
		VyC  =   ix      + iy*nxVy     + nVxTot              ; // VyC
		VyE  =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE
		VyN  =   ix      + iy*nxVy     + nVxTot    + nxVy    ; // VyN
		PS   =   ix-1    + (iy-1)*nxN + nVxTot+nVyTot        ; // PS
		PN   =   ix-1    + (iy  )*nxN + nVxTot+nVyTot        ; // PN

		NormalN = ix-1+1    + (iy  +1)*nxEC ;
		NormalS = ix-1+1    + (iy-1+1)*nxEC ;
		ShearE  = ix      + iy*nxS    ;
		ShearW  = ix-1    + iy*nxS    ;


		// Special cases for periodic BC
		if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				if (UPPER_TRI) {
					shift = 5;
				}
				VxSW += nxVx-1		; // VxSW
				VxNW += nxVx-1     	; // VxNW
				VyW  += nxVy-2  	; // VyW
				PS   += nxN    		; // PS
				PN   += nxN    		; // PN

				NormalN += nxN			;
				NormalS += nxN			;
				ShearW  += nxS-1 		;

				order[ 0] =  1; // VxSW
				order[ 1] =  0; // VxSE
				order[ 2] =  3; // VxNW
				order[ 3] =  2; // VxNE
				order[ 4] =  4; // VyS
				order[ 5] =  7; // VyW
				order[ 6] =  5; // VyC
				order[ 7] =  6; // VyE
				order[ 8] =  8; // VyN
				order[ 9] =  9; // PS
				order[10] = 10; // PN
			}
			else if (ix==nxVy-3) {
				if (UPPER_TRI) {
					shift = 7;
				}
				order[ 0] =  0; // VxSW
				order[ 1] =  1; // VxSE
				order[ 2] =  2; // VxNW
				order[ 3] =  3; // VxNE
				order[ 4] =  4; // VyS
				order[ 5] =  6; // VyW
				order[ 6] =  7; // VyC
				order[ 7] =  5; // VyE
				order[ 8] =  8; // VyN
				order[ 9] =  9; // PS
				order[10] = 10; // PN
			}
		}


		// =====================================================================
		//                        fill Jloc, Vloc, bloc
		// =====================================================================
		// Get Viscosities
		// ================
		EtaN    = Physics->eta[NormalN];
		EtaS    = Physics->eta[NormalS];
		EtaE    = Physics->etaShear[ShearE];
		EtaW    = Physics->etaShear[ShearW];


		Jloc[order[ 0]] = VxSW;
		Jloc[order[ 1]] = VxSE;
		Jloc[order[ 2]] = VxNW;
		Jloc[order[ 3]] = VxNE;
		Jloc[order[ 4]] = VyS;
		Jloc[order[ 5]] = VyW;
		Jloc[order[ 6]] = VyC;
		Jloc[order[ 7]] = VyE;
		Jloc[order[ 8]] = VyN;
		Jloc[order[ 9]] = PS;
		Jloc[order[10]] = PN;


		// Fill Vloc: list of coefficients
		// ================================
		Vloc[order[ 0]] =  EtaW/dy/dx; // VxSW
		Vloc[order[ 1]] = -EtaE/dy/dx; // VxSE
		Vloc[order[ 2]] = -EtaW/dy/dx; // VxNW
		Vloc[order[ 3]] =  EtaE/dy/dx; // VxNE
		Vloc[order[ 4]] =  2.0 * EtaS/dy/dy; // VyS
		Vloc[order[ 5]] =  EtaW/dx/dx; 		 //VyW
		Vloc[order[ 6]] = -2.0 * EtaN/dy/dy   -2.0 * EtaS/dy/dy   -1.0 * EtaE/dx/dx   -1.0 * EtaW/dx/dx; // VyC
		Vloc[order[ 7]] =  EtaE/dx/dx; // VyE
		Vloc[order[ 8]] =  2.0 * EtaN/dy/dy; //VyN
		Vloc[order[ 9]] =  1.0/dy; // PS
		Vloc[order[10]] = -1.0/dy; // PN

		bloc = - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );




	}
	else if (Stencil==P)
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
		VxW     =   ix   + (iy+1)*nxVx              ; // VxW
		VxE     =   ix+1 + (iy+1)*nxVx              ; // VxE
		VyS     =   ix+1 + iy*(nxVy)      + nVxTot  ; // VyS
		VyN     =   ix+1 + (iy+1)*(nxVy)  + nVxTot  ; // VyN


		if (SetupType==SimpleShearPeriodic) {

			if (ix==nxN-1) {
				if (UPPER_TRI) {
					shift = 4;
				}
				order[0] = 1;
				order[1] = 0;
				order[2] = 2;
				order[3] = 3;
			}
		}


		Jloc[order[0]] = VxW;
		Jloc[order[1]] = VxE;
		Jloc[order[2]] = VyS;
		Jloc[order[3]] = VyN;


		// =====================================================================
		//                               locV
		// =====================================================================
		// Fill Vloc: list of coefficients
		// ================================
		Vloc[order[0]] = -1.0/dx;
		Vloc[order[1]] =  1.0/dx;
		Vloc[order[2]] = -1.0/dy;
		Vloc[order[3]] =  1.0/dy;



		bloc = 0;


	}

	else if (Stencil==T) {

		nLoc = 5; // Maximum number of non zeros for Stokes on the staggered grid
		if (UPPER_TRI) {
			shift = 2;
		}
		else {
			shift = 0;
		}
		TS =  ix 	+ (iy-1)*(nxC+2);
		TW = (ix-1) +  iy   *(nxC+2);
		TC =  ix 	+  iy   *(nxC+2);
		TE = (ix+1) +  iy   *(nxC+2);
		TN =  ix 	+ (iy+1)*(nxC+2);


		if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				if (UPPER_TRI) {
					shift = 1;
				}
				TW  += (nxC+2)-2  ; // VyW

				order[ 0] =  0; // TS
				order[ 1] =  3; // TW
				order[ 2] =  1; // TC
				order[ 3] =  2; // TE
				order[ 4] =  4; // TN

			}
			else if (ix==nxC-3) {
				if (UPPER_TRI) {
					shift = 3;
				}
				order[ 0] =  0; // TS
				order[ 1] =  2; // TW
				order[ 2] =  3; // TC
				order[ 3] =  1; // TE
				order[ 4] =  4; // TN
			}
		}

		Jloc[order[0]] = TS;
		Jloc[order[1]] = TW;
		Jloc[order[2]] = TC;
		Jloc[order[3]] = TE;
		Jloc[order[4]] = TN;


		kN = (2*Physics->k[TN]*Physics->k[TC])/(Physics->k[TN]+Physics->k[TC]);
		kS = (2*Physics->k[TS]*Physics->k[TC])/(Physics->k[TS]+Physics->k[TC]);
		kW = (2*Physics->k[TW]*Physics->k[TC])/(Physics->k[TW]+Physics->k[TC]);
		kE = (2*Physics->k[TE]*Physics->k[TC])/(Physics->k[TE]+Physics->k[TC]);

		Vloc[order[0]] =  -kS/dy/dy; // TS
		Vloc[order[1]] =  -kW/dx/dx; // TW
		Vloc[order[2]] =  -(-kW/dx/dx -kE/dx/dx -kN/dy/dy -kS/dy/dy) + 1.0/dt; // TC
		Vloc[order[3]] =  -kE/dx/dx; // TE
		Vloc[order[4]] =  -kN/dy/dy; // TN


		bloc = + Physics->T[TC]/dt;
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

	// Fill right hand side with the local right hand side
	EqSystem->b[iEq] += bloc;


	//printListi(Jloc, 11);


	// Deal with Neumann
	int Iloc, IBC, IC;
	if (Stencil==Vx) { // for Vx equation
		IC = 2;
	}
	else if (Stencil==Vy) {
		IC = 6;
	}
	else if (Stencil==T) {
		IC = 2;
	}


		for (i=0; i<nLoc; i++) {
			//printf ("%i ",Jloc[order[i]]);
			Iloc = Numbering->map[Jloc[order[i]]];

			if (Iloc < 0) { // if Boundary node
				IBC = abs(Iloc) - 1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC];
					EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC];
				}
				else  if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Vloc[order[IC]]  += -Vloc[order[i]]; // +1 to VxC
					//EqSystem->b[iEq] += -Vloc[i] * 2*BC->value[IBC];
					EqSystem->b[iEq] += -Vloc[order[i]] * 2*BC->value[IBC];

				}
				else if (BC->type[IBC]==NeumannGhost) { // NeumannGhost
					Vloc[order[IC]] += Vloc[order[i]]; // +1 to VxC
					if (Stencil==Vx) {
						if 		(i==0) { // VxS
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * dy;
						}
						else if (i==4) { // VxN
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * dy;
						}
					}


					else if (Stencil==Vy) {
						if 		(i==5) { // VyW
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * dx;
						}
						else if (i==7) { // VyE
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * dx;
						}
					}

				}
				else {
					printf("error: unknown boundary type\n");
					exit(0);
				}
			}
		}









	//printf("\n");


	J = 0;
	for (i=0; i<nLoc; i++) {
		if (Numbering->map[Jloc[i]] >=0) { // if free
			if (i>=shift) {
				EqSystem->J[I+J] = Numbering->map[Jloc[i]];
				EqSystem->V[I+J] = Vloc[i];
				J++;
			}
		}

	}



}

void EqSystem_solve(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering)
{
	//int i;
	INIT_TIMER
	TIC


	// Reinitialize Pressure
	int iCell;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = 0;
	}



	if (UPPER_TRI) {
		pardisoSolveSymmetric(EqSystem, Solver, Grid, Physics, BC, Numbering);
	}
	else {
		printf("No solver function for assymmetric matrices\n");
		exit(0);
		//err = pardisoSolveAssymmetric(EqSystem->I ,EqSystem->J, EqSystem->V ,EqSystem->x , EqSystem->b, EqSystem->nEq);
	}

	TOC


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










void EqSystem_initSolver (EqSystem* EqSystem, Solver* Solver)
{

	//int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n
	printf("===== Init Solver =====\n");
	INIT_TIMER
	TIC
	int i;


	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0;
	}



	Solver->mtype = -2;        /* Real symmetric matrix */

	Solver->nrhs = 1;          /* Number of right hand sides. */

	/* Internal solver memory pointer pt,                  */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
	/* or void *pt[64] should be OK on both architectures  */
	//void    *pt[64];

	/* Pardiso control parameters. */
	//int      iparm[64];
	//double   dparm[64];
	//int      maxfct, mnum, phase, error, msglvl, solver;

	int phase;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters.                                */
	/* -------------------------------------------------------------------- */

	int error = 1;
	int solver = 0;// use sparse direct solver
	Solver->msglvl = 0;         // Print statistical information
	pardisoinit (Solver->pt,  &Solver->mtype, &solver, Solver->iparm, Solver->dparm, &error);

	if (error != 0)
	{
		if (error == -10 )
			printf("No license file found \n");
		if (error == -11 )
			printf("License is expired \n");
		if (error == -12 )
			printf("Wrong username or hostname \n");
		exit(0);
	}
	else {
		printf("[PARDISO]: License check was successful ... \n");
	}





	// Numbers of processors, value of OMP_NUM_THREADS
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	printf("Number of procs: %i\n", num_procs);


	// Solver options
	Solver->iparm[2]  = num_procs;

	Solver->iparm[32] = 0; /* compute determinant */
	Solver->iparm[7]  = 1; /* Max numbers of iterative refinement steps. */

	Solver->iparm[27] = 1; // 0: sequential reordering, 1: parallel reordering in METIS
	//Solver->iparm[50] = 0; // 0: openMP, 1:MPI
	//Solver->iparm[51] = 0; // number of compute nodes for MPI


	Solver->iparm[29] = 100; //size of supernodes, default 80

	Solver->maxfct = 1;		// Maximum number of numerical factorizations.
	Solver->mnum   = 1;     // Which factorization to use.

	Solver->msglvl = 0;     // Print statistical information
	error  = 0;         	// Initialize error flag



	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */

	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] += 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] += 1;
	}




	/* -------------------------------------------------------------------- */
	/*  .. pardiso_chk_matrix(...)                                          */
	/*     Checks the consistency of the given matrix.                      */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */
	if  (DEBUG) {
		printf("--  chkmatrix\n");
		pardiso_chkmatrix  (&Solver->mtype, &EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &error);
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
		pardiso_chkvec (&EqSystem->nEq, &Solver->nrhs, EqSystem->b, &error);
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
		pardiso_printstats (&Solver->mtype, &EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &Solver->nrhs, EqSystem->b, &error);
		if (error != 0) {
			printf("\nERROR right hand side: %d", error);
			exit(1);
		}
	}
	/* -------------------------------------------------------------------- */
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */

	phase = 11;

	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error, Solver->dparm);

	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors  = %d", Solver->iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d\n", Solver->iparm[18]);









	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] -= 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] -= 1;
	}

	TOC
	printf("Solver initialization: %.3f s\n", toc);



}


void pardisoSolveSymmetric(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering)
{



	INIT_TIMER
	int i, phase;
	double   	ddum;              // Double dummy
	int      	idum;              // Integer dummy.
	int 		error;


	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0;
	}

	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */


	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] += 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] += 1;
	}



	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */

	if (TIMER) {
		TIC
	}

	phase = 22;

	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);

	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//printf("Factorization completed ...\n ");

	if (TIMER) {
		TOC
		printf("Phase 22 - Numerical factorization: %.3f s\n", toc);
	}



	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	if (TIMER) {
		TIC
	}


	phase = 33;

	// Solve full system Vx, Vy, P


	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);




	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}







	if (TIMER) {
		TOC
		printf("Phase 33 - Back substitution: %.3f s\n", toc);
	}
	if  (DEBUG) {
		printf("\nThe solution of the system is: ");

		for (i = 0; i < EqSystem->nEq; i++) {
			printf("\n x [%d] = % f", i, EqSystem->x[i] );
		}
	}






	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] -= 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] -= 1;
	}







}




void EqSystem_computeNormResidual(EqSystem* EqSystem)
{
	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq;
	int J,i;
	EqSystem->normResidual = 0;


	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		Residual[iEq] = EqSystem->b[iEq];
	}

	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		//Residual = EqSystem->b[iEq];
		for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {

			J = EqSystem->J[i];
			Residual[iEq] += - (EqSystem->V[i]*EqSystem->x[J]);
			if (UPPER_TRI) {
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);
			}
			//printf("%.3e\n",EqStokes.x[iEq] );
		}

	}
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		EqSystem->normResidual += Residual[iEq]*Residual[iEq];
	}
	EqSystem->normResidual = sqrt(EqSystem->normResidual);

	// compute the norm of b
	compute norm_b = 0;
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		norm_b += EqSystem->b[iEq]*EqSystem->b[iEq];
	}
	norm_b = sqrt(norm_b);

	// Normalize the residual
	EqSystem->normResidual /= norm_b; // Normalize the norm of the residual by the norm of the right hand side



	free(Residual);
}

