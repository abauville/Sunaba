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
	EqSystem->x0 = (compute*) malloc( EqSystem->nEq * sizeof(compute));
	EqSystem->dx = (compute*) malloc( EqSystem->nEq * sizeof(compute));
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
	free(EqSystem->x0);
	free(EqSystem->dx);
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
			fill_J_V_local(Type, ix, iy, I, iEq, EqSystem, Grid, Numbering, Physics, BC);
		}

	}
	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		if (!EqSystem->penaltyMethod) {
			for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
				EqSystem->J[EqSystem->I[i]] = i;
				EqSystem->V[EqSystem->I[i]] = 0.0;
			}
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
	/*
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

	 */
	/*
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
	 */
}





void fill_J_V_local(int Type, int ix, int iy,int I, int iEq, EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC)
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

	int BCtype = BC->SetupType;

	compute EtaN, EtaS, EtaW, EtaE;

	compute dx = Grid->dx;
	compute dy = Grid->dy;

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

	nxN = nxC;

	nxS = nxC+1;

	int Jloc[11];
	compute Vloc[11];
	compute bloc;

	if (Type==0)
	{
		// =========================================================================
		// =========================================================================
		//
		//                                   Vx
		//
		// =========================================================================
		// =========================================================================

		if (EqSystem->penaltyMethod) {
			nLoc = 9; // Maximum number of non zeros for Stokes on the staggered grid
		}
		else {
			nLoc = 11; // Maximum number of non zeros for Stokes on the staggered grid
		}

		if (UPPER_TRI) {
			shift = 2;
		}
		else {
			shift = 0;
		}

		if (BCtype==0) {

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
			Jloc[ 5] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySW
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

			EtaN    = Physics->etaShear[ShearN];
			EtaS    = Physics->etaShear[ShearS];
			EtaE    = Physics->eta[NormalE];
			EtaW    = Physics->eta[NormalW];


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

			bloc = - Physics->g[0] * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );

			if (EqSystem->penaltyMethod) {
				// contribution of PE
				Vloc[3] += -EqSystem->penaltyFac/dx/dx;// VxE
				Vloc[2] +=  EqSystem->penaltyFac/dx/dx;// VxC
				Vloc[8] += -EqSystem->penaltyFac/dy/dx;// VyNE
				Vloc[6] +=  EqSystem->penaltyFac/dy/dx;// VySE

				// contribution of PW
				Vloc[2] +=  EqSystem->penaltyFac/dx/dx;// VxC
				Vloc[1] += -EqSystem->penaltyFac/dx/dx;// VxW
				Vloc[7] +=  EqSystem->penaltyFac/dy/dx;// VyNW
				Vloc[5] += -EqSystem->penaltyFac/dy/dx;// VySW
			}

		}
		// Adjustment for the case where boundary conditions are periodic horizontally
		else if (BCtype==1) {
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
				// Get Viscosities
				// ================
				NormalE = ix      + (iy-1)*nxN;
				ShearN  = ix      + iy*nxS    ;
				ShearS  = ix      + (iy-1)*nxS;

				EtaN    = Physics->etaShear[ShearN];
				EtaS    = Physics->etaShear[ShearS];
				EtaE    = Physics->eta[NormalE];

				NormalW = ix-1    + (iy-1)*nxN + nxN; // **
				EtaW    = Physics->eta[NormalW]; 	// **

				bloc = - Physics->g[0]  * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );

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

				if (EqSystem->penaltyMethod) {
					// contribution of PE
					Vloc[2] += -EqSystem->penaltyFac/dx/dx;// VxE
					Vloc[1] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[8] += -EqSystem->penaltyFac/dx/dy;// VyNE
					Vloc[6] +=  EqSystem->penaltyFac/dx/dy;// VySE

					// contribution of PW
					Vloc[1] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[3] += -EqSystem->penaltyFac/dx/dx;// VxW
					Vloc[7] +=  EqSystem->penaltyFac/dx/dy;// VyNW
					Vloc[5] += -EqSystem->penaltyFac/dx/dy;// VySW
				}




			}


			else if (ix==nxVx-2) {
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
				// Get Viscosities
				// ================
				NormalE = ix      + (iy-1)*nxN;
				NormalW = ix-1    + (iy-1)*nxN;
				ShearN  = ix      + iy*nxS    ;
				ShearS  = ix      + (iy-1)*nxS;

				EtaN    = Physics->etaShear[ShearN];
				EtaS    = Physics->etaShear[ShearS];
				EtaE    = Physics->eta[NormalE];
				EtaW    = Physics->eta[NormalW];
				bloc = - Physics->g[0]  * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );

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

				if (EqSystem->penaltyMethod) {
					// contribution of PE
					Vloc[1] += -EqSystem->penaltyFac/dx/dx;// VxE
					Vloc[3] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[7] += -EqSystem->penaltyFac/dx/dy;// VyNE
					Vloc[5] +=  EqSystem->penaltyFac/dx/dy;// VySE

					// contribution of PW
					Vloc[3] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[2] += -EqSystem->penaltyFac/dx/dx;// VxW
					Vloc[8] +=  EqSystem->penaltyFac/dx/dy;// VyNW
					Vloc[6] += -EqSystem->penaltyFac/dx/dy;// VySW
				}




			}

			else {
				// Same as BCType == 0

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
				Jloc[ 5] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySW
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

				EtaN    = Physics->etaShear[ShearN];
				EtaS    = Physics->etaShear[ShearS];
				EtaE    = Physics->eta[NormalE];
				EtaW    = Physics->eta[NormalW];


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

				bloc = - Physics->g[0] * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );

				if (EqSystem->penaltyMethod) {
					// contribution of PE
					Vloc[3] += -EqSystem->penaltyFac/dx/dx;// VxE
					Vloc[2] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[8] += -EqSystem->penaltyFac/dy/dx;// VyNE
					Vloc[6] +=  EqSystem->penaltyFac/dy/dx;// VySE

					// contribution of PW
					Vloc[2] +=  EqSystem->penaltyFac/dx/dx;// VxC
					Vloc[1] += -EqSystem->penaltyFac/dx/dx;// VxW
					Vloc[7] +=  EqSystem->penaltyFac/dy/dx;// VyNW
					Vloc[5] += -EqSystem->penaltyFac/dy/dx;// VySW
				}
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

		if (EqSystem->penaltyMethod) {
			nLoc = 9; // Maximum number of non zeros for Stokes on the staggered grid
		}
		else {
			nLoc = 11; // Maximum number of non zeros for Stokes on the staggered grid
		}
		if (UPPER_TRI) {
			shift = 6;
		}
		else {
			shift = 0;
		}
		if (BCtype==0) {
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


			EtaN    = Physics->eta[NormalN];
			EtaS    = Physics->eta[NormalS];
			EtaE    = Physics->etaShear[ShearE];
			EtaW    = Physics->etaShear[ShearW];


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

			bloc = - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );


			if (EqSystem->penaltyMethod) {
				// contribution of PN
				Vloc[8] += -EqSystem->penaltyFac/dy/dy;// VyN
				Vloc[6] +=  EqSystem->penaltyFac/dy/dy;// VyC
				Vloc[3] += -EqSystem->penaltyFac/dx/dy;// VxNE
				Vloc[2] +=  EqSystem->penaltyFac/dx/dy;// VxNW

				// contribution of PS
				Vloc[6] +=  EqSystem->penaltyFac/dy/dy;// VyC
				Vloc[4] += -EqSystem->penaltyFac/dy/dy;// VyS
				Vloc[1] +=  EqSystem->penaltyFac/dx/dy;// VxSE
				Vloc[0] += -EqSystem->penaltyFac/dx/dy;// VxSW
			}

		}



		else if (BCtype==1) {
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

				ShearE  = ix      + iy*nxS    ;
				NormalN = ix-1    + (iy  )*nxN + nxN; //**
				NormalS = ix-1    + (iy-1)*nxN + nxN; //**
				ShearW  = ix-1    + iy*nxS   + nxS-1 ; // **

				EtaN    = Physics->eta[NormalN];
				EtaS    = Physics->eta[NormalS];
				EtaE    = Physics->etaShear[ShearE];
				EtaW    = Physics->etaShear[ShearW];




				bloc 	= - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );

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

				if (EqSystem->penaltyMethod) {
					// contribution of PN
					Vloc[8] += -EqSystem->penaltyFac/dy/dy;// VyN
					Vloc[5] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[2] += -EqSystem->penaltyFac/dx/dy;// VxNE
					Vloc[3] +=  EqSystem->penaltyFac/dx/dy;// VxNW

					// contribution of PS
					Vloc[5] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[4] += -EqSystem->penaltyFac/dy/dy;// VyS
					Vloc[0] +=  EqSystem->penaltyFac/dx/dy;// VxSE
					Vloc[1] += -EqSystem->penaltyFac/dx/dy;// VxSW
				}


			}
			else if (ix==nxVy-3) {
				if (UPPER_TRI) {
					shift = 7;
				}

				// Fill Jloc: list of all J indices (including Dirichlet)
				// ================================================================
				Jloc[ 0] =   ix      + (iy  )*nxVx - 1                   ; // VxSW
				Jloc[ 1] =   ix      + (iy  )*nxVx                       ; // VxSE
				Jloc[ 2] =   ix      + (iy+1)*nxVx - 1                   ; // VxNW
				Jloc[ 3] =   ix      + (iy+1)*nxVx                       ; // VxNE
				Jloc[ 4] =   ix      + iy*nxVy     + nVxTot    - nxVy    ; // VyS
				Jloc[ 5] =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE **
				Jloc[ 6] =   ix      + iy*nxVy     + nVxTot   - 1        ; // VyW **
				Jloc[ 7] =   ix      + iy*nxVy     + nVxTot              ; // VyC **
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


				EtaN    = Physics->eta[NormalN];
				EtaS    = Physics->eta[NormalS];
				EtaE    = Physics->etaShear[ShearE];
				EtaW    = Physics->etaShear[ShearW];


				Vloc[ 0] =  EtaW/dy/dx;
				Vloc[ 1] = -EtaE/dy/dx;
				Vloc[ 2] = -EtaW/dy/dx;
				Vloc[ 3] =  EtaE/dy/dx;
				Vloc[ 4] =  2.0 * EtaS/dy/dy;
				Vloc[ 5] =  EtaE/dx/dx; //**
				Vloc[ 6] =  EtaW/dx/dx; //**
				Vloc[ 7] = -2.0 * EtaN/dy/dy   -2.0 * EtaS/dy/dy   -1.0 * EtaE/dx/dx   -1.0 * EtaW/dx/dx; //**
				Vloc[ 8] =  2.0 * EtaN/dy/dy;
				Vloc[ 9] =  1.0/dy;
				Vloc[10] = -1.0/dy;



				bloc 	= - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );



				if (EqSystem->penaltyMethod) {


					// contribution of PN
					Vloc[8] += -EqSystem->penaltyFac/dy/dy;// VyN
					Vloc[7] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[3] += -EqSystem->penaltyFac/dx/dy;// VxNE
					Vloc[2] +=  EqSystem->penaltyFac/dx/dy;// VxNW

					// contribution of PS
					Vloc[7] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[4] += -EqSystem->penaltyFac/dy/dy;// VyS
					Vloc[1] +=  EqSystem->penaltyFac/dx/dy;// VxSE
					Vloc[0] += -EqSystem->penaltyFac/dx/dy;// VxSW
				}



			}

			else {
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


				EtaN    = Physics->eta[NormalN];
				EtaS    = Physics->eta[NormalS];
				EtaE    = Physics->etaShear[ShearE];
				EtaW    = Physics->etaShear[ShearW];


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

				bloc = - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );


				if (EqSystem->penaltyMethod) {
					// contribution of PN
					Vloc[8] += -EqSystem->penaltyFac/dy/dy;// VyN
					Vloc[6] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[3] += -EqSystem->penaltyFac/dx/dy;// VxNE
					Vloc[2] +=  EqSystem->penaltyFac/dx/dy;// VxNW

					// contribution of PS
					Vloc[6] +=  EqSystem->penaltyFac/dy/dy;// VyC
					Vloc[4] += -EqSystem->penaltyFac/dy/dy;// VyS
					Vloc[1] +=  EqSystem->penaltyFac/dx/dy;// VxSE
					Vloc[0] += -EqSystem->penaltyFac/dx/dy;// VxSW
				}

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
		if (BCtype==0) {
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

		}

		else if (BCtype==1) {

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
			else {
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

			}
		}

		bloc = 0;


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

	EqSystem->b[iEq] += bloc;

	J = 0;
	for (i=0; i<nLoc; i++) {
		if (Numbering->map[Jloc[i]] !=-1) { // if not Dirichlet
			if (i>=shift) {
				EqSystem->J[I+J] = Numbering->map[Jloc[i]];
				EqSystem->V[I+J] = Vloc[i];
				J++;
			}
		}
		else { // if it is a Dirichlet

			// Loop through the Dirichlet list to find the index inside the list of the equation corresponding to
			// Jloc[i]. The index is then used to retrieve the boundary value in BC_Value_Dir
			iDir = 0;
			while (BC->listDir[iDir] !=Jloc[i]) {

				iDir++;
				if (iDir==BC->nDir) {
					fprintf(stderr,"Couldn't find the Dirichlet equations");
					exit(0);
				}
			}


			EqSystem->b[iEq] += -Vloc[i] * BC->valueDir[iDir];


		}
	}



}

void EqSystem_solve(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering)
{
	//int i;
	INIT_TIMER
	TIC

	if (!EqSystem->penaltyMethod) {
		// Reinitialize Pressure
		int iCell;
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			Physics->P[iCell] = 0;
		}
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


	if (EqSystem->penaltyMethod) {
		Solver->mtype =  -2;        /* Real symmetric matrix */
	}
	else {
		Solver->mtype = -2;        /* Real symmetric matrix */
	}
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

	printf("===== Enter the solver function =====\n");


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

	if (!EqSystem->penaltyMethod) {
		// Solve full system Vx, Vy, P


		pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
				&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
				Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);

		/*
		Physics_set_VxVyP_FromSolution(Physics, Grid, BC, Numbering, EqSystem);
		EqSystem_computePressureAndUpdateRHS(EqSystem, Grid, Numbering, Physics, BC);
		printf("maxDivVel=%.3e\n", EqSystem->maxDivVel);
		 */
	}
	else {

		// Vx, Vy and update P through an Uzawa loop
		EqSystem->maxDivVel = 1;
		compute tolerance = 1.0E-10;
		int maxUzawa = 10;
		int it = 0;

		for (i = 0; i < Grid->nCTot; ++i) {
			Physics->P[i] = 0;
		}


		EqSystem->b0 = (compute*) malloc(EqSystem->nEq*sizeof(compute));
		for (i = 0; i < EqSystem->nEq; ++i) {
			EqSystem->b0[i] = EqSystem->b[i];
		}

		// Initialize Vx, Vy, P
		int i;
		for (i = 0; i < Grid->nVxTot; ++i) {
			Physics->Vx[i] = 0;
		}
		for (i = 0; i < Grid->nVyTot; ++i) {
			Physics->Vy[i] = 0;
		}
		for (i = 0; i < Grid->nCTot; ++i) {
			Physics->P[i] = 0;
		}


		printf("Uzawa iterations\n");
		while (EqSystem->maxDivVel>tolerance && it<maxUzawa) {
			//while (it<maxUzawa) {

			phase = 33;
			for (i=0; i<EqSystem->nEq; i++) {
				//EqSystem->x[i] = 0;
			}
			TIC


			EqSystem_computePressureAndUpdateRHS(EqSystem, Grid, Numbering, Physics, BC);
			pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
					&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
					Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);

			Physics_set_VxVyP_FromSolution(Physics, Grid, BC, Numbering, EqSystem);
			if (it==0)
				EqSystem->maxDivVel = 1;

			TOC
			printf("Solve:%.2fs\n",toc);
			/*
			printf("===== RHS =====\n");
			for (i=0; i<EqSystem->nEq; i++) {
				printf("RHS[%i] = %.2e\n", i, EqSystem->b[i]);
			}
			printf("\n");


			 */

			printf("maxDivVel=%.3e\n", EqSystem->maxDivVel);
			/*
			printf("=== Pressure ===\n");
			int ix, iy;
			for (iy = 0; iy < Grid->nyC; ++iy) {
				for (ix = 0; ix < Grid->nxC; ++ix) {
					printf("%.3f  ",Physics->P[ix+iy*Grid->nxC]);
				}
				printf("\n");
			}
			 */



			it++;
		}
		free(EqSystem->b0);

	}


	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}







	printf("Solve completed\n");
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





	/*
	phase = -1;                 // Release internal memory.
		pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
				 &EqSystem->nEq, &ddum, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
				 Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);
	 */


}



void EqSystem_computePressureAndUpdateRHS(EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC)
{
	int ix, iy;
	int iEq;
	int VxC, VxE, VxW, VyNE, VySE, VyNW, VySW;
	int VyC, VyN, VyS, VxNE, VxSE, VxNW, VxSW;
	compute divVelE, divVelW, divVelN, divVelS;
	compute divVel;
	EqSystem->maxDivVel = 0;
	compute K = EqSystem->penaltyFac;
	compute dum, dum1, dum2;
	//printf("=== divVel ===\n");
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			VxE = ix+1 + (iy+1)*Grid->nxVx;
			VxW = ix   + (iy+1)*Grid->nxVx;
			VyN = ix+1 + (iy+1)*Grid->nxVy;
			VyS = ix+1 + (iy  )*Grid->nxVy;
			divVel = (Physics->Vx[VxE]-Physics->Vx[VxW])/Grid->dx + (Physics->Vy[VyN]-Physics->Vy[VyS])/Grid->dy;
			dum1 = (Physics->Vx[VxE]-Physics->Vx[VxW])/Grid->dx ;
			dum2 = (Physics->Vy[VyN]-Physics->Vy[VyS])/Grid->dy;
			//Physics->P[ix+iy*Grid->nxC] += EqSystem->penaltyFac*divVel;
			Physics->P[ix+iy*Grid->nxC] += K*divVel;
			//printf("ix=%i, iy=%i,  , VyN=%i, VyS=%i, VxE=%i, VxW=%i, Physics->Vx[VxE]=%.2f, Physics->Vx[VxW]=%.2f, Physics->Vy[VyN]=%.2f, Physics->Vy[VyS] = %.2f, dum1 = %.2e, dum2 = %.2e, sum = %.2e, P=%.4e\n",ix, iy, VyN, VyS, VxE, VxW, Physics->Vx[VxE], Physics->Vx[VxW], Physics->Vy[VyN], Physics->Vy[VyS], dum1, dum2, dum1+dum2, Physics->P[ix+iy*Grid->nxC]);
			//printf("%.2e   ", divVel);
			if (fabs(divVel)>EqSystem->maxDivVel)
				EqSystem->maxDivVel = fabs(divVel);


		}
		//printf("\n");
		Physics->P[Grid->nCTot-1] = 0; //Add a dirichlet condition
	}

	int IX[2], IY[2];
	int i;
	compute Coeff[2];
	int PW, PE, PN, PS;
	//printf("=== divVel loc ===\n");
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];

		if (iEq<EqSystem->VyEq0) { // Vx equations
			if (BC->isNeu[ix+iy*Grid->nxVx]==false) { // If Free equation (i.e. not Neumann equation)
				VxC 	= (ix  )   	+ (iy  ) * Grid->nxVx;
				VxE 	= (ix+1)   	+ (iy  ) * Grid->nxVx;
				VxW 	= (ix-1)   	+ (iy  ) * Grid->nxVx;
				VyNE 	= (ix+1) 	+ (iy  ) * Grid->nxVy;
				VySE 	= (ix+1) 	+ (iy-1) * Grid->nxVy;
				VyNW 	= (ix  ) 	+ (iy  ) * Grid->nxVy;
				VySW 	= (ix  ) 	+ (iy-1) * Grid->nxVy;

				PW = (ix-1) + (iy-1)*Grid->nxC;
				PE = (ix  ) + (iy-1)*Grid->nxC;

				if (BC->SetupType==1) {
					if (ix==0) {
						VxW += Grid->nxVx-1;
						PW  += Grid->nxC;
					}
				}

				divVelE = (Physics->Vx[VxE] - Physics->Vx[VxC])/Grid->dx + (Physics->Vy[VyNE] - Physics->Vy[VySE])/Grid->dy;
				divVelW = (Physics->Vx[VxC] - Physics->Vx[VxW])/Grid->dx + (Physics->Vy[VyNW] - Physics->Vy[VySW])/Grid->dy;

				//printf("ix = %i, iy = %i, divVelE = %.2e , divVelW = %.2e\n", ix, iy, divVelE, divVelW);



				//EqSystem->b[iEq-1] -= (Physics->P[PE] - Physics->P[PW])/Grid->dx; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				//EqSystem->b[iEq-1] += (divVelE-divVelW)/Grid->dx*EqSystem->penaltyFac; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				//EqSystem->b[iEq] += (divVelE-divVelW)/Grid->dx*K; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				EqSystem->b[iEq] = EqSystem->b0[iEq] + (Physics->P[PE] - Physics->P[PW])/Grid->dx; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				//printf("iEq=%i, PE=%.2e, PW=%.2e", iEq, Physics->P[PE], Physics->P[PW]);
				//printf("iEq=%i, ix=%i, iy=%i, VxW=%i, b=%.4e, dum=%.4e, DivVelE*K=%.4e, DivVelW*K=%.4e, PE=%.4e, PW=%.4e\n",iEq, ix, iy, VxW, EqSystem->b[iEq], dum, divVelE*K, divVelW*K, Physics->P[PE], Physics->P[PW]);
			}
		}


		else{
			if (BC->isNeu[ix+iy*Grid->nxVy+Grid->nVxTot]==false) { // If Free equation (i.e. not Neumann equation)
				VyC 	= (ix  )   	+ (iy  ) * Grid->nxVy;
				VyN 	= (ix  )   	+ (iy+1) * Grid->nxVy;
				VyS 	= (ix  )   	+ (iy-1) * Grid->nxVy;
				VxNE 	= (ix  ) 	+ (iy+1) * Grid->nxVx;
				VxSE 	= (ix  ) 	+ (iy  ) * Grid->nxVx;
				VxNW 	= (ix-1) 	+ (iy+1) * Grid->nxVx;
				VxSW 	= (ix-1) 	+ (iy  ) * Grid->nxVx;

				if (BC->SetupType==1) {
					if (ix==0) {
						VxNW += Grid->nxVx-1;
						VxSW += Grid->nxVx-1;
					}
				}

				PN = (ix-1) + (iy  )*Grid->nxC;
				PS = (ix-1) + (iy-1)*Grid->nxC;

				divVelN = (Physics->Vy[VyN] - Physics->Vy[VyC])/Grid->dy    +    (Physics->Vx[VxNE] - Physics->Vx[VxNW])/Grid->dx;

				divVelS = (Physics->Vy[VyC] - Physics->Vy[VyS])/Grid->dy    +    (Physics->Vx[VxSE] - Physics->Vx[VxSW])/Grid->dx;

				//printf("ix = %i, iy = %i, divVelN = %.2e , divVelS = %.2e , VyC=%i, VyN=%i, VyS=%i, VxNE=%i, VxNW=%i, VxSE=%i, VxSW=%i, Physics->Vy[VyN]=%.2f , Physics->Vy[VyC]=%.2f, Physics->Vx[VxNE]=%.2f, Physics->Vx[VxNW]=%.2f   \n", ix, iy, divVelN, divVelS, VyC, VyN, VyS, VxNE, VxNW, VxSE, VxSW, Physics->Vy[VyN], Physics->Vy[VyC] , Physics->Vx[VxNE], Physics->Vx[VxNW]);


				//EqSystem->b[iEq-1] -= (Physics->P[PE] - Physics->P[PW])/Grid->dx; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				//EqSystem->b[iEq-1] += (divVelN-divVelS)/Grid->dy*EqSystem->penaltyFac; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				//EqSystem->b[iEq] += (divVelN-divVelS)/Grid->dy*K; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
				EqSystem->b[iEq] = EqSystem->b0[iEq] + (Physics->P[PN] - Physics->P[PS])/Grid->dy; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
			}
		}
	}
}


void EqSystem_computeNormResidual(EqSystem* EqSystem)
{
	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq;
	int I,J,i;
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
		}

	}
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		EqSystem->normResidual += Residual[iEq]*Residual[iEq];
	}
	EqSystem->normResidual = sqrt(EqSystem->normResidual);

	free(Residual);
}

