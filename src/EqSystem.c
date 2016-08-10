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


	// Init variables
	// ===============================
	int nLoc;
	int i, J;



	int SetupType = BC->SetupType;





	int Jloc[13];
	compute Vloc[13];
	compute bloc;

	int shift = 0;




	// Reinitialize b
	/*
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->J[i] = -1;
	}

	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->b[i] = 0.0;
	}
*/




	// Fill J, V and EqSystem->b for Free and Dirichlet nodes
	// ===============================================
	if (DEBUG) {
		printf("Start Filling loop\n");
		printf("nEq = %i, EqSystem->nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	}
	int iEq, ix, iy, I;
	int IC = 0;
	int Iloc, IBC;
	StencilType Stencil;
	//int INumMap;

	int order[13] = {0,1,2,3,4,5,6,7,8,9,10,11,12};

//#pragma omp parallel for private(iEq, I, ix, iy, i, Stencil, order, nLoc, IC, Jloc, Vloc, bloc, shift, J,  Iloc, IBC) schedule(static,32)
	for (iEq=0; iEq<EqSystem->nEq; iEq++) {

		I = EqSystem->I[iEq];
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];

		i = 1;
		while (iEq>=Numbering->subEqSystem0[i]) {
			i++;
		}
		Stencil = Numbering->Stencil[i-1];


		// Reinitialize order
		for (i=0;i<11;i++) {
			order[i] = i;
		}


		//printf("iEq = %i, Stencil #%i\n",iEq,Stencil);
		// Call the required Stencil function and fill Jloc, Vloc, bloc, etc...
		LocalStencil_Call(Stencil, order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &IC);





		// ===========================================
		// Fill the right hand side and apply BC
		// ===========================================

		// Fill right hand side with the local right hand side
		EqSystem->b[iEq] = bloc;

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
					EqSystem->b[iEq] += -Vloc[order[i]] * 2.0*BC->value[IBC];

				}
				else if (BC->type[IBC]==NeumannGhost) { // NeumannGhost
					Vloc[order[IC]] += Vloc[order[i]]; // +1 to VxC



					switch (Stencil) {
					case Stencil_Stokes_Darcy_Momentum_x:
					case Stencil_Stokes_Momentum_x:
						if 		(i==0) { // VxS
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						}
						else if (i==4) { // VxN
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						break;

					case Stencil_Stokes_Darcy_Momentum_y:
					case Stencil_Stokes_Momentum_y:
						if 		(i==5) { // VyW
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						}
						else if (i==7) { // VyE
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxS-1];
						}
						break;

					case Stencil_Stokes_Continuity:
						break;

					case Stencil_Stokes_Darcy_Continuity:
						break;

					case Stencil_Stokes_Darcy_Darcy:
					case Stencil_Poisson:
					case Stencil_Heat:

						if (i==0) { // S
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						} else if (i==1) { // W
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						} else if (i==3) { // E
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxEC-2];
						} else if (i==4) { // N
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyEC-2];
						}

						// For the moment only 0 gradient is implement
						// This section should be filled to account for a given gradient
						break;
					}

				}
				else {
					printf("error: unknown boundary type\n");
					exit(0);
				}

			}
		}


		// ===========================================
		// Fill J and V global
		// ===========================================

		// Jloc gets the numbering with Dirichlet;
		// When values are inputted in Jsparse they are transformed to NumMap[Jloc[i]]
		// i.e. numbering without dirichlet
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


		/*
		printf("iEq = %i, IX = %i, IY = %i\n",iEq, ix, iy);
		int i;
		for (i = 0; i < 11; ++i) {
			printf("Jloc[%i] = %i\n",i, Numbering->map[Jloc[i]]);
		}
		*/
		/*
		for (i = 0; i < 4; ++i) {
			printf("Jloc[%i] = %i\n",i, Jloc[i]);
		}
		*/




	} // end of the equation loop




	printf("nEq = %i, nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
			EqSystem->J[EqSystem->I[i]] = i;
			EqSystem->V[EqSystem->I[i]] = 0.0;
		}
	}
	//printf("nEq = %i, nRow = %i\n", EqSystem->nEq, EqSystem->nRow);

#if (DEBUG)
	EqSystem_check(EqSystem);
#endif



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
	printf("%.2e  \n", EqSystem->V[0]);
	for (i=0;i<EqSystem->nnz;i++) {
		printf("%.2e  ", EqSystem->V[i]);
	}
	printf(" \n");

	printf(" ===== b=====\n");
	for (i=0;i<EqSystem->nEq;i++) {
		printf("%.2e  ", EqSystem->b[i]);
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
	}


	TOC

	printf("Direct solve: %.2f s\n", toc);



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
	 	EqSystem->b[i] = 0;
	}
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->V[i] = 0;
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
		printf("\nThe solution of the system is: \n");

		for (i = 0; i < EqSystem->nEq; i++) {
			printf(" x [%d] = %.2e\n", i, EqSystem->x[i] );
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


#pragma omp parallel for private(iEq, i, J) schedule(static,32)
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		Residual[iEq] = EqSystem->b[iEq];
		for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
			J = EqSystem->J[i];
			Residual[iEq] += - (EqSystem->V[i]*EqSystem->x[J]);
			/*
			if (UPPER_TRI) {
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}*/
		}
	}

	if (UPPER_TRI) {

#pragma omp parallel for private(iEq, i, J) schedule(static,32)
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
				J = EqSystem->J[i];
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}
		}

	}

	compute norm_b = 0;
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		EqSystem->normResidual += Residual[iEq]*Residual[iEq];
		norm_b += EqSystem->b[iEq]*EqSystem->b[iEq];
	}
	norm_b = sqrt(norm_b);
	EqSystem->normResidual = sqrt(EqSystem->normResidual);


	// Normalize the residual
	EqSystem->normResidual /= norm_b; // Normalize the norm of the residual by the norm of the right hand side



	free(Residual);
}

