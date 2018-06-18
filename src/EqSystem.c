/*
 * eqSystem.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"





void EqSystem_Memory_allocateI (EqSystem* EqSystem)
{
	EqSystem->I 			= (int*) malloc((EqSystem->nEq+1)  * sizeof(int));
}

void EqSystem_Memory_allocate(EqSystem* EqSystem)
{
	EqSystem->J = (int*)     malloc(EqSystem->nnz * sizeof(int));
	EqSystem->V = (compute*) malloc(EqSystem->nnz * sizeof(compute));
	EqSystem->b = (compute*) malloc( EqSystem->nEq * sizeof(compute));
#if (PENALTY_METHOD)
	EqSystem->b0= (compute*) malloc( EqSystem->nEq * sizeof(compute));
#endif
	EqSystem->x = (compute*) malloc( EqSystem->nEq * sizeof(compute));
	EqSystem->S = (compute*) malloc( EqSystem->nEq * sizeof(compute));

	int i;
	for (i = 0; i < EqSystem->nEq; ++i) {
		EqSystem->S[i] = 1.0;
	}

	EqSystem->penaltyFac = 1e5;

}

void EqSystem_Memory_free(EqSystem* EqSystem, Solver* Solver)
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
	free(EqSystem->S);
}

void EqSystem_assemble(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering, bool updateScale, Numerics* Numerics)
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



	// Fill J, V and EqSystem->b for Free and Dirichlet nodes
	// ===============================================
	if (DEBUG) {
		printf("Start Filling loop\n");
		printf("nEq = %i, EqSystem->nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	}
	int iEq, ix, iy, I;
	int Ic = 0;
	int Iloc, IBC;
	StencilType Stencil;
	//int INumMap;

	int order[13] = {0,1,2,3,4,5,6,7,8,9,10,11,12};

	compute scale;



#pragma omp parallel for private(iEq, I, ix, iy, i, Stencil, order, nLoc, Ic, Jloc, Vloc, bloc, shift, J,  Iloc, IBC, scale) OMP_SCHEDULE
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
		for (i=0;i<13;i++) {
			order[i] = i;
		}

		// Call the required Stencil function and fill Jloc, Vloc, bloc, etc...
		LocalStencil_Call(Stencil, order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift, &nLoc, &Ic, Numerics);


		// ===========================================
		// Fill the right hand side and apply BC
		// ===========================================

		// Fill right hand side with the local right hand side
		EqSystem->b[iEq] = bloc;
		for (i=0; i<nLoc; i++) {
			Iloc = Numbering->map[Jloc[order[i]]];

			if (Iloc < 0) { // if Boundary node
				IBC = abs(Iloc) - 1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC];
				}
				else  if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Vloc[order[Ic]]  += -Vloc[order[i]]; // +1 to VxC
					EqSystem->b[iEq] += -Vloc[order[i]] * 2.0*BC->value[IBC];

				}
				else if (BC->type[IBC]==Neumann) { // Neumann

					switch (Stencil) {
					case Stencil_Stokes_Darcy_Momentum_x:
					case Stencil_Stokes_Momentum_x:

						//printf("Mx, Ic = %i, i = %i\n",Ic, i);
						if 		(i==1) { // VxW
							Vloc[order[Ic]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//+ BC->value[IBC]/Grid->dx;//-Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						}
						else if (i==3) { // VxE
							Vloc[order[Ic]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						else if (i==5) { // VySW
							Vloc[order[6]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						else if (i==6) { // VySE
							Vloc[order[5]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						else if (i==7) { // VyNW
							Vloc[order[8]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						else if (i==8) { // VyNE
							Vloc[order[7]]  += + Vloc[order[i]]; // +1 to VxCp
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						break;

					case Stencil_Stokes_Darcy_Momentum_y:
					case Stencil_Stokes_Momentum_y:
						//printf("My, Ic = %i, i = %i\n",Ic, i);
						if 		(i==4) { // VyS
							Vloc[order[Ic]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//+ BC->value[IBC]/Grid->dy;//-Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						}
						else if (i==8) { // VyN
							Vloc[order[Ic]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC]/Grid->dy;//+Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxS-1];
						}
						else if (i==0) { // VxSW
							Vloc[order[1]]  += + Vloc[order[i]];
							EqSystem->b[iEq] += 0.0;
						}
						else if (i==1) { // VxSE
							Vloc[order[0]]  += + Vloc[order[i]];
							EqSystem->b[iEq] += 0.0;
						}
						else if (i==2) { // VxNW
							Vloc[order[3]]  += + Vloc[order[i]];
							EqSystem->b[iEq] += 0.0;
						}
						else if (i==3) { // VxNE
							Vloc[order[2]]  += + Vloc[order[i]];
							EqSystem->b[iEq] += 0.0;
						}
						break;

					case Stencil_Stokes_Continuity:
					case Stencil_Stokes_Darcy_Continuity:
						//printf("C., Ic = %i, i = %i\n",Ic, i);
						if 		(i==0) { // VxW
							Vloc[order[1]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//+ BC->value[IBC];///Grid->dx;//-Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						}
						else if (i==1) { // VxE
							Vloc[order[0]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC];///Grid->dx;// - BC->value[IBC]/Grid->dx;//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}

						if 		(i==2) { // VyS
							Vloc[order[3]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//+ BC->value[IBC];///Grid->dy;//-Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						}
						else if (i==3) { // VyN
							Vloc[order[2]]  += + Vloc[order[i]]; // +1 to VxC
							EqSystem->b[iEq] += 0.0;//- BC->value[IBC];///Grid->dy;//+Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxS-1];
						}
						break;

					case Stencil_Stokes_Darcy_Darcy:
					case Stencil_Poisson:
					case Stencil_Heat:
						break;
					}

				}
				else if (BC->type[IBC]==NeumannGhost) { // NeumannGhost
					Vloc[order[Ic]] += Vloc[order[i]]; // +1 to VxC



					switch (Stencil) {
					case Stencil_Stokes_Darcy_Momentum_x:
					case Stencil_Stokes_Momentum_x:
						if 		(i==0) { // VxS
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dy;
							//EqSystem->b[iEq] += + BC->value[IBC]/Grid->DYEC[0];//-Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						}
						else if (i==4) { // VxN
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dy;
							//EqSystem->b[iEq] += - BC->value[IBC]/Grid->DYEC[Grid->nyS-1];//+Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
						break;

					case Stencil_Stokes_Darcy_Momentum_y:
					case Stencil_Stokes_Momentum_y:
						if 		(i==5) { // VyW
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dx;
							//EqSystem->b[iEq] += + BC->value[IBC]/Grid->DXEC[0];//-Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						}
						else if (i==7) { // VyE
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dx;
							//EqSystem->b[iEq] += - BC->value[IBC]/Grid->DXEC[Grid->nxS-1];//+Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxS-1];
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
							//EqSystem->b[iEq] += +BC->value[IBC] / Grid->DYEC[0];
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						} else if (i==1) { // W
							//EqSystem->b[iEq] += +BC->value[IBC] / Grid->DXEC[0];
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						} else if (i==3) { // E
							//EqSystem->b[iEq] += -BC->value[IBC] / Grid->DXEC[Grid->nxEC-2];
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxEC-2];
						} else if (i==4) { // N
							//EqSystem->b[iEq] += -BC->value[IBC] / Grid->DYEC[Grid->nyEC-2];
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyEC-2];							
						}

						// For the moment only 0 gradient is implement
						// This section should be filled to account for a given gradient
						break;
					}

				}
				else if (BC->type[IBC]==Infinity) { // NeumannGhost

					switch (Stencil) {
					case Stencil_Stokes_Darcy_Momentum_x:
					case Stencil_Stokes_Momentum_x:
					case Stencil_Stokes_Darcy_Momentum_y:
					case Stencil_Stokes_Momentum_y:
					case Stencil_Stokes_Continuity:
					case Stencil_Stokes_Darcy_Continuity:
					case Stencil_Stokes_Darcy_Darcy:
					case Stencil_Poisson:
						printf("error: infinity like BC not implemented yet for the stencil #%i",Stencil);
						break;

					case Stencil_Heat:
						//printf("koko!! =======================================================================\n");

						if (i==0) { // S
							Vloc[order[Ic]]  += Vloc[order[i]] * BC->DeltaL/(BC->DeltaL+Grid->DYEC[0]); // +1 to VxC
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0]/(BC->DeltaL+Grid->DYEC[0]);
						} else if (i==1) { // W
							Vloc[order[Ic]]  += Vloc[order[i]] * BC->DeltaL/(BC->DeltaL+Grid->DXEC[0]); // +1 to VxC
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0]/(BC->DeltaL+Grid->DXEC[0]);
						} else if (i==3) { // E
							Vloc[order[Ic]]  += Vloc[order[i]] * BC->DeltaL/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]); // +1 to VxC
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxEC-2]/(BC->DeltaL+Grid->DXEC[Grid->nxEC-2]);
						} else if (i==4) { // N
							Vloc[order[Ic]]  += Vloc[order[i]] * BC->DeltaL/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]); // +1 to VxC
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyEC-2]/(BC->DeltaL+Grid->DYEC[Grid->nyEC-2]);
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



		// ===========================================
			// Compute the scaling factor
			// ===========================================
			if (updateScale) {
				//printf("in There, WTF!!!!\n");

				if (Ic == -1) { // 0 in the diagonal
					EqSystem->S[iEq] = 1.0;
				} else {
#if (DARCY)
					scale = 1.0;//1.0/sqrt(fabs(Vloc[order[Ic]]));
#else
					scale = 1.0/sqrt(fabs(Vloc[order[Ic]]));
#endif
					//printf("iEq = %i, Vloc = %.2e, scale = %.2e\n",iEq, Vloc[order[Ic]], scale );
					if (scale<1e-8) {
						EqSystem->S[iEq] = 1.0;
					} else {
						EqSystem->S[iEq] = scale;
					}
				}
			}
	} // end of the equation loop

#if (!PENALTY_METHOD)
	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
			EqSystem->J[EqSystem->I[i]] = i;
			EqSystem->V[EqSystem->I[i]] = 0.0;
		}
	}
#endif
}

























void EqSystem_check(EqSystem* EqSystem)
{
	int i, j, I;
	// Check
/*
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

	printf(" ===== b =====\n");
	for (i=0;i<EqSystem->nEq;i++) {
		printf("%.2e  ", EqSystem->b[i]);
	}
	printf(" \n");

	printf(" ===== S =====\n");
	for (i=0;i<EqSystem->nEq;i++) {
		printf("%.2e  ", EqSystem->S[i]);
	}
	printf(" \n");

	*/


	/*
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
	*/

	/*
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

	printf("===== V per row=====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*.2e",10,EqSystem->V[I+j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("===== RHS =====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("RHS[%i] = %.2e\n", i, EqSystem->b[i]);
	}
	printf("\n");

	printf("=====  x  =====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("x[%i] = %.2e\n", i, EqSystem->x[i]);
	}
	printf("\n");

}






void EqSystem_solve(EqSystem* EqSystem, Solver* Solver, BC* BC, Numbering* Numbering, Model* Model)
{
	//int i;
	INIT_TIMER
	TIC
	if (UPPER_TRI) {

#if (PENALTY_METHOD)
		if (EqSystem == &(Model->EqStokes)) {
			pardisoSolveSymmetric_Penalty(EqSystem, Solver, BC, Numbering, Model);
		} else {
			pardisoSolveSymmetric(EqSystem, Solver, BC, Numbering, Model);
		}
#else
		pardisoSolveSymmetric(EqSystem, Solver, BC, Numbering, Model);
#endif
	}
	else {
		printf("No solver function for assymmetric matrices\n");
		exit(0);
	}

	TOC
	printf("Direct solve: %.3f s\n", toc);

}










void EqSystem_initSolver (EqSystem* EqSystem, Solver* Solver)
{
	//int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n
	printf("===== Init Solver =====\n");
	INIT_TIMER
	TIC
	int i;

	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0.0;
	 	EqSystem->b[i] = 0.0;
#if (PENALTY_METHOD)
	 	EqSystem->b0[i] = 0.0;
#endif
	}
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->V[i] = 0.0;
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
	Solver->iparm[0] = 0;
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
	printf("phase 11, nEqIni = %i, nEq = %i, nnz= %i\n", EqSystem->nEqIni, EqSystem->nEq, EqSystem->nnz);
	//exit(0);
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
	printf("\nNumber of factorization MFLOPS = %d", Solver->iparm[18]);
	printf("\nPeak memory usage during analysis and symbolic factorizaiton = %d MB\n", Solver->iparm[14]/1024);



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


void pardisoSolveSymmetric(EqSystem* EqSystem, Solver* Solver, BC* BC, Numbering* Numbering, Model* Model)
{


	INIT_TIMER
	int i, phase;
	double   	ddum;              // Double dummy
	int      	idum;              // Integer dummy.
	int 		error;

	EqSystem_scale(EqSystem);
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
	EqSystem_unscale(EqSystem);


}



void pardisoSolveStokesAndUpdatePlasticity(EqSystem* EqSystem, Solver* Solver, BC* BC, Numbering* Numbering, Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics		= &(Model->Physics);
	MatProps* MatProps		= &(Model->MatProps);
	Numerics* Numerics		= &(Model->Numerics);
#if NON_LINEAR_VISU
		Visu* Visu = &(Model->Visu);
#endif

	INIT_TIMER
	int i, phase;
	double   	ddum;              // Double dummy
	int      	idum;              // Integer dummy.
	int 		error;




	compute* b_VE = (compute*) malloc(EqSystem->nEq * sizeof(compute));
	compute* NonLin_x0 = (compute*) malloc(EqSystem->nEq * sizeof(compute));
	compute* NonLin_b0 = (compute*) malloc(EqSystem->nEq * sizeof(compute));
	compute* NonLin_dx = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	// ===== get EffStrainRate =====
	// ===== get EffStrainRate =====
	int iEq, iy, ix, iCell;
	
	




	int Method = Numerics->yieldComputationType;
	bool useVEPGuess = true;
		// initial guess
	compute minL = 1e30;
	if (Numerics->timeStep>0 && useVEPGuess) {
		Physics_Eta_ZandLambda_updateGlobal(Model);
	}
	printf("minL = %.2e\n",minL);



	EqSystem_assemble(EqSystem, Grid, BC, Physics, Numbering, true, Numerics);
	EqSystem_scale(EqSystem);

	for (iEq=0; iEq<EqSystem->nEq; iEq++) {
		b_VE[iEq] = EqSystem->b[iEq];
	}

	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0.0;
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

	//if (TIMER) {
		TIC
	//}



	
	if (Method == 0) {

		phase = 22;
		pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
				&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
				Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);

		if (Numerics->timeStep>0 && useVEPGuess) {
			EqSystem_ApplyRHSPlasticity(Model, b_VE);
		}

		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			exit(2);
		}
		//printf("Factorization completed ...\n ");

		//if (TIMER) {
			TOC
			printf("Phase 22 - Numerical factorization: %.3f s\n", toc);
		//}
		phase = 33;
	} else {
		phase = 23;
	}
	
	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	//if (TIMER) {
		TIC
	//}


	int Counter = 0;
	EqSystem->normResidual = 1e100;
	compute tol = Numerics->absoluteTolerance;
	int maxCounter = Numerics->maxNonLinearIter;

	Numerics->lsLastRes = 1e100;

	

	while (EqSystem->normResidual>tol && Counter<maxCounter) {
		
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			NonLin_x0[iEq] = EqSystem->x[iEq]; 
		}
		pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
					&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
					Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);

#pragma omp parallel for private(iEq) OMP_SCHEDULE
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			NonLin_b0[iEq] = EqSystem->b[iEq]; 
			NonLin_dx[iEq] = EqSystem->x[iEq] - NonLin_x0[iEq]; 
		}
		for (i = 0; i < EqSystem->nEq+1; i++) {
			EqSystem->I[i] -= 1;
		}
		for (i = 0; i < EqSystem->nnz; i++) {
			EqSystem->J[i] -= 1;
		}
		
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			if (isnan(EqSystem->x[iEq]) ) {
				printf("nan in solution\n");
				exit(0);
			}
		}
		

		int iLS = 0;
		Numerics->lsGlob = 1.0;
		Counter++;
		Numerics->minRes = 1E100;	
		Numerics->oldRes = EqSystem->normResidual;
		// Line Search

		

		while (iLS < Numerics->nLineSearch+1) {
			//printf("iLs = %i, Numerics->nLineSearch = %i\n", iLS, Numerics->nLineSearch);
#pragma omp parallel for private(iEq) OMP_SCHEDULE
			for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
				EqSystem->x[iEq] = NonLin_x0[iEq] + Numerics->lsGlob*(NonLin_dx[iEq]);
				EqSystem->b[iEq] = NonLin_b0[iEq];
			}

			

			// Unscale the solution vector
	#pragma omp parallel for private(i) OMP_SCHEDULE
			for (i=0; i<EqSystem->nEq; ++i) {
				EqSystem->x[i] *= EqSystem->S[i];
				EqSystem->b[i] /= EqSystem->S[i];
			}
			Physics_Velocity_retrieveFromSolution(Model);
			Physics_P_retrieveFromSolution(Model);
			// Re-scale the solution vector
	#pragma omp parallel for private(i) OMP_SCHEDULE
			for (i=0; i<EqSystem->nEq; ++i) {
				EqSystem->x[i] /= EqSystem->S[i];
				EqSystem->b[i] *= EqSystem->S[i];
			}


			for (iy = 0; iy<Grid->nyS; iy++) {
				for (ix = 0; ix<Grid->nxS; ix++) {
					int iNode = ix + iy*Grid->nxS;
					if (isnan(Physics->ZShear[iNode])) {
						printf("Zshear nan found\n");
					}
				}
			}
			for (i = 0; i<Grid->nVyTot; ++i) {
				if (isnan(Physics->Vy[i])) {
					printf("EqSys, nan in Vy\n");
					
				}
			}


#if (1)
			//if (Numerics->yieldComputationType==2) {
			//	Physics_Eta_computeLambda_FromParticles_updateGlobal(Model, false);
			//} else {
			Physics_Eta_ZandLambda_updateGlobal(Model);
			//} 

			

			// ===== Apply the correction to the right hand side vector =====
			
			if (Method == 0) {
				EqSystem_ApplyRHSPlasticity(Model, b_VE);
			} else {
				EqSystem_assemble(EqSystem, Grid, BC, Physics, Numbering, false, Numerics);
				EqSystem_scale(EqSystem);
			}
#endif
			
			// ===== Apply the correction to the right hand side vector =====


			// Do stuff =====================================
			EqSystem_computeNormResidual(EqSystem);
			printf("LS: backSubs %i: a = %.3f,  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", Counter-1, Numerics->lsGlob, fabs(EqSystem->normResidual-Numerics->oldRes), EqSystem->normResidual);

			if (EqSystem->normResidual<Numerics->minRes) {
				Numerics->minRes = EqSystem->normResidual;
				Numerics->lsBestGlob = Numerics->lsGlob;
			}
			iLS++;
			if (iLS<Numerics->nLineSearch) {
				Numerics->lsGlob = Numerics->lsGlob/2.0;
			} else {
				if (Numerics->lsGlob == Numerics->lsBestGlob) {
					break;
				} else {
					Numerics->lsGlob = Numerics->lsBestGlob;
				}
			}

			if (iLS == 1 && Numerics->minRes<Numerics->lsLastRes) {
				break;
			}
			if (EqSystem->normResidual>1e10) {
				break;
			}
			if (isnan(EqSystem->normResidual) || isinf(EqSystem->normResidual)) {
				printf("\n\n\n\n error: Something went wrong. The norm of the residual is NaN\n");
				for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
					if (isnan(EqSystem->b[iEq])) {
						printf("nan in b!\n");
					}
				}

				EqSystem_computeNormResidual(EqSystem);
				printf("With b_VE: backSubs %i: a = %.3f,  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", Counter-1, Numerics->lsGlob, fabs(EqSystem->normResidual-Numerics->oldRes), EqSystem->normResidual);

				break;
			}
			
		} // end of line search
		Numerics->lsLastRes = EqSystem->normResidual;

		for (i = 0; i < EqSystem->nEq+1; i++) {
			EqSystem->I[i] += 1;
		}
		for (i = 0; i < EqSystem->nnz; i++) {
			EqSystem->J[i] += 1;
		}
		
		// anti-Numerics->stalling
		if (fabs(EqSystem->normResidual-Numerics->oldRes)<EqSystem->normResidual*Numerics->relativeTolerance) {
			break;
			Numerics->stalling = true;
			Numerics->stallingCounter++;
		} else {
			Numerics->stalling = false;
			Numerics->stallingCounter = 0;
		}

#if NON_LINEAR_VISU

		Visu->update = true;
		Visu->updateGrid = false;
		Visu_main(Model);
		if (glfwWindowShouldClose(Visu->window))
			break;
			
#endif
		
	}

	Numerics->lsGlob = 1.0;

	free(b_VE);
	free(NonLin_x0);
	free(NonLin_b0);
	free(NonLin_dx);


	// =========================================================
	// 				Apply the plastic correction
	// =========================================================



	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}




	//if (TIMER) {
		TOC
		printf("Phase 33 - Back substitution and plastic corr: %.3f s\n", toc);
	//}
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

	EqSystem_computeNormResidual(EqSystem);
	printf("backSubs %i: |F|/|b|: %.2e\n", Counter-1, EqSystem->normResidual);

	EqSystem_unscale(EqSystem);

}




void EqSystem_computeNormResidual(EqSystem* EqSystem)
{
	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq;
	int J,i;
	EqSystem->normResidual = 0.0;


#pragma omp parallel for private(iEq, i, J) OMP_SCHEDULE
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

//#pragma omp parallel for private(iEq, i, J) OMP_SCHEDULE
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
				J = EqSystem->J[i];
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}
		}

	}

	compute norm_b = 0.0;
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		EqSystem->normResidual += Residual[iEq]*Residual[iEq];
		norm_b += EqSystem->b[iEq]*EqSystem->b[iEq];
	}
	norm_b = sqrt(norm_b);
	EqSystem->normResidual = sqrt(EqSystem->normResidual);


	// Normalize the residual
	EqSystem->normResidual /= norm_b; // Normalize the norm of the residual by the norm of the right hand side
	EqSystem->norm_b = norm_b; // Normalize the norm of the residual by the norm of the right hand side


	free(Residual);
}



void EqSystem_scale(EqSystem* EqSystem) {

	int i,j; // matrix index
	int locNNZ;
	int I, J; // index in I, J, V


	// Scale b
	for (i=0; i<EqSystem->nEq; ++i) {
		EqSystem->b[i] *= EqSystem->S[i];
	}


	//printf("===========   Scaling   ================\n");

	// Scale A
	int C = 0;
	for (i=0; i<EqSystem->nEq; ++i) {
		I = EqSystem->I[i];
		locNNZ = (EqSystem->I[i+1] - EqSystem->I[i]);
		for (J = 0; J < locNNZ; ++J) {
			j = EqSystem->J[I+J];
			EqSystem->V[I + J] *=  EqSystem->S[i] * EqSystem->S[j];
			if (I+J!=C) {
				printf("I+J = %i, C = %i\n",I+J, C);
				exit(0);
			}
			C++;
		}
	}

}




void EqSystem_ApplyRHSPlasticity(Model* Model, compute* b_VE) {

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics		= &(Model->Physics);
	Numbering* Numbering	= &(Model->NumStokes);
	EqSystem* EqSystem		= &(Model->EqStokes);

	int nxEC = Grid->nxEC;
	int nxS = Grid->nxS;
	int iCell, iNode;
	compute dxC = Grid->dx;
	compute dyC = Grid->dy;
	int iEq, ix, iy, i;
	StencilType Stencil;


	compute dVxdy, dVydx, dVxdx, dVydy;

	compute* Exx_VE_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy_VE_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute dt = Physics->dt;
#pragma omp parallel for private(iy,ix, iCell, dVxdx, dVydy) OMP_SCHEDULE
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			compute G = Physics->G[iCell];
			Exx_VE_CellGlobal[iCell] = 0.5*(dVxdx-dVydy) + Physics->sigma_xx_0[iCell]/(2.0*G*dt);
		}
	}

#pragma omp parallel for private(iy,ix, iNode, dVxdy, dVydx) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;

			compute G = Physics->GShear[iNode];
			Exy_VE_NodeGlobal[iNode] = 0.5*(dVxdy+dVydx) + Physics->sigma_xy_0[iNode]/(2.0*G*dt);
		}
	}



#pragma omp parallel for private(iEq, iy, ix, i, Stencil) OMP_SCHEDULE
		for (iEq=0; iEq<EqSystem->nEq; iEq++) {

			ix = Numbering->IX[iEq];
			iy = Numbering->IY[iEq];

			i = 1;
			while (iEq>=Numbering->subEqSystem0[i]) {
				i++;
			}
			Stencil = Numbering->Stencil[i-1];
		

			if (Stencil==Stencil_Stokes_Momentum_x)		{
				int NormalPeriod = 0;
				int NormalE = ix+1   + (iy)*nxEC;
				int NormalW = ix     + (iy)*nxEC + NormalPeriod;
				int ShearN = ix      + iy*nxS;
				int ShearS = ix      + (iy-1)*nxS;


				compute Tau_p_xxE,Tau_p_xxW, Tau_p_xyN, Tau_p_xyS;


				Tau_p_xxE = 2.0 * (1.0 - Physics->Lambda[NormalE]    ) * Physics->Z[NormalE]*Exx_VE_CellGlobal[NormalE];
				Tau_p_xxW = 2.0 * (1.0 - Physics->Lambda[NormalW]    ) * Physics->Z[NormalW]*Exx_VE_CellGlobal[NormalW];
				Tau_p_xyN = 2.0 * (1.0 - Physics->LambdaShear[ShearN]) * Physics->ZShear[ShearN]*Exy_VE_NodeGlobal[ShearN];
				Tau_p_xyS = 2.0 * (1.0 - Physics->LambdaShear[ShearS]) * Physics->ZShear[ShearS]*Exy_VE_NodeGlobal[ShearS];

				EqSystem->b[iEq]  = b_VE[iEq] + EqSystem->S[iEq] * (  ( Tau_p_xxE  -   Tau_p_xxW)/dxC  +  ( Tau_p_xyN  -  Tau_p_xyS)/dyC );

				if (isnan(EqSystem->b[iEq])) {
					printf("Momentum x b is nan, b_VE[iEq] = %.2e S = %.2e, Tau_p_xxE = %.2e, Tau_p_xxW = %.2e, Tau_p_xyN = %.2e, Tau_p_xyS = %.2e, dxC = %.2e, dyC = %.2e \n", b_VE[iEq], EqSystem->S[iEq], Tau_p_xxE, Tau_p_xxW, Tau_p_xyN, Tau_p_xyS, dxC, dyC);
				}

			}
			else if (Stencil==Stencil_Stokes_Momentum_y) 	{
				int NormalN = ix      + (iy+1)*nxEC ;
				int NormalS = ix      + (iy)*nxEC ;
				int ShearE  = ix      + iy*nxS    ;
				int ShearW  = ix-1    + iy*nxS    ;

				compute Tau_p_yyN,Tau_p_yyS, Tau_p_xyE, Tau_p_xyW;


				Tau_p_yyN = - 2.0 * (1.0 - Physics->Lambda[NormalN]    ) * Physics->Z[NormalN]*Exx_VE_CellGlobal[NormalN]; // i.e. -Tau_xx
				Tau_p_yyS = - 2.0 * (1.0 - Physics->Lambda[NormalS]    ) * Physics->Z[NormalS]*Exx_VE_CellGlobal[NormalS]; // i.e. -Tau_xx
				Tau_p_xyE =   2.0 * (1.0 - Physics->LambdaShear[ShearE]) * Physics->ZShear[ShearE]*Exy_VE_NodeGlobal[ShearE];
				Tau_p_xyW =   2.0 * (1.0 - Physics->LambdaShear[ShearW]) * Physics->ZShear[ShearW]*Exy_VE_NodeGlobal[ShearW];

				EqSystem->b[iEq] = b_VE[iEq] + EqSystem->S[iEq] * ( ( Tau_p_yyN  -   Tau_p_yyS)/dyC  +  ( Tau_p_xyE  -  Tau_p_xyW)/dxC );

				if (isnan(EqSystem->b[iEq])) {
					printf("Momentum x b is nan, b_VE[iEq] = %.2e S = %.2e, Tau_p_yyN = %.2e, Tau_p_yyS = %.2e, Tau_p_xyE = %.2e, Tau_p_xyW = %.2e, dxC = %.2e, dyC = %.2e \n", b_VE[iEq], EqSystem->S[iEq], Tau_p_yyN, Tau_p_yyS, Tau_p_xyE, Tau_p_xyW, dxC, dyC);
				}
			}
			else if (Stencil==Stencil_Stokes_Continuity) 	{
				// do nothing
			}
			// re-scale solution and right hand side
		}



		free(Exx_VE_CellGlobal);
		free(Exy_VE_NodeGlobal);

}


































#if (PENALTY_METHOD)
void pardisoSolveSymmetric_Penalty(EqSystem* EqSystem, Solver* Solver, BC* BC, Numbering* Numbering, Model* Model)
{



	INIT_TIMER
	int i, phase;
	double   	ddum;              // Double dummy
	int      	idum;              // Integer dummy.
	int 		error;


	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0.0;
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


		EqSystem->maxDivVel = 1;
		compute tolerance = 1.0E-10;
		int maxUzawa = 10;
		int it = 0;

		for (i = 0; i < Grid->nCTot; ++i) {
			Physics->P[i] = 0.0;
		}

		for (i = 0; i < EqSystem->nEq; ++i) {
			EqSystem->b0[i] = EqSystem->b[i];
		}

		// Initialize Vx, Vy, P
		//int i;
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

			TIC
			EqSystem_computePressureAndUpdateRHS(EqSystem, Grid, Numbering, Physics, BC);
			pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
					&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
					Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);

			Physics_Velocity_retrieveFromSolution(Model);
			if (it==0)
				EqSystem->maxDivVel = 1;

			TOC
			printf("Solve:%.2fs\n",toc);

			printf("maxDivVel=%.3e\n", EqSystem->maxDivVel);


			it++;
		}



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

#endif







void EqSystem_unscale(EqSystem* EqSystem) {
	int i;
	for (i=0; i<EqSystem->nEq; ++i) {
		EqSystem->x[i] *= EqSystem->S[i];
		EqSystem->b[i] /= EqSystem->S[i];
	}
}

compute EqSystem_maxDivVel(Model *Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	int ix, iy;
	compute dx, dy, divV;
	compute maxDivV = 0.0;
	for (iy = 1; iy < Grid->nyEC - 1; ++iy)
	{
		for (ix = 1; ix < Grid->nxEC - 1; ++ix)
		{
			dx = Grid->DXS[ix - 1];
			dy = Grid->DYS[iy - 1];
			divV = (Physics->Vx[ix + iy * Grid->nxVx] - Physics->Vx[ix - 1 + iy * Grid->nxVx]) / dx;
			divV += (Physics->Vy[ix + iy * Grid->nxVy] - Physics->Vy[ix + (iy - 1) * Grid->nxVy]) / dy;

			maxDivV = fmax(maxDivV, divV);

		}
	}

	return maxDivV;

}


#if (PENALTY_METHOD)
void EqSystem_computePressureAndUpdateRHS(EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC)
{
	int ix, iy;
	int iEq;

	compute divV, dx, dy;;
	EqSystem->maxDivVel = 0.0;
	compute K = EqSystem->penaltyFac;
	//printf("=== divVel ===\n");
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			dx = Grid->DXS[ix - 1];
			dy = Grid->DYS[iy - 1];
			divV = (Physics->Vx[ix + iy * Grid->nxVx] - Physics->Vx[ix - 1 + iy * Grid->nxVx]) / dx;
			divV += (Physics->Vy[ix + iy * Grid->nxVy] - Physics->Vy[ix + (iy - 1) * Grid->nxVy]) / dy;

			EqSystem->maxDivVel = fmax(EqSystem->maxDivVel, divV);
			Physics->P[ix+iy*Grid->nxC] += K*divV;
		}
	}
	//Physics->P[Grid->nCTot-1] = 0; //Add a dirichlet condition

	int i;
	int PW, PE, PN, PS;
	StencilType Stencil;

	//printf("=== divVel loc ===\n");
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];

		i = 1;
		while (iEq>=Numbering->subEqSystem0[i]) {
			i++;
		}
		Stencil = Numbering->Stencil[i-1];

		if (Stencil==Stencil_Stokes_Momentum_x)		{//if (iEq<EqSystem->VyEq0) { // Vx equations
			//if (BC->isNeu[ix+iy*Grid->nxVx]==false) { // If Free equation (i.e. not Neumann equation)
				PW = (ix-1) + (iy-1)*Grid->nxC;
				PE = (ix  ) + (iy-1)*Grid->nxC;
				// Note: should be ok for periodic as well
				EqSystem->b[iEq] = EqSystem->b0[iEq] + (Physics->P[PE] - Physics->P[PW])/Grid->dx; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
			//}
		} else if (Stencil==Stencil_Stokes_Momentum_y) {
			//if (BC->isNeu[ix+iy*Grid->nxVy+Grid->nVxTot]==false) { // If Free equation (i.e. not Neumann equation)
				PN = (ix-1) + (iy  )*Grid->nxC;
				PS = (ix-1) + (iy-1)*Grid->nxC;
				EqSystem->b[iEq] = EqSystem->b0[iEq] + (Physics->P[PN] - Physics->P[PS])/Grid->dy; // /!\ the minus one is there to take care of the fortran indexing (used by Pardiso)
			//}
		}
	}
}
#endif
