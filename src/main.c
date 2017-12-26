/*
 ============================================================================
 Name        : main.c
 Author      : Arthur Bauville
 Version     :
 Copyright   : 
 Description :
 ============================================================================
 */

#include "stokes.h"



int main(int argc, char *argv[]) {

	printf("\n\n\n\n\n\n");
	printf("               ============================\n"
			"               ============================\n");
	printf("\n\n\n\n\n\nBeginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());

#if (DARCY)
	int iy, ix, iCell;
#endif


	// Declare structures
	// =================================
	
	Model Model;
	// Declare structures
	// =================================
	// General
	Grid* Grid 				= &(Model.Grid);
	MatProps* MatProps 		= &(Model.MatProps);
	Particles* Particles 	= &(Model.Particles);
	Physics* Physics 		= &(Model.Physics);
	Char* Char 				= &(Model.Char);
	// Stokes: Conservation of momentum + continuity (+- Darcy)
	Numbering* NumStokes 	= &(Model.NumStokes);
	BC* BCStokes 			= &(Model.BCStokes);
	EqSystem* EqStokes				= &(Model.EqStokes);
	Solver* SolverStokes 	= &(Model.SolverStokes);

	IC* ICDarcy 			= &(Model.ICDarcy);

	// Heat conservation
	Numbering* NumThermal 	= &(Model.NumThermal);
	IC* ICThermal 			= &(Model.ICThermal);
	BC* BCThermal 			= &(Model.BCThermal);
	EqSystem* EqThermal  	= &(Model.EqThermal);
	Solver* SolverThermal 	= &(Model.SolverThermal);
	// Numerics
	Numerics* Numerics 		= &(Model.Numerics);
	// Visu
#if (VISU)
	Visu* Visu 				= &(Model.Visu);
#endif
	// Input/Output
	Input* Input 			= &(Model.Input);
	Output* Output 			= &(Model.Output);



/*
// ===================================================================
// For reference
	// Declare structures
	// =================================
	// General
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);
	// Stokes: Conservation of momentum + continuity (+- Darcy)
	Numbering* NumStokes 	= &(Model->NumStokes);
	BC* BCStokes 			= &(Model->BCStokes);
	EqSystem* EqStokes				= &(Model->EqStokes);
	Solver* SolverStokes 	= &(Model->SolverStokes);

	IC* ICDarcy 			= &(Model->ICDarcy);

	// Heat conservation
	Numbering* NumThermal 	= &(Model->NumThermal);
	IC* ICThermal 			= &(Model->ICThermal);
	BC* BCThermal 			= &(Model->BCThermal);
	EqSystem* EqThermal  	= &(Model->EqThermal);
	Solver* SolverThermal 	= &(Model->SolverThermal);
	// Numerics
	Numerics* Numerics 		= &(Model->Numerics);
	// Visu
#if (VISU)
	Visu* Visu 				= &(Model->Visu);
#endif
	// Input/Output
	Input* Input 			= &(Model->Input);
	Output* Output 			= &(Model->Output);
// For reference
// ===================================================================
*/





	INIT_TIMER

	strncpy(Input->currentFolder, argv[0],strlen(argv[0])-strlen("StokesFD") );
	//Input->currentFolder[strlen(argv[0])-strlen("StokesFD")] = '\0';
	Input->currentFolder[strlen(argv[0])-strlen("StokesFD")] = '\0';
	printf("Input->currentFolder = %s\n",Input->currentFolder);

	if (argc < 2) {
		strcpy(Input->inputFile,INPUT_FILE);
	} else {
		strcpy(Input->inputFile,argv[1]);
	}

	printf("Reading input\n");
	Input_read(&Model);
#if (LINEAR_VISCOUS)
	if (Numerics->maxNonLinearIter>1) {
		printf("error: you requested %i non linear iterations, however, they are switched off due to LINEAR_VISCOUS==true\n", Numerics->maxNonLinearIter);
		exit(0);
	}
#endif

#if (!HEAT)
	if (BCThermal->TB!=1.0 || BCThermal->TT!=1.0) {
		printf("TB = %.3e, TT = %.3e\n",BCThermal->TB, BCThermal->TT);
		printf("error: you specified non default thermal boundary conditions, however, the heat equation is switched off due to HEAT==false ");
		exit(0);
	}
#endif


#if (VISU)
	printf("Reading Visu input\n");
	Input_readVisu(&Model);
#endif

	printf("Reading input over\n");

	if (Output->nTypes>0) {
		Output_writeInputCopyInOutput(Output, Input);
	}


	Output->counter = 0;
	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                       		 USER INPUT      		   	              	  //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	printf("nTimesteps = %i\n",Numerics->nTimeSteps);
	Grid->isPeriodic = false;
	if (BCStokes->SetupType==Stokes_SimpleShear) {
		Grid->isPeriodic = true;
		Grid->isFixed 	= true;
	}

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Model);
	Physics->dt = 1.0; //i.e. 0.1*Char.time/Char.time
	Numerics->dtPrevTimeStep = 1.0; //i.e. 0.1*Char.time/Char.time
	Physics->epsRef = fabs(BCStokes->backStrainRate);

	if (Physics->epsRef == 0)
		Physics->epsRef = 1E0;

	Physics->maxVx = (Grid->xmax-Grid->xmin)/Physics->epsRef;
	Physics->maxVy = (Grid->ymax-Grid->ymin)/Physics->epsRef;

	// Init grid and particles
	// =================================
	Grid->nCTot  = Grid->nxC*Grid->nyC;

	Grid->nxEC = Grid->nxC+2;
	Grid->nyEC = Grid->nyC+2;
	Grid->nECTot = Grid->nxEC*Grid->nyEC;

	Grid->nxVx 	= Grid->nxC+1; 		Grid->nyVx	= Grid->nyC+2;
	Grid->nxVy 	= Grid->nxC+2;		Grid->nyVy	= Grid->nyC+1;
	Grid->nxS 	= Grid->nxC+1;		Grid->nyS	= Grid->nyC+1;
	Grid->nSTot  = Grid->nxS*Grid->nyS;

	Grid->nVxTot = Grid->nxVx*Grid->nyVx;
	Grid->nVyTot = Grid->nxVy*Grid->nyVy;

	Grid->xmax_ini = Grid->xmax;
	Grid->xmin_ini = Grid->xmin;
	Grid->ymax_ini = Grid->ymax;
	Grid->ymin_ini = Grid->ymin;

	Grid->dx = (Grid->xmax-Grid->xmin)/Grid->nxC;
	Grid->dy = (Grid->ymax-Grid->ymin)/Grid->nyC;

	
	NumThermal->nSubEqSystem 	= 1;
	NumThermal->Stencil[0] = Stencil_Heat;
	EqThermal->nEqIni 		= Grid->nECTot;

#if (DARCY)
	NumStokes->nSubEqSystem 	= 4;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Darcy_Momentum_x; 	// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Darcy_Momentum_y; 	// Vy
	NumStokes->Stencil[2] 	= Stencil_Stokes_Darcy_Darcy;	   	// Pf
	NumStokes->Stencil[3] 	= Stencil_Stokes_Darcy_Continuity; 	// Pc
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot + Grid->nECTot + Grid->nECTot;
#else
#if (PENALTY_METHOD)
	NumStokes->nSubEqSystem 	= 2;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Momentum_x;		// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Momentum_y; 		// Vy
	//NumStokes->Stencil[2]	= Stencil_Stokes_Continuity;		// P
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot;// + Grid->nECTot;
#else
	NumStokes->nSubEqSystem 	= 3;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Momentum_x;		// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Momentum_y; 		// Vy
	NumStokes->Stencil[2]	= Stencil_Stokes_Continuity;		// P
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot + Grid->nECTot;
#endif
#endif


	Particles->nPC 	= Particles->nPCX * Particles->nPCY;
	Particles->n 	= Grid->nCTot*Particles->nPC;

#if (VISU)
	Visu->ntri   	= 2;//Grid->nxC*Grid->nyC*2;//2;//Grid->nxC*Grid->nyC*2;
	Visu->ntrivert 	= Visu->ntri*3;
	Visu->nParticles = Particles->n+ (int) (Particles->n*0.1); // overallocate 5% of the number of particles
#endif

	


	Numerics->oneMoreIt = false;


	if (DEBUG) {
		if (Grid->nCTot>200) {
			printf("error: The system size exceeds the maximum allowed for debugging\n");
			exit(0);
		}
	}

	if (Grid->isPeriodic && Grid->nxC%2!=0) {
		printf("error: When using the periodic boundaries nxC must be even (because of node 'coloring'\n");
		exit(0);
	}


//======================================================================================================
//======================================================================================================
//
//                          				INITIALIZATION
//
	Numerics->timeStep = -1;
	Numerics->itNonLin = -1;

	// Other variables
	// =================================
	int iEq, iLS;

	//Init Grid
	// =================================
	Grid_Memory_allocate(Grid);
	Grid_init(&Model);
	// Init Physics
	// =================================
	printf("Init Physics\n");
	Physics_Memory_allocate	(&Model);

	// Init Numerics
	// =================================
	printf("Init Numerics\n");
	Numerics_init(Numerics);

	// Initialize Particles
	// =================================
	printf("Particles: Init Particles\n");
	Particles_Memory_allocate	(Particles, Grid);
	Particles_initCoord			(Particles, Grid);
	Particles_updateLinkedList	(Particles, Grid, Physics); // in case a ridiculous amount of noise is put on the particle
	Input_assignPhaseToParticles(&Model);
	Particles_initPassive		(Particles, Grid, Physics);

	// Initialize Physics
	// =================================
	printf("Physics: Init Physics\n");
	Interp_All_Particles2Grid_Global	(&Model);
	Physics_Rho_updateGlobal(&Model);

	Physics_Phase_updateGlobal					(&Model);
	//Physics_dt_update(&Model);
	
#if (HEAT)
	IC_T(Physics, Grid, ICThermal, BCThermal);
	Interp_All_Particles2Grid_Global	(&Model);
#endif
#if (DARCY)
	IC_phi(Physics, Grid, Numerics, ICDarcy, MatProps, Particles);
	Interp_All_Particles2Grid_Global	(&Model);
	memcpy(Physics->phi, Physics->phi0, Grid->nECTot * sizeof(compute));
#endif

	Physics_Rho_updateGlobal	(&Model);

	Physics_P_initToLithostatic (&Model);

	Physics_Eta_init(&Model);
	Physics_Eta_smoothGlobal(&Model);
	Physics_dt_update(&Model);
	Physics_Eta_init(&Model);
	Physics_Eta_smoothGlobal(&Model);

#if (DEBUG)
	Physics_check(&Model);
#endif


	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes			(BCStokes , Grid, Physics, EqStokes);
#if (HEAT)
	BC_initThermal			(BCThermal, Grid, Physics, EqThermal);
#endif
	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_Memory_allocateI		(EqStokes);
	Numbering_Memory_allocate(NumStokes, EqStokes, Grid);
	Numbering_init			(BCStokes, Grid, EqStokes, NumStokes, Physics, Numerics);
	printf("EqSystem: init Stokes\n");
	EqSystem_Memory_allocate	(EqStokes );

	printf("Number of Unknowns for Stokes: %i \n", EqStokes->nEq);

#if (HEAT)
	printf("Numbering: init Thermal\n");
	EqSystem_Memory_allocateI		(EqThermal);
	Numbering_Memory_allocate(NumThermal, EqThermal, Grid);
	Numbering_init			(BCThermal, Grid, EqThermal, NumThermal, Physics);
	printf("EqSystem: init Thermal\n");
	EqSystem_Memory_allocate	(EqThermal);

	printf("Number of Unknowns for Heat: %i \n", EqThermal->nEq);
#endif





	


#if (VISU)
	printf("Visu: Init Visu\n");
	Visu_Memory_allocate(Visu, Grid);
	Visu_init(Visu, Grid, Particles, Char, Input);
#endif

	printf("koko, nEq=%i, nVxTot = %i, nVyTot = %i, nECTot = %i\n", EqStokes->nEq, Grid->nVxTot, Grid->nVyTot, Grid->nECTot);
	// Init Solvers
	// =================================
	printf("EqStokes: Init Solver\n");
	EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, false, Numerics); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (EqStokes, SolverStokes);

#if (HEAT)
	printf("EqThermal: Init Solver\n");
	EqSystem_assemble(EqThermal, Grid, BCThermal, Physics, NumThermal, false, Numerics); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (EqThermal, SolverThermal);
#endif



//
//                          				INITIALIZATION
//
//======================================================================================================
//======================================================================================================




	// Update Cell Values with Part
	// =================================
	Interp_All_Particles2Grid_Global(&Model);
	
	Physics_Rho_updateGlobal(&Model);
	
	Physics_P_initToLithostatic 			(&Model);
	
	// Update BC
	// =================================
	printf("BC: Update\n");
	BCStokes->counter = 0;
	BC_updateStokes_Vel(BCStokes, Grid, Physics, true);
#if (DARCY)
	BC_updateStokesDarcy_P(BCStokes, Grid, Physics, true);
#endif

#if (DARCY)
	Physics_Perm_updateGlobal(&Model);
#endif
	Physics_Rho_updateGlobal(&Model);

	

	compute stressFacIni = Numerics->dt_stressFac;

//======================================================================================================
//======================================================================================================
//
//                          				TIME LOOP
//

	Numerics->timeStep = 0;
	Physics->time = 0;

	double timeStepTic;


	compute* NonLin_x0 = (compute*) malloc(EqStokes->nEq * sizeof(compute));
	compute* NonLin_dx = (compute*) malloc(EqStokes->nEq * sizeof(compute));
	
	//printf("Numerics->maxTime = %.2e, Physics->time = %.2e\n",Numerics->maxTime,Physics->time);
	while(Numerics->timeStep!=Numerics->nTimeSteps && Physics->time <= Numerics->maxTime) {
		printf("\n\n\n          ========  Time step %i, t= %3.2e yrs  ========   \n"
					 "              ===================================== \n\n",Numerics->timeStep, Physics->time*Char->time/(3600*24*365));
		printf("time0 = %.2e\n", Physics->time*Model.Char.time);
		Numerics->dtPrevTimeStep = Physics->dt;
		Numerics->dtAlphaCorr = Numerics->dtAlphaCorrIni;

		printf("dt = %.2e yrs, maxVx = %.2e cm/yr, maxVy = %.2e cm/yr, dx/maxVx = %.2e yrs, dy/maxVy = %.2e yrs", Physics->dt*Char->time/(3600*24*365), (Physics->maxVx*Char->length/Char->time) / (0.01/(3600*24*365)), (Physics->maxVy*Char->length/Char->time) / (0.01/(3600*24*365)), (Grid->dx/Physics->maxVx) * Char->time/(3600*24*365), (Grid->dy/Physics->maxVy) * Char->time/(3600*24*365));
#if (VISU)
		timeStepTic = glfwGetTime();
#endif
		Numerics->itNonLin = -1;


#if (HEAT)
		// save the value from the previous time step
		if (Numerics->itNonLin == -1) {
			for (i = 0; i < Grid->nECTot; ++i) {
				Physics->T0[i] = Physics->T[i];
			}
		}
#endif

		/*
		if(Physics->time*Char->time > 5e6 * (3600*24*365.25)) {
			//Physics->g[1] = 0.0;
			BCStokes->backStrainRate = 0.0;//
		}
		*/
		
		

		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										NON-LINEAR ITERATION											//
		Numerics->stalling = false;
		Numerics->stallingCounter = 0;
		EqStokes->normResidual = 1.0;
		Numerics->normRes0 = 1.0;
		Numerics->normResRef = 1.0;
		Numerics->cumCorrection_fac = 0.0;
		Numerics->lsLastRes = 1E15;
		Numerics->lsGlob = 1.00;
		Numerics->lsBestRes = 1e15;
		Numerics->lsBestGlob = 1.0;
#if (VISU)
		Visu->nonLinItisOver = false;
#endif

		Numerics->itNonLin = 0;

		Numerics->lsLastRes = 1E100;

		Numerics->dt_stressFac = stressFacIni;
		Numerics->oneMoreIt = true;
		//Physics_dt_update(&Model);
#if (PLASTIC_CORR_RHS)
		while(Numerics->oneMoreIt) {
#else
		while( ( (( (EqStokes->normResidual > Numerics->absoluteTolerance ) && Numerics->itNonLin<Numerics->maxNonLinearIter ) || Numerics->itNonLin<Numerics->minNonLinearIter)  || Numerics->cumCorrection_fac<=0.999   ) || Numerics->oneMoreIt) {
#endif
			printf("\n\n  ==== Non linear iteration %i ==== \n",Numerics->itNonLin);
			Numerics->oneMoreIt = false;

			// =====================================================================================//
			//																						//
			// 										COMPUTE STOKES									//
			//if (Numerics->itNonLin<=0 || Physics->dt<1e-2 || Physics->dt>1e2) {
			if (Physics->dt<0.5 || Physics->dt>2.0) {
				printf("Rescale\n");
				Char_rescale(&Model, NonLin_x0);
			}
			//Physics_dt_update(&Model);
			//Physics_Eta_Simple_updateGlobal(&Model);
			memcpy(NonLin_x0, EqStokes->x, EqStokes->nEq * sizeof(compute));
			EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, true, Numerics);
			//EqSystem_scale(EqStokes);
			//EqSystem_solve(EqStokes, SolverStokes, BCStokes, NumStokes, &Model);
			pardisoSolveStokesAndUpdatePlasticity(EqStokes, SolverStokes, BCStokes, NumStokes, &Model);
			//EqSystem_unscale(EqStokes);

			//EqSystem_computeNormResidual(EqStokes);
			//printf("afterSol:  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", fabs(EqStokes->normResidual-Numerics->oldRes), EqStokes->normResidual);
			//if (Numerics->itNonLin<=0) {
				Physics_Velocity_retrieveFromSolution(&Model);
				Physics_P_retrieveFromSolution(&Model);
				//Physics_dt_update(&Model);
				//Char_rescale(&Model, NonLin_x0);
				
			//}
			/*
			printf("EqStokes->x[100] = %.2e, scaledEqStokes->x[100] = %.2e\n", EqStokes->x[100], EqStokes->x[100]*Char->velocity);
			printf("Physics->Vx[100] = %.2e, Physics->Vx[100] = %.2e\n", Physics->Vx[100], Physics->Vx[100]*Char->velocity);
			printf("EqStokes->V[100] = %.2e\n", EqStokes->V[100]);
			printf("Z[50] = %.2e, scaledZ[50] = %.2e\n", Physics->Z[50], Physics->Z[50]*Char->viscosity);
			printf("Zshear[50] = %.2e, scaledZshear[50] = %.2e\n", Physics->ZShear[50], Physics->ZShear[50]*Char->viscosity);
			printf("G[50] = %.2e, scaledG[50] = %.2e\n", Physics->G[50], Physics->G[50]*Char->stress);
			printf("Physics->dt = %.2e\n",Physics->dt);
			*/

			// 										COMPUTE STOKES									//
			//																						//
			// =====================================================================================//


#if (VISCOSITY_TYPE==1)
			printf("/!\\ /!\\ LINEAR_VISCOUS==true, Non-linear iterations are ineffective/!\\ \n");
			Physics_Velocity_retrieveFromSolution(&Model);
			Physics_P_retrieveFromSolution(&Model);
			Physics_Rho_updateGlobal(&Model);
			Physics_Eta_updateGlobal(&Model);

			break;
#elif (VISCOSITY_TYPE==2)
			Physics_Velocity_retrieveFromSolution(&Model);
			break;
#else // VISCOSITY_TYPE==0



			// =====================================================================================//
			//																						//
			// 										LINE SEARCH										//
			// Compute dx
			for (iEq = 0; iEq < EqStokes->nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes->x[iEq] - NonLin_x0[iEq];
			}


			Numerics->minRes = 1E100;
			Numerics->lsGlob = 1.0;
			Numerics->lsState = -1;
			iLS = 0;
			Numerics->oldRes = EqStokes->normResidual;



#if (HEAT)
			// =====================================================================================//
			//																						//
			// 										COMPUTE HEAT									//
			//TIC
			Physics_Velocity_retrieveFromSolution(&Model);
			Physics_P_retrieveFromSolution(&Model);


#if (DARCY)
			Physics_Phi_updateGlobal(&Model);
			Physics_Perm_updateGlobal(&Model);
#endif

			Physics_Rho_updateGlobal(&Model);
			Physics_Eta_updateGlobal(&Model);
			printf("Heat assembly and solve\n");
			EqSystem_assemble(EqThermal, Grid, BCThermal, Physics, NumThermal, true, Numerics);

			EqSystem_scale(EqThermal);
			EqSystem_solve(EqThermal, SolverThermal, BCThermal, NumThermal, &Model);
			EqSystem_unscale(EqThermal);
			Physics_T_retrieveFromSolution(&Model);

			//TOC
			printf("Temp Assembly+Solve+Interp: %.3f s\n", toc);


			// 										COMPUTE HEAT									//
			//																						//
			// =====================================================================================//
#endif
#if (!PLASTIC_CORR_RHS)
			while (iLS < Numerics->nLineSearch+1) {
#pragma omp parallel for private(iEq) OMP_SCHEDULE
				for (iEq = 0; iEq < EqStokes->nEq; ++iEq) {
					EqStokes->x[iEq] = NonLin_x0[iEq] + Numerics->lsGlob*(NonLin_dx[iEq]);
				}

				Physics_Velocity_retrieveFromSolution(&Model);
				Physics_P_retrieveFromSolution(&Model);


#if (DARCY)
				Physics_Phi_updateGlobal(&Model);
				Physics_Perm_updateGlobal(&Model);
#endif
				Physics_Rho_updateGlobal(&Model);
				Physics_Eta_Simple_updateGlobal(&Model);
				//Physics_Eta_updateGlobal(&Model);
				//Physics_Eta_FromParticles_updateGlobal(&Model);
				//Physics_Eta_smoothGlobal(&Model);


#if (DEBUG)
				Physics_check(&Model);
#endif
				EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, false, Numerics);

				EqSystem_computeNormResidual(EqStokes);

				printf("a = %.3f,  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", Numerics->lsGlob, fabs(EqStokes->normResidual-Numerics->oldRes), EqStokes->normResidual);

				if (EqStokes->normResidual<Numerics->minRes) {
					Numerics->minRes = EqStokes->normResidual;
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

				if (Numerics->minRes<Numerics->lsLastRes) {
					break;
				}

				if (EqStokes->normResidual>1e10) {
					break;
				}

				if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
					printf("\n\n\n\n error: Something went wrong. The norm of the residual is NaN\n");
					break;
				}



				if (Numerics->lsGoingDown || Numerics->lsGoingUp) { // if the time step is changing there is no need to check other values of lsGlob (in most cases 1.0 is the best),
					// the residual goes up only because the time step is changeing. Once it stabilizes it will go down again.
					//break;
				}


			} // end of line search
#endif

			// 		   								LINE SEARCH										//
			//																						//
			// =====================================================================================//

			Numerics->cumCorrection_fac += Numerics->lsBestGlob;
			Numerics->lsLastRes = EqStokes->normResidual;

/*
			if (Numerics->lsState == -2) {
				//printf("Break!!\n");
				break;
			}



			// anti-Numerics->stalling
			if (fabs(EqStokes->normResidual-Numerics->oldRes)<EqStokes->normResidual*Numerics->relativeTolerance) {
				break;
				Numerics->stalling = true;
				Numerics->stallingCounter++;
			} else {
				Numerics->stalling = false;
				Numerics->stallingCounter = 0;
			}
			*/
			/*
			if (fabs(EqStokes->normResidual-Numerics->oldRes)<Numerics->absoluteTolerance*1e-4) {
				break;
			}
			*/
			



#if NON_LINEAR_VISU
		Visu->update = true;
		Visu->updateGrid = false;
		Visu_main(&Model);
		if (glfwWindowShouldClose(Visu->window))
			break;
#endif


#if (PLASTIC_CORR_RHS)
	compute dtAdv 	= Numerics->CFL_fac_Stokes*Grid->dx/(Physics->maxVx); // note: the min(dx,dy) is the char length, so = 1
	dtAdv 	= fmin(dtAdv,  Numerics->CFL_fac_Stokes*Grid->dy/(Physics->maxVy));
	printf("dtAdv = %.2e, dt = %.2e, lsGlob = %.2e\n", dtAdv, Physics->dt, Numerics->lsGlob);
	if (dtAdv<Physics->dt && Physics->dt>Numerics->dtMin) {
		compute dtOld = Physics->dt;
		Physics_dt_update(&Model);
		if (Physics->dt!=Physics->dt) {
			Numerics->oneMoreIt = true;
		} else {
			Numerics->oneMoreIt = false;
		}
	} else {
		Numerics->oneMoreIt = false;
	}
	
#endif


			


			if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
				printf("\n\n\n\nerror: Something went wrong. The norm of the residual is NaN\n");
				break;
			}

			if (EqStokes->normResidual>1e10) {
				break;
			}

			/*
			// Break if it already went down auite reasonnably and is Numerics->stalling
			if (fabs(EqStokes->normResidual-Numerics->oldRes)<1e-8 && EqStokes->normResidual<100.0*Numerics->absoluteTolerance) {
				//break;
				//compute stressFacOld = Numerics->dt_stressFac;
				//Numerics->dt_stressFac /= 2.0;
				compute dtOld = Physics->dt;
				Physics->dt /= 2.0;
				//printf("too many iter: stressFac: %.2e --> %.2e\n",stressFacOld,Numerics->dt_stressFac);
				printf("too many iter: dt: %.2e --> %.2e\n",dtOld,Physics->dt);
			}

			if (Numerics->itNonLin>1 && fmod(Numerics->itNonLin,20)==0 && EqStokes->normResidual>10.0*Numerics->absoluteTolerance) {
				//compute stressFacOld = Numerics->dt_stressFac;
				//Numerics->dt_stressFac /= 2.0;
				compute dtOld = Physics->dt;
				Physics->dt /= 2.0;
				//printf("too many iter: stressFac: %.2e --> %.2e\n",stressFacOld,Numerics->dt_stressFac);
				printf("too many iter: dt: %.2e --> %.2e\n",dtOld,Physics->dt);
			}
			*/
			if (EqStokes->normResidual>Numerics->absoluteTolerance*10.0) {
			//	Numerics->oneMoreIt = true;
			}
			
			
			if (Numerics->lsGoingDown) {
				//pNumerics->oneMoreIt = true;
				printf("going down\n");
			}
			

			//Numerics->oneMoreIt = false; // for some reasons it stalls sometime
#endif







			Numerics->itNonLin++;
		} // end of non-linear loop


#if (VISU)
		Visu->nonLinItisOver = true;
#endif


		if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
			printf("\n\n\n\nerror: Something went wrong. The norm of the residual is NaN\n");
			break;
		}
		if (EqStokes->normResidual>1e10) {
			break;
		}



#if (VISU)
		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);
#endif



		// 										NON-LINEAR ITERATION 											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//



#if (VISU)
		timeStepTic = glfwGetTime();
#endif

#if (VISCOSITY_TYPE==0)

		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 									INTERPOLATION FROM CELL TO PARTICLES								//

		// update stress on the particles
		// =============================
		Physics_Dsigma_updateGlobal  (&Model);
		//Physics_Sigma0_updateGlobal_fromGrid(&Model);
#if (ADV_INTERP)
		Interp_Stresses_Grid2Particles_Global(&Model);
#endif

#if (DARCY)
		Interp_Phi_Grid2Particles_Global	(&Model);
#endif


#if (HEAT)
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DT[i] = Physics->T[i] - Physics->T0[i];
		}
		Interp_Temperature_Grid2Particles_Global(&Model);

#endif
#if (STRAIN_SOFTENING)
		Interp_Strain_Grid2Particles_Global(&Model);
#endif

		// 									INTERPOLATION FROM CELL TO PARTICLES								//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//

#endif




		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 												OUTPUT AND VISU											//

		// Output
		// =================

		if (Output->nTypes>0 || Output->nPartTypes>0) {
			bool writeOutput = false;
			if (Output->useTimeFrequency) {
				printf("Output->counter*Output->timeFrequency = %.2e, tim = %.2e\n", Output->counter*Output->timeFrequency, Physics->time);
				if (Physics->time>Output->counter*Output->timeFrequency) {
					writeOutput = true;
				} else if (Output->saveFirstStep && Numerics->timeStep == 0) {
					writeOutput = true;
				}
			} else {
				if (((Numerics->timeStep) % Output->frequency)==0) {
					if (Numerics->timeStep>0 || Output->saveFirstStep) {
						writeOutput = true;
					}
				}
			}
			if (writeOutput) {
				printf("Write output ...\n");
				Output_modelState(&Model);
				Output_data(&Model);
				Output_particles(&Model);
				Output->counter++;
			}
		}


#if VISU
		Visu->update = true;
		if (!Grid->isFixed) {
			Visu->updateGrid = true;
		}
		Visu_main(&Model);
		if (glfwWindowShouldClose(Visu->window))
			break;
#endif


		// 												OUTPUT AND VISU											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//











		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 							ADVECTION AND INTERPOLATION	FROM PARTICLES TO CELL							//

#if (ADV_INTERP)	
		// Advect Particles
		// =============================
		printf("Particles: Advect\n");
		Particles_advect(Particles, Grid, Physics);

		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(Particles, Grid, Physics);


		// Inject particles
		// =================================
		if (Grid->isFixed) {
			Particles_injectAtTheBoundaries(Particles, Grid, Physics, MatProps);
		}

		printf("Particle injection\n");
		Particles_injectOrDelete(Particles, Grid);



		// Advect the box and update Particles position if needed
		// =============================
		switch (BCStokes->SetupType) {
		case Stokes_PureShear:
		case Stokes_Sandbox:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				Grid_updatePureShear(&Model);
				Particles_teleportInsideTheDomain(Particles, Grid, Physics);
			}
			break;
		case Stokes_SimpleShear:
			Particles_Periodicize(Particles, Grid);
			break;
		case Stokes_FixedLeftWall:
			break;
		case Stokes_CornerFlow:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				printf("error: For the corner flow setup, Grid->fixedBox must be true. Correct the input file");
				exit(0);
			}
			break;
		case Stokes_WindTunnel:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				printf("error: For the WindTunnel setup, Grid->fixedBox must be true. Correct the input file");
				exit(0);
			}
			break;
		default:
			break;
		}





		

		Particles_switchStickyAir			(Particles, Grid, Physics, Numerics, MatProps, BCStokes);
		Particles_surfaceProcesses			(&Model);
		// Update the Phase matrix
		// =================================
		Physics_Phase_updateGlobal					(&Model);

#if (VISCOSITY_TYPE==0)
		// Update the Physics on the Cells
		// =================================
		printf("Physics: Interp from particles to grid\n");
		Interp_All_Particles2Grid_Global(&Model);
#endif
#if (CRANK_NICHOLSON_VEL || INERTIA)
		if (Numerics->timeStep>0) {
			Physics_Velocity_advectEulerian(&Model);
		} else {
			#if (CRANK_NICHOLSON_VEL || INERTIA)
			Physics_VelOld_POld_updateGlobal(&Model);
			#endif
		}
#endif

		Physics_Rho_updateGlobal(&Model);
		



#if (DARCY)
		compute dx, dy;
		for (iy = 1; iy < Grid->nyEC-1; ++iy) {
			for (ix = 1; ix < Grid->nxEC-1; ++ix) {
				iCell = ix + iy*Grid->nxEC;
				dx = Grid->DXS[ix-1];
				dy = Grid->DYS[iy-1];
				Physics->divV0[iCell]  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/dx;
				Physics->divV0[iCell] += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/dy;
			}
		}
#endif

		// Update BC
		// =================================
		printf("BC: Update\n");
		BCStokes->counter = 0;
		BC_updateStokes_Vel(BCStokes, Grid, Physics, true);
#if (DARCY)
		BC_updateStokesDarcy_P(BCStokes, Grid, Physics, true);
#endif
#if (HEAT)
		BCThermal->counter = 0;
		BC_updateThermal(BCThermal, Grid, Physics, true);
#endif

#else // i.e. if (!ADV_INTERP)
		Physics_Sigma0_updateGlobal_fromGrid(&Model);
#endif // if (ADV_INTERP)



		// 							ADVECTION AND INTERPOLATION FROM PARTICLES TO CELL 							//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//
		printf("timeN = %.2e\n", Physics->time*Model.Char.time);

		Physics->time += Physics->dtAdv;
		Physics->dtAdv0 = Physics->dtAdv;
		Numerics->timeStep++;

		Physics_dt_update(&Model);
		/*
		// Compute the Visco-elastic effective viscosity
		// =================================
		Physics_dt_update(&Model);
		int i;
		for (i=0;i<Grid->nECTot;++i) {
			Physics->P[i] = 1e100;
#if (DARCY)
			Physics->Pc[i] = 1e100;
			Physics->Pf[i] = 1e100;
#endif
		}
		Physics_Eta_updateGlobal(&Model);
		*/
		
		
		Physics_Eta_Simple_updateGlobal(&Model);

#if (VISU)
		timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("interp+adv+visu timestep took: %.2f\n",toc);
#endif



	}

	printf("Simulation successfully completed\n");

//
//                          								END OF TIME LOOP													//
//																																//
//==============================================================================================================================//
//==============================================================================================================================//







	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //

	free(NonLin_x0);
	free(NonLin_dx);
	// Free memory
	printf("Free Physics->..\n");
	Physics_Memory_free(&Model);
	printf("Free NumStokes->..\n");
	Numbering_Memory_free(NumStokes);

	printf("Free EqStokes->..\n");
	EqSystem_Memory_free(EqStokes, SolverStokes);
	printf("Free BCStokes->..\n");
	BC_Memory_free(BCStokes);
#if (HEAT)
	printf("Free NumThermal->..\n");
	Numbering_Memory_free(NumThermal);
	printf("Free EqThermal->..\n");
	EqSystem_Memory_free(EqThermal,SolverThermal);
	printf("Free BCThermal->..\n");
	BC_Memory_free(BCThermal);
#endif
	printf("Free Particles->..\n");
	Particles_Memory_free(Particles, Grid);
	printf("Free Numerics->..\n");
	Numerics_Memory_free(Numerics);
	printf("Free Grid->..\n");
	Grid_Memory_free(Grid);
	printf("Free Output->..\n");
	Output_free(Output);


#if VISU
	// Quit glfw
	printf("Quit GFLW...\n");
	glfwDestroyWindow(Visu->window);
	glfwTerminate();
	printf("Free Visu->..\n");
	Visu_Memory_free(Visu);
#endif

printf("Memory freed successfully\n");

	//return EXIT_SUCCESS;
	return 1;

	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//

}




