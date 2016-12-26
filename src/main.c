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


int main(void) {

	printf("\n\n\n\n\n\n");
	printf("               ============================\n"
			"               ============================\n");
	printf("\n\n\n\n\n\nBeginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());

	int i;
#if (DARCY)
	int iy, ix, iCell;
#endif


	//exit(0);
	//int C = 0;
	//int ix, iy;
	// Declare structures
	// =================================
	// General
	Grid 		Grid;
	MatProps 	MatProps;
	Particles 	Particles;
	Physics 	Physics;
	Char 		Char;

	// Stokes: Conservation of momentum + continuity
	Numbering 	NumStokes;
	BC 			BCStokes;
	EqSystem 	EqStokes;
	Solver 		SolverStokes;


	// Heat conservation
	Numbering 	NumThermal;
	BC 			BCThermal;
	EqSystem  	EqThermal;
	Solver 		SolverThermal;

	// Numerics
	Numerics 	Numerics;


	// Visu
#if (VISU)
	Visu 		Visu;
#endif

	// Input
	Input 		Input;

	// Output
	Output 		Output;

	INIT_TIMER

	strcpy(Input.inputFile,"./Setups/Test/input.json");
	//strcpy(Input.inputFile,"/Users/abauville/JAMSTEC/StokesFD/Setups/Test/input.json");
	//strcpy(Input.inputFile,"/home/abauvill/mySoftwares/StokesFD/Setups/Test/input.json");

	printf("Reading input\n");
	Input_read(&Input, &Grid, &Numerics, &Physics, &MatProps, &Particles, &Char, &BCStokes, &BCThermal);

#if (LINEAR_VISCOUS)
	if (Numerics.maxNonLinearIter>1) {
		printf("error: you requested %i non linear iterations, however, they are switched off due to LINEAR_VISCOUS==true\n", Numerics.maxNonLinearIter);
		exit(0);
	}
#endif

#if (!HEAT)
	if (BCThermal.TB!=1.0 || BCThermal.TT!=1.0) {
		printf("TB = %.3e, TT = %.3e\n",BCThermal.TB, BCThermal.TT);
		printf("error: you specified non default thermal boundary conditions, however, the heat equation is switched off due to HEAT==false ");
		exit(0);
	}
#endif


#if (VISU)
	printf("Reading Visu input\n");
	Input_readVisu(&Input, &Visu);
#endif

	printf("Reading input over\n");


	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                       		 USER INPUT      		   	              	  //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	printf("nTimesteps = %i\n",Numerics.nTimeSteps);
	Grid.isPeriodic = false;
	if (BCStokes.SetupType==Stokes_SimpleShear) {
		Grid.isPeriodic = true;
		Grid.isFixed 	= true;
	}


	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);



	// Set characteristic quantities
	// =================================

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes, &BCThermal);

	//printf("Eta0[1] = %.3e", MatProps.eta0[1]);

	//Numerics.etaMin = 1E-5;
	//Numerics.etaMax = 1E3;
	Physics.epsRef = fabs(BCStokes.backStrainRate);

	printf("max backStrainRate = %.3e\n",BCStokes.backStrainRate);
	if (Physics.epsRef == 0)
		Physics.epsRef = 1E0;


	//BCThermal.SetupType = BCStokes.SetupType;

	Physics.maxV = (Grid.xmax-Grid.xmin)/Physics.epsRef;

	// Init grid and particles
	// =================================
	Grid.nCTot  = Grid.nxC*Grid.nyC;

	Grid.nxEC = Grid.nxC+2;
	Grid.nyEC = Grid.nyC+2;
	Grid.nECTot = Grid.nxEC*Grid.nyEC;


	Grid.nxVx 	= Grid.nxC+1; 		Grid.nyVx	= Grid.nyC+2;
	Grid.nxVy 	= Grid.nxC+2;		Grid.nyVy	= Grid.nyC+1;
	Grid.nxS 	= Grid.nxC+1;		Grid.nyS	= Grid.nyC+1;
	Grid.nSTot  = Grid.nxS*Grid.nyS;

	Grid.nVxTot = Grid.nxVx*Grid.nyVx;
	Grid.nVyTot = Grid.nxVy*Grid.nyVy;




	EqThermal.nEqIni 		= Grid.nECTot;






#if (DARCY)
	NumStokes.nSubEqSystem 	= 4;
	NumStokes.Stencil[0] 	= Stencil_Stokes_Darcy_Momentum_x; 	// Vx
	NumStokes.Stencil[1] 	= Stencil_Stokes_Darcy_Momentum_y; 	// Vy
	NumStokes.Stencil[2] 	= Stencil_Stokes_Darcy_Darcy;	   	// Pf
	NumStokes.Stencil[3] 	= Stencil_Stokes_Darcy_Continuity; 	// Pc
	EqStokes.nEqIni  	 	= Grid.nVxTot + Grid.nVyTot + Grid.nECTot + Grid.nECTot;
#else
	NumStokes.nSubEqSystem 	= 3;
	NumStokes.Stencil[0] 	= Stencil_Stokes_Momentum_x;		// Vx
	NumStokes.Stencil[1] 	= Stencil_Stokes_Momentum_y; 		// Vy
	NumStokes.Stencil[2]	= Stencil_Stokes_Continuity;		// P
	EqStokes.nEqIni  	 	= Grid.nVxTot + Grid.nVyTot + Grid.nECTot;
#endif









	NumThermal.nSubEqSystem 	= 1;
	NumThermal.Stencil[0] = Stencil_Heat;

	Particles.nPC 	= Particles.nPCX * Particles.nPCY;
	Particles.n 	= Grid.nCTot*Particles.nPC;

#if (VISU)
	Visu.ntri   	= Grid.nxC*Grid.nyC*2;//2;//Grid.nxC*Grid.nyC*2;
	Visu.ntrivert 	= Visu.ntri*3;
	Visu.nParticles = Particles.n+ (int) (Particles.n*0.1); // overallocate 5% of the number of particles
	//Visu.particleMeshRes = 6;
#endif

	printf("xmin = %.3f, ymin = %.3f\n", Grid.xmin, Grid.ymin);




	Grid.xmax_ini = Grid.xmax;
	Grid.xmin_ini = Grid.xmin;
	Grid.ymax_ini = Grid.ymax;
	Grid.ymin_ini = Grid.ymin;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;


	if (DEBUG) {
		if (Grid.nCTot>200) {
			printf("error: The system size exceeds the maximum allowed for debugging\n");
			exit(0);
		}
	}

	if (Grid.isPeriodic && Grid.nxC%2!=0) {
		printf("error: When using the periodic boundaries nxC must be even (because of node 'coloring'\n");
		exit(0);
	}


//======================================================================================================
//======================================================================================================
//
//                          				INITIALIZATION
//
	Numerics.timeStep = -1;
	Numerics.itNonLin = -1;

	// Other variables
	// =================================
	int iEq, iLS;

	//Init Grid
	// =================================
	Grid_allocateMemory(&Grid);
	Grid_init(&Grid, &Numerics);
	// Init Physics
	// =================================
	printf("Init Physics\n");
	Physics_allocateMemory	(&Physics, &Grid);

	// Init Numerics
	// =================================
	printf("Init Numerics\n");
	Numerics_init(&Numerics);

	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes			(&BCStokes , &Grid, &Physics, &EqStokes);
#if (HEAT)
	BC_initThermal			(&BCThermal, &Grid, &Physics, &EqThermal);
#endif

	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_allocateI		(&EqStokes);
	Numbering_allocateMemory(&NumStokes, &EqStokes, &Grid);
	Numbering_init			(&BCStokes, &Grid, &EqStokes, &NumStokes, &Physics);
	printf("EqSystem: init Stokes\n");
	EqSystem_allocateMemory	(&EqStokes );

	printf("Number of Unknowns for Stokes: %i \n", EqStokes.nEq);

#if (HEAT)
	printf("Numbering: init Thermal\n");
	EqSystem_allocateI		(&EqThermal);
	Numbering_allocateMemory(&NumThermal, &EqThermal, &Grid);
	Numbering_init			(&BCThermal, &Grid, &EqThermal, &NumThermal, &Physics);
	printf("EqSystem: init Thermal\n");
	EqSystem_allocateMemory	(&EqThermal);

	printf("Number of Unknowns for Heat: %i \n", EqThermal.nEq);
#endif

	// Initialize Particles
	// =================================
	printf("Particles: Init Particles\n");
	Particles_allocateMemory	(&Particles, &Grid);
	Particles_initCoord			(&Particles, &Grid);
	Particles_updateLinkedList	(&Particles, &Grid, &Physics); // in case a ridiculous amount of noise is put on the particle
	Particles_initPassive		(&Particles, &Grid);


	printf("In input assignphase\n");
	printf("Char.time = %.2e, Char.length = %.2e, Char.mass = %.2e\n", Char.time, Char.length, Char.mass);
	Input_assignPhaseToParticles(&Input, &Particles, &Grid, &Char);

#if (HEAT)
	// Get Init P to litho
	for (i = 0; i < Grid.nECTot; ++i) {
		//Physics.DT[i] = BCThermal.TB;//BCThermal.TB;
		Physics.DT[i] = BCThermal.TB;//BCThermal.TB;
	}
	Physics_interpTempFromCellsToParticle	(&Grid, &Particles, &Physics, &BCStokes,  &MatProps);
#endif
	Physics_getPhase					(&Physics, &Grid, &Particles, &MatProps, &BCStokes);
	Physics_interpFromParticlesToCell	(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);
	Physics_computeRho(&Physics, &Grid, &MatProps);
	Physics_initPToLithostatic 			(&Physics, &Grid);
	Physics_initEta(&Physics, &Grid, &MatProps);
#if (DEBUG)
	Physics_check(&Physics, &Grid, &Char);
#endif

	//Physics_computeEta					(&Physics, &Grid, &Numerics);
	// Init Solvers
	// =================================
	printf("EqStokes: Init Solver\n");
	EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes, false); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (&EqStokes, &SolverStokes);

#if (HEAT)
	printf("EqThermal: Init Solver\n");
	EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal, false); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (&EqThermal, &SolverThermal);
#endif
#if (VISU)
	// Init GLFW
	// =======================================

	Visu_allocateMemory(&Visu, &Grid);
	Visu_init(&Visu, &Grid, &Particles, &Char);
#endif


//
//                          				INITIALIZATION
//
//======================================================================================================
//======================================================================================================










//======================================================================================================
//======================================================================================================
//
//                     					INITIAL TEMPERATURE DISTRIBUTION
//

#if (HEAT)

	printf("EqThermal: compute the initial temperature distribution\n");



	Physics.dt = (3600*24*365.25 * 80E6)/Char.time; // initial value is really high to set the temperature profile. Before the advection, dt is recomputed to satisfy CFL
	Physics_computeRho(&Physics, &Grid, &MatProps);
	//Physics_computeThermalProps(&Physics, &Grid, &MatProps);
	EqSystem_assemble						(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal, false); // dummy assembly to give the EqSystem initSolvers
	//printf("P0 = %.2e\n", Physics.P[0]);
	EqSystem_solve							(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);


	// Add some random noise on the temperature
	srand(time(NULL));

	for (i = 0; i < EqThermal.nEq; ++i) {
		//EqThermal.x[i] += EqThermal.x[i]*(0.5 - (rand() % 1000)/1000.0)*0.2;
	}



	Physics_get_T_FromSolution				(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal, &Numerics);
	/*
	for (i = 0; i < Grid.nECTot; ++i) {
		Physics.DT[i] = Physics.T[i];
	}
	*/
	Physics_interpTempFromCellsToParticle	(&Grid, &Particles, &Physics, &BCStokes,  &MatProps);
	//Physics_interpFromParticlesToCell	 	(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);

	//Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
	Physics.dt = Numerics.dLmin/(3*min(MatProps.k,MatProps.nPhase)); // CFL condition, to get a reasonnable time step for the first computation of T



#endif

//
//                    					 INITIAL TEMPERATURE DISTRIBUTION
//
//======================================================================================================
//======================================================================================================









//======================================================================================================
//======================================================================================================
//
//                     					INITIAL DARCY STUFF
//

#if (DARCY)


	Physics_initPhi(&Physics, &Grid, &MatProps, &Numerics);
	Physics_interpPhiFromCellsToParticle	(&Grid, &Particles, &Physics);



#endif

//
//                    					 INITIAL DARCY STUFF
//
	//======================================================================================================
	//======================================================================================================


	// Update Cell Values with Part
	// =================================
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);
	Physics_computeRho(&Physics, &Grid, &MatProps);
	Physics_initPToLithostatic 			(&Physics, &Grid);

//Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);




	// Update BC
	// =================================
	printf("BC: Update\n");
	BCStokes.counter = 0;
	BC_updateStokes_Vel(&BCStokes, &Grid, &Physics, true);
#if (DARCY)
	BC_updateStokesDarcy_P(&BCStokes, &Grid, &Physics, true);
#endif




#if (DARCY)
		//Physics_computePhi(&Physics, &Grid, &Numerics, &BCStokes);
		Physics_computePerm(&Physics, &Grid, &Numerics, &MatProps);
#endif


		Physics_computeRho(&Physics, &Grid, &MatProps);
		//Physics_computePlitho(&Physics, &Grid);
		//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);




		// Initial viscosity
		// =======================================================




		TIC
		//EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes, true);
		TOC
		printf("Stokes Assembly: %.3f s\n", toc);

		// 								Assemble Stokes
		// ==========================================================================




//======================================================================================================
//======================================================================================================
//
//                          				TIME LOOP
//

	Numerics.timeStep = 0;
	Physics.dt = 1.0;//dtmax*1000;// pow(10,(log10(dtmin)+log10(dtmax))/2);
	Physics.time = 0;

	double timeStepTic;

	Physics.maxV = 1e2;
	while(Numerics.timeStep!=Numerics.nTimeSteps) {
		printf("\n\n\n          ========  Time step %i, t= %.2e Myrs  ========   \n"
				     "              ===================================== \n\n",Numerics.timeStep, Physics.time*Char.time/(3600*24*365*1e6));
#if (VISU)
		timeStepTic = glfwGetTime();
#endif
		Numerics.itNonLin = -1;
		// ==========================================================================
		// 							Solve the heat conservation

#if (HEAT)
		TIC
		printf("Heat assembly and solve\n");
		EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal, false);
		EqSystem_solve(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);
		Physics_get_T_FromSolution(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal, &Numerics);
		Physics_interpTempFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes, &MatProps);
		TOC
		printf("Temp Assembly+Solve+Interp: %.3f s\n", toc);



#endif

		// 							Solve the heat conservation
		// ==========================================================================





		// ==========================================================================
		// 								Assemble Stokes

		//Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
		//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);






		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										NON-LINEAR ITERATION											//
		Numerics.itNonLin = 0;
		EqStokes.normResidual = 1.0;
		Numerics.normRes0 = 1.0;
		Numerics.normResRef = 1.0;
		Numerics.cumCorrection_fac = 0.0;
		Numerics.lsLastRes = 1E15;
		Numerics.lsGlob = 1.00;
		Numerics.lsBestRes = 1e15;
		Numerics.lsBestGlob = 1.0;

		#if (!LINEAR_VISCOUS)
		compute* NonLin_x0 = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		compute* NonLin_dx = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		//compute* Sigma_xy0 = (compute*) malloc(Grid.nSTot * sizeof(compute));
		//compute* Sigma_xx0 = (compute*) malloc(Grid.nECTot * sizeof(compute));
		compute* EtaNonLin0 = (compute*) malloc(Grid.nECTot * sizeof(compute));
		compute* KhiNonLin0 = (compute*) malloc(Grid.nECTot * sizeof(compute));
#if (DARCY)
		compute* KhiBNonLin0 = (compute*) malloc(Grid.nECTot * sizeof(compute));
#endif

		//compute* EtaShearNonLin0 = (compute*) malloc(Grid.nSTot * sizeof(compute));
		//compute* KhiShearNonLin0 = (compute*) malloc(Grid.nSTot * sizeof(compute));

#endif


		//memcpy(Sigma_xx0, Physics.sigma_xx_0, Grid.nECTot * sizeof(compute));
		//memcpy(Sigma_xy0, Physics.sigma_xy_0, Grid.nSTot * sizeof(compute));

		Numerics.lsLastRes = 1E100;
		compute oldRes = EqStokes.normResidual;
		while((( (EqStokes.normResidual > Numerics.absoluteTolerance ) && Numerics.itNonLin<Numerics.maxNonLinearIter ) || Numerics.itNonLin<Numerics.minNonLinearIter)  || Numerics.cumCorrection_fac<=0.999) {
			printf("\n\n  ==== Non linear iteration %i ==== \n",Numerics.itNonLin);
			TIC
			Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
			TOC
			printf("update dt: %.3f s\n", toc);
/*
			memcpy(Sigma_xx0, Physics.sigma_xx_0, Grid.nECTot * sizeof(compute));
			memcpy(Sigma_xy0, Physics.sigma_xy_0, Grid.nSTot * sizeof(compute));
*/





			// =====================================================================================//
			//																						//
			// 										COMPUTE STOKES									//

			// update Dt
			//printf("####### before dt = %.2e\n", Physics.dt);

			//printf("####### dt = %.2e\n", Physics.dt);


			//Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
			//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);

			// Save X0
			//Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
			//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);
			memcpy(EtaNonLin0, Physics.eta, Grid.nECTot * sizeof(compute));
			memcpy(KhiNonLin0, Physics.khi, Grid.nECTot * sizeof(compute));
			//memcpy(EtaShearNonLin0, Physics.etaShear, Grid.nSTot * sizeof(compute));
			//memcpy(KhiShearNonLin0, Physics.khiShear, Grid.nSTot * sizeof(compute));
#if (DARCY)

			memcpy(KhiBNonLin0, Physics.khi_b, Grid.nECTot * sizeof(compute));
#endif

			memcpy(NonLin_x0, EqStokes.x, EqStokes.nEq * sizeof(compute));
			int i;
			/*
			for (i=0; i<EqStokes.nEq; ++i) {
				NonLin_x0[i] = EqStokes.x[i]*EqStokes.S[i];
			}
			*/





			//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);
			// Solve: A(X0) * X = b
			EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes, true);
			EqSystem_scale(&EqStokes);
			//EqSystem_check(&EqStokes);
			/*
			for (i=0; i<EqStokes.nEq; ++i) {
				NonLin_x0[i] /= EqStokes.S[i];
			}
			*/

			EqSystem_solve(&EqStokes, &SolverStokes, &Grid, &Physics, &BCStokes, &NumStokes);

			EqSystem_unscale(&EqStokes);

			// 										COMPUTE STOKES									//
			//																						//
			// =====================================================================================//


#if (LINEAR_VISCOUS)
			printf("/!\\ /!\\ LINEAR_VISCOUS==true, Non-linear iterations are ineffective/!\\ \n");
			Physics_get_VxVy_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_get_P_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes, &Numerics);
			Physics_computeRho(&Physics, &Grid);
			Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);
			break;
#else
			// =====================================================================================//
			//																						//
			// 										LINE SEARCH										//
			// Compute dx
			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes.x[iEq] - NonLin_x0[iEq];
			}


			Numerics.minRes = 1E100;
			Numerics.lsGlob = 1.0;
			Numerics.lsState = -1;
			iLS = 0;
			compute oldRes = EqStokes.normResidual;

			while (iLS < Numerics.nLineSearch+1) {
				//printf("== Line search %i:  ", iLS);

				// X1 = X0 + a*(X-X0)
				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					EqStokes.x[iEq] = NonLin_x0[iEq] + Numerics.lsGlob*(NonLin_dx[iEq]);
				}



				/*
				for (i=0;i<Grid.nSTot;++i) {
					Physics.etaShear[i] = EtaShearNonLin0[i] ;
					Physics.khiShear[i] = KhiShearNonLin0[i] ;
				}
				*/




				// Update the stiffness matrix
				Physics_get_VxVy_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_get_P_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes, &Numerics);


#if (DARCY)
				Physics_computePhi(&Physics, &Grid, &Numerics);
				Physics_computePerm(&Physics, &Grid, &Numerics, &MatProps);
#endif

#if (DEBUG)
				printf("before computeEta\n");
				//Physics_check(&Physics, &Grid, &Char);
#endif

				TIC
				Physics_computeRho(&Physics, &Grid, &MatProps);
				Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);
				TOC
				printf("Compute Rho, Stress, Eta: %.3f s\n", toc);
#if (DEBUG)
				printf("after computeEta\n");
				Physics_check(&Physics, &Grid, &Char);
#endif

				TIC
				EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes, false);
				TOC
				printf("Assembly: %.3f s\n", toc);

				// compute the norm of the  residual:
				// F = b - A(X1) * X1
				TIC
				EqSystem_computeNormResidual(&EqStokes);
				TOC
				printf("Compute Res: %.3f s\n", toc);
				// update the best globalization factor and break if needed
				//int Break = Numerics_updateBestGlob(&Numerics, &EqStokes, &iLS);
				//Numerics_LineSearch_chooseGlob(&Numerics, &EqStokes);
				/*
				if (Numerics.lsState < 0) {
					//printf("Break!!\n");
					break;
				}
				*/

				printf("a = %.3f,  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", Numerics.lsGlob, fabs(EqStokes.normResidual-oldRes), EqStokes.normResidual);

				if (EqStokes.normResidual<Numerics.minRes) {
					Numerics.minRes = EqStokes.normResidual;
					Numerics.lsBestGlob = Numerics.lsGlob;
				}
				iLS++;
				if (iLS<Numerics.nLineSearch) {
					Numerics.lsGlob = Numerics.lsGlob/2.0;
				} else {
					if (Numerics.lsGlob == Numerics.lsBestGlob) {
						break;
					} else {
						Numerics.lsGlob = Numerics.lsBestGlob;
					}
				}

				//printf("minRes = %.2e, lastRes = %.2e\n",Numerics.minRes, Numerics.lsLastRes);
				if (Numerics.minRes<Numerics.lsLastRes) {
					break;
				}


			}

			Numerics.cumCorrection_fac += Numerics.lsBestGlob;
			Numerics.lsLastRes = EqStokes.normResidual;

			if (Numerics.lsState == -2) {
				//printf("Break!!\n");
				break;
			}

			if (fabs(EqStokes.normResidual-oldRes)<EqStokes.normResidual/500.0) {
				break;
			}



#if NON_LINEAR_VISU


		Visu.update = true;
		Visu.updateGrid = false;
		Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps, &EqStokes, &EqThermal, &NumStokes, &NumThermal);
		if (glfwWindowShouldClose(Visu.window))
			break;
#endif







			// 		   								LINE SEARCH										//
			//																						//
			// =====================================================================================//
#endif



			Numerics.itNonLin++;
		} // end of non-linear loop




#if (!LINEAR_VISCOUS)
		free(EtaNonLin0);
		free(KhiNonLin0);
		//free(EtaShearNonLin0);
		//free(KhiShearNonLin0);
#if (DARCY)
		free(KhiBNonLin0);
#endif
		free(NonLin_x0);
		free(NonLin_dx);
#endif
#if (VISU)
		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);
#endif

		if (isnan(EqStokes.normResidual)) {
			printf("\n\n\n\nerror: Something went wrong. The norm of the residual is NaN\n");
			break;
		}



		// 										NON-LINEAR ITERATION 											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//










#if VISU


		Visu.update = true;
		if (~Grid.isFixed) {
			Visu.updateGrid = true;
		}
		Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps, &EqStokes, &EqThermal, &NumStokes, &NumThermal);
		if (glfwWindowShouldClose(Visu.window))
			break;
#endif


		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										ADVECTION AND INTERPOLATION										//

#if(HEAT)
		printf("####### bef end dt = %.2e\n", Physics.dt);
		if (Numerics.maxNonLinearIter==1) {
			Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
		}
		printf("####### end dt = %.2e\n", Physics.dt);
#endif


		//Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
		//printf("dt = %.3e, dtAdv = %.3e\n",Physics.dt,Physics.dtAdv);
		// update stress on the particles
		// =============================
		Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
		//Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);



		/*
		for (i = 0; i < 100; ++i) {
				Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes, &MatProps);
		}
		*/

		Physics_interpStressesFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal, &MatProps);

#if (DARCY)
		Physics_interpPhiFromCellsToParticle	(&Grid, &Particles, &Physics);


		compute dx, dy;
		for (iy = 1; iy < Grid.nyEC-1; ++iy) {
			for (ix = 1; ix < Grid.nxEC-1; ++ix) {
				iCell = ix + iy*Grid.nxEC;
				dx = Grid.DXS[ix-1];
				dy = Grid.DYS[iy-1];
				Physics.divV0[iCell]  = (  Physics.Vx[ix+iy*Grid.nxVx] - Physics.Vx[ix-1+ iy   *Grid.nxVx]  )/dx;
				Physics.divV0[iCell] += (  Physics.Vy[ix+iy*Grid.nxVy] - Physics.Vy[ix  +(iy-1)*Grid.nxVy]  )/dy;
			}
		}
#endif






		// Output
		// =================
		/*
		printf("Write output ...\n");
		Output_modelState(&Output, &Grid, &Physics, &Char, &Numerics);
		printf("Success1...\n");
		Output_data(&Output, &Grid, &Physics, &Char, &Numerics);
		printf("Success2!!!\n");
		*/








		// Advect Particles
		// =============================
		printf("Particles: Advect\n");
		Particles_advect(&Particles, &Grid, &Physics);

		// Advect the box and update Particles position if needed
		// =============================
		switch (BCStokes.SetupType) {
		case Stokes_PureShear:
		case Stokes_Sandbox:
			if (Grid.isFixed) {
				Particles_deleteIfOutsideTheDomain(&Particles, &Grid);
				Particles_injectAtTheBoundaries(&Particles, &Grid);
			} else {
				Grid_updatePureShear(&Grid, &BCStokes, &Numerics, Physics.dt);
				Particles_teleportInsideTheDomain(&Particles, &Grid, &Physics);
			}
			break;
		case Stokes_SimpleShear:
			Particles_Periodicize(&Particles, &Grid);
			break;
		case Stokes_FixedLeftWall:
			break;
		case Stokes_CornerFlow:
			if (Grid.isFixed) {
				Particles_deleteIfOutsideTheDomain(&Particles, &Grid);
				Particles_injectAtTheBoundaries(&Particles, &Grid);
			} else {
				printf("error: For the corner flow setup, Grid.fixedBox must be true. Correct the input file");
				exit(0);
			}
			break;
		default:
			break;
		}


		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(&Particles, &Grid, &Physics);
		Particles_injectOrDelete(&Particles, &Grid);

		// Update the Phase matrix
		// =================================
		Physics_getPhase					(&Physics, &Grid, &Particles, &MatProps, &BCStokes);

		// Update the Physics on the Cells
		// =================================
		printf("Physics: Interp from particles to cell\n");
		Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);


		// Update BC
		// =================================
		printf("BC: Update\n");
		BCStokes.counter = 0;
		BC_updateStokes_Vel(&BCStokes, &Grid, &Physics, true);
#if (DARCY)
		BC_updateStokesDarcy_P(&BCStokes, &Grid, &Physics, true);
#endif
#if (HEAT)
		BCThermal.counter = 0;
		BC_updateThermal(&BCThermal, &Grid, &Physics, true);
#endif


		// 										ADVECTION AND INTERPOLATION 									//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//





		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);
		Physics.time += Physics.dt;

		Numerics.timeStep++;










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

	// Free memory
	printf("Free Physics...\n");
	Physics_freeMemory(&Physics, &Grid);
	printf("Free NumStokes...\n");
	Numbering_freeMemory(&NumStokes);

	printf("Free EqStokes...\n");
	EqSystem_freeMemory(&EqStokes, &SolverStokes);
	printf("Free BCStokes...\n");
	BC_freeMemory(&BCStokes);
#if (HEAT)
	printf("Free NumThermal...\n");
	Numbering_freeMemory(&NumThermal);
	printf("Free EqThermal...\n");
	EqSystem_freeMemory(&EqThermal,&SolverThermal);
	printf("Free BCThermal...\n");
	BC_freeMemory(&BCThermal);
#endif
	printf("Free Particles...\n");
	Particles_freeMemory(&Particles, &Grid);
	printf("Free Numerics...\n");
	Numerics_freeMemory(&Numerics);
	printf("Free Grid...\n");
	Grid_freeMemory(&Grid);


#if VISU
	// Quit glfw
	printf("Quit GFLW...\n");
	glfwDestroyWindow(Visu.window);
	glfwTerminate();
	printf("Free Visu...\n");
	Visu_freeMemory(&Visu);
#endif

printf("Memory freed successfully\n");

	return EXIT_SUCCESS;

	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//

}




