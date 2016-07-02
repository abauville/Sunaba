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

	INIT_TIMER
	strcpy(Input.inputFile,"/Users/abauville/Work/StokesFD/Setups/Test/input.json");
	//strcpy(Input.inputFile,"/Users/abauville/JAMSTEC/StokesFD/Setups/Test/input.json");
	//strcpy(Input.inputFile,"/home/abauvill/mySoftwares/StokesFD/Setups/Test/input.json");

	printf("Reading input\n");
	Input_read(&Input, &Grid, &Numerics, &Physics, &MatProps, &Particles, &Char, &BCStokes, &BCThermal);




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

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;



	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);



	// Set characteristic quantities
	// =================================

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes, &BCThermal);

	printf("Eta0[1] = %.3e", MatProps.eta0[1]);

	Numerics.etaMin = 1E-5;
	Numerics.etaMax = 1E3;
	Physics.epsRef = fabs(BCStokes.backStrainRate);

	printf("max backStrainRate = %.3e\n",BCStokes.backStrainRate);
	if (Physics.epsRef == 0)
		Physics.epsRef = 1E0;




	// Determine dtmin and dt max based on maxwell times on materials

	/*
	for (i=0;i<MatProps.nPhase;++i) {
		if (MatProps.maxwellTime[i]<dtmin) {
			dtmin = MatProps.maxwellTime[i];
		}
		if (MatProps.maxwellTime[i]>dtmax) {
			dtmax = MatProps.maxwellTime[i];
		}
		printf("maxwell time = %.2e\n", MatProps.maxwellTime[i]);
	}
	*/
	Numerics.dtmax = 1.0;
	Numerics.dtmin = 1E-100;

	BCThermal.SetupType = BCStokes.SetupType;

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

	EqStokes.nEqIni  		= Grid.nxVx*Grid.nyVx + Grid.nxVy*Grid.nyVy + Grid.nxC*Grid.nyC;
	EqThermal.nEqIni 		= (Grid.nxC+2)*(Grid.nyC+2);
	NumStokes.nSubEqSystem 	= 3;
	NumStokes.Stencil[0] = Vx;
	NumStokes.Stencil[1] = Vy;
	NumStokes.Stencil[2] = P;

	NumThermal.nSubEqSystem 	= 1;
	NumThermal.Stencil[0] = T;

	Particles.nPC 	= Particles.nPCX * Particles.nPCY;
	Particles.n 	= Grid.nCTot*Particles.nPC;

#if (VISU)
	Visu.ntri   	= 2;//Grid.nxC*Grid.nyC*2;
	Visu.ntrivert 	= Visu.ntri*3;
	Visu.nParticles = Particles.n+ (int) (Particles.n*0.1); // overallocate 5% of the number of particles
	//Visu.particleMeshRes = 6;
#endif

	printf("xmin = %.3f, ymin = %.3f\n", Grid.xmin, Grid.ymin);


	Grid.xmax_ini = Grid.xmax;
	Grid.xmin_ini = Grid.xmin;
	Grid.ymax_ini = Grid.ymax;
	Grid.ymin_ini = Grid.ymin;



	Numerics.dLmin = fmin(Grid.dx,Grid.dy);



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

	// Init Physics
	// =================================
	printf("Init Physics\n");
	Physics_allocateMemory	(&Physics, &Grid);

	// Init Numerics
	// =================================
	printf("Init Physics\n");
	Numerics_init(&Numerics);

	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes			(&BCStokes , &Grid, &EqStokes);
#if (HEAT)
	BC_initThermal			(&BCThermal, &Grid, &EqThermal);
#endif


	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_allocateI		(&EqStokes);
	Numbering_allocateMemory(&NumStokes, &EqStokes, &Grid);
	Numbering_init			(&BCStokes, &Grid, &EqStokes, &NumStokes);
	printf("EqSystem: init Stokes\n");
	EqSystem_allocateMemory	(&EqStokes );


#if (HEAT)
	printf("Numbering: init Thermal\n");
	EqSystem_allocateI		(&EqThermal);
	Numbering_allocateMemory(&NumThermal, &EqThermal, &Grid);
	Numbering_init			(&BCThermal, &Grid, &EqThermal, &NumThermal);
	printf("EqSystem: init Thermal\n");
	EqSystem_allocateMemory	(&EqThermal);
#endif


	// Initialize Particles
	// =================================
	printf("Particles: Init Particles\n");
	Particles_allocateMemory	(&Particles, &Grid);
	Particles_initCoord			(&Particles, &Grid);
	Particles_updateLinkedList	(&Particles, &Grid, &Physics); // in case a ridiculous amount of noise is put on the particle
	Particles_initPassive		(&Particles, &Grid);

	printf("In input assignphase\n");
	Input_assignPhaseToParticles(&Input, &Particles, &Grid);


	// Get Init P to litho
	Physics_interpFromParticlesToCell	(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);
	Physics_initPToLithostatic 			(&Physics, &Grid);

	//Physics_computeEta					(&Physics, &Grid, &Numerics);
	// Init Solvers
	// =================================
	printf("EqStokes: Init Solver\n");
	EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (&EqStokes, &SolverStokes);


#if (HEAT)
	printf("EqThermal: Init Solver\n");
	EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (&EqThermal, &SolverThermal);
#endif

#if (VISU)
	// Init GLFW
	// =======================================

	Visu_allocateMemory(&Visu, &Grid);
	Visu_init(&Visu, &Grid, &Particles);
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
	Physics.dt = (3600*24*365.25 * 100E6)/Char.time; // initial value is really high to set the temperature profile. Before the advection, dt is recomputed to satisfy CFL

	EqSystem_assemble						(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal); // dummy assembly to give the EqSystem initSolvers
	printf("P0 = %.2e\n", Physics.P[0]);
	EqSystem_solve							(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);


	/*
	// Add some random noise on the temperature
	int i;
	srand(time(NULL));
	for (i = 0; i < EqThermal.nEq; ++i) {
		EqThermal.x[i] += EqThermal.x[i]*(0.5 - (rand() % 1000)/1000.0)*0.1;
	}
	*/



	Physics_get_T_FromSolution				(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);

	Physics_interpTempFromCellsToParticle	(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal);
	//Physics_interpFromParticlesToCell	 	(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);





#endif

//
//                    					 INITIAL TEMPERATURE DISTRIBUTION
//
//======================================================================================================
//======================================================================================================

















//======================================================================================================
//======================================================================================================
//
//                          				TIME LOOP
//

	Numerics.timeStep = 0;
	Physics.dt = 1.0;//dtmax*1000;// pow(10,(log10(dtmin)+log10(dtmax))/2);
	Physics.time = 0;
	Numerics.itNonLin = -1;
	double timeStepTic;

	while(Numerics.timeStep!=Numerics.nTimeSteps) {
		printf("\n\n\n          ========  Time step %i  ========   \n"
				     "       ===================================== \n\n",Numerics.timeStep);
		timeStepTic = glfwGetTime();

		// ==========================================================================
		// 							Solve the heat conservation

#if (HEAT)
		TIC
		printf("Heat assembly and solve\n");
		EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal);
		EqSystem_solve(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);
		Physics_get_T_FromSolution(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);
		Physics_interpTempFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal);
		TOC
		printf("Temp Assembly+Solve+Interp: %.2fs\n", toc);



#endif

		// 							Solve the heat conservation
		// ==========================================================================




		// ==========================================================================
		// 								Assemble Stokes

		Physics_computeEta(&Physics, &Grid, &Numerics);
		TIC
		EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);
		TOC
		printf("Stokes Assembly: %.2fs\n", toc);

		// 								Assemble Stokes
		// ==========================================================================










		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										NON-LINEAR ITERATION											//
		Numerics.itNonLin = 0;
		EqStokes.normResidual = 1.0;
		Numerics.normRes0 = 1.0;
		Numerics.normResRef = 1.0;
		compute* NonLin_x0 = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		compute* NonLin_dx = (compute*) malloc(EqStokes.nEq * sizeof(compute));


		while(( (EqStokes.normResidual > Numerics.absoluteTolerance ) && Numerics.itNonLin!=Numerics.maxNonLinearIter ) || Numerics.itNonLin<Numerics.minNonLinearIter) {
			printf("\n\n  ==== Non linear iteration %i ==== \n",Numerics.itNonLin);


#if VISU

			// Update only if user input are received
			Visu.update = false;
			Visu.updateGrid = false;
			Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps);
			if (glfwWindowShouldClose(Visu.window))
				break;

#endif


			// =====================================================================================//
			//																						//
			// 										COMPUTE STOKES									//

			// update Dt
			Physics_updateDt(&Physics, &Numerics);

			// Save X0
			memcpy(NonLin_x0, EqStokes.x, EqStokes.nEq * sizeof(compute));

			// Solve: A(X0) * X = b
			EqSystem_solve(&EqStokes, &SolverStokes, &Grid, &Physics, &BCStokes, &NumStokes);


			// 										COMPUTE STOKES									//
			//																						//
			// =====================================================================================//







			// =====================================================================================//
			//																						//
			// 										LINE SEARCH										//

			// Compute dx
			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes.x[iEq] - NonLin_x0[iEq];
			}

			Numerics.glob[Numerics.nLineSearch] = 1.0/Numerics.nLineSearch; // this is the best value
			Numerics.minRes = 1E100;
			for (iLS = 0; iLS < Numerics.nLineSearch+1; ++iLS) {
				printf("== Line search %i:  ", iLS);

				// X1 = X0 + a*(X-X0)
				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					EqStokes.x[iEq] = NonLin_x0[iEq] + Numerics.glob[iLS]*(NonLin_dx[iEq]);
				}

				// Update the stiffness matrix
				Physics_get_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_computeEta(&Physics, &Grid, &Numerics);
				EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);

				// compute the norm of the  residual:
				// F = b - A(X1) * X1
				EqSystem_computeNormResidual(&EqStokes);

				// update the best globalization factor and break if needed
				int Break = Numerics_updateBestGlob(&Numerics, &EqStokes, iLS);
				if (Break==1) {
					break;
				}
			}

			// 		   								LINE SEARCH										//
			//																						//
			// =====================================================================================//







			// =====================================================================================//
			// 				     	Blowing up check: if the residual is too large					//

			/*
			// wipe up the solution vector and start the iteration again with 0 everywhere initial guess
			if (Numerics.timeStep>1 && Numerics.minRes>10.0) {
				for (i=0; i<EqStokes.nEq; i++) {
					EqStokes.x[i] = 0;
				}
				Physics_get_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Numerics.itNonLin = 0;
				printf("/! /!  Warning  /! /! : The residual is larger than the tolerance. The non linear iterations might be diverging. Wiping up the solution and starting the iteration again\n");
			}
			*/

			// 						Blowing up check: if the residual is too large					//
			// =====================================================================================//


			Numerics.itNonLin++;
		} // end of non-linear loop



		free(NonLin_x0);
		free(NonLin_dx);
		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);


		// 										NON-LINEAR ITERATION 											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//












		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										ADVECTION AND INTERPOLATION										//

#if(HEAT)
		if (Numerics.maxNonLinearIter==1) {
			Physics_updateDt(&Physics, &Numerics);
		}
#endif

		// update stress on the particles
		// =============================
		Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
		Physics_interpStressesFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal, &MatProps);



		// Advect Particles
		// =============================
		printf("Particles: Advect\n");
		Particles_advect(&Particles, &Grid, &Physics);

		// Advect the box and update Particles position if needed
		// =============================
		switch (BCStokes.SetupType) {
		case PureShear:
		case Sandbox:
			Grid_updatePureShear(&Grid, &BCStokes, Physics.dt);
			Particles_teleportInsideTheDomain(&Particles, &Grid, &Physics);
			break;
		case SimpleShearPeriodic:
			Particles_Periodicize(&Particles, &Grid, &BCStokes);
			break;
		case FixedLeftWall:
			break;

		default:
			break;
		}


		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(&Particles, &Grid, &Physics);


		// Update the Physics on the Cells
		// =================================
		printf("Physics: Interp from particles to cell\n");
		Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);



		// Update BC
		// =================================
		printf("BC: Update\n");
		BC_updateStokes(&BCStokes, &Grid);
#if (HEAT)
		BC_updateThermal(&BCThermal, &Grid);
#endif


		// 										ADVECTION AND INTERPOLATION 									//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//





		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);
		Physics.time += Physics.dt;

		Numerics.timeStep++;
#if VISU
		Visu.update = true;
		Visu.updateGrid = true;
		Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps);
		if (glfwWindowShouldClose(Visu.window))
			break;
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

	// Free memory
	printf("Free Physics...\n");
	Physics_freeMemory(&Physics);
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
	printf("Memory freed successfully\n");



#if VISU
	// Quit glfw
	glfwDestroyWindow(Visu.window);
	glfwTerminate();
	Visu_freeMemory(&Visu);
#endif

	printf("SUCCESS!\n\n\n");

	return EXIT_SUCCESS;

	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//

}

