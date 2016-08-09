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

	//strcpy(Input.inputFile,"./Setups/Test/input.json");
	strcpy(Input.inputFile,"/Users/abauville/JAMSTEC/StokesFD/Setups/Test/input.json");
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





	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);



	// Set characteristic quantities
	// =================================

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes, &BCThermal);

	//printf("Eta0[1] = %.3e", MatProps.eta0[1]);

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
	//Numerics.dtmax = 1.0;
	//Numerics.dtmin = 1E-100;

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

	printf("Number of Velocity - Pressure unknowns: %i \n", Grid.nVxTot + Grid.nVyTot + Grid.nCTot);


	Grid.xmax_ini = Grid.xmax;
	Grid.xmin_ini = Grid.xmin;
	Grid.ymax_ini = Grid.ymax;
	Grid.ymin_ini = Grid.ymin;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	printf("###########nxC = %i\n",Grid.nxC);

	if (DEBUG) {
		if (Grid.nCTot>200) {
			printf("error: The system size exceeds the maximum allowed for debugging\n");
			exit(0);
		}
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
	Grid_init(&Grid, &Input, &Numerics);
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
	Numbering_init			(&BCStokes, &Grid, &EqStokes, &NumStokes, &Physics);
	printf("EqSystem: init Stokes\n");
	EqSystem_allocateMemory	(&EqStokes );

#if (HEAT)
	printf("Numbering: init Thermal\n");
	EqSystem_allocateI		(&EqThermal);
	Numbering_allocateMemory(&NumThermal, &EqThermal, &Grid);
	Numbering_init			(&BCThermal, &Grid, &EqThermal, &NumThermal, &Physics);
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
	printf("Char.time = %.2e, Char.length = %.2e, Char.mass = %.2e\n", Char.time, Char.length, Char.mass);
	Input_assignPhaseToParticles(&Input, &Particles, &Grid, &Char);


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
	Physics.dtT = (3600*24*365.25 * 100E6)/Char.time; // initial value is really high to set the temperature profile. Before the advection, dt is recomputed to satisfy CFL
	Physics_computeRho(&Physics, &Grid);
	EqSystem_assemble						(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal); // dummy assembly to give the EqSystem initSolvers
	printf("P0 = %.2e\n", Physics.P[0]);
	EqSystem_solve							(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);


	// Add some random noise on the temperature
	srand(time(NULL));
	for (i = 0; i < EqThermal.nEq; ++i) {
		EqThermal.x[i] += EqThermal.x[i]*(0.5 - (rand() % 1000)/1000.0)*0.1;
	}



	Physics_get_T_FromSolution				(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);

	Physics_interpTempFromCellsToParticle	(&Grid, &Particles, &Physics, &BCStokes,  &MatProps);
	//Physics_interpFromParticlesToCell	 	(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);

	Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
	Physics.dtT = Numerics.dLmin/(3*min(MatProps.k,MatProps.nPhase)); // CFL condition, to get a reasonnable time step for the first computation of T



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


	Physics_initPhi(&Physics, &Grid);


	Physics_interpPhiFromCellsToParticle	(&Grid, &Particles, &Physics);



#endif

//
//                    					 INITIAL DARCY STUFF
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
		Physics_interpTempFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes, &MatProps);
		TOC
		printf("Temp Assembly+Solve+Interp: %.2fs\n", toc);



#endif

		// 							Solve the heat conservation
		// ==========================================================================




		// ==========================================================================
		// 								Assemble Stokes


#if (DARCY)
		memcpy(Physics.phi0,Physics.phi, Grid.nECTot*sizeof(compute));
		memcpy(Physics.Pc0,Physics.Pc, Grid.nECTot*sizeof(compute));
		printf("***********phi = %.2e\n",Physics.phi[150]);

		Physics_computePhi(&Physics, &Grid, &Numerics, &BCStokes);
		Physics_computePerm(&Physics, &Grid, &Numerics, &BCStokes);
		Physics_computeRho(&Physics, &Grid);
		printf("***********phi = %.2e\n",Physics.phi[150]);
#endif

		Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes);

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

/*
#if VISU

			// Update only if user input are received
			//Visu.paused = true;
			//Visu.update = true;

			//Visu.update = false;
			Visu.updateGrid = false;
			Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps);
			if (glfwWindowShouldClose(Visu.window))
				break;

#endif
*/


			// =====================================================================================//
			//																						//
			// 										COMPUTE STOKES									//

			// update Dt
			Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);

			// Save X0
			memcpy(NonLin_x0, EqStokes.x, EqStokes.nEq * sizeof(compute));

			// Solve: A(X0) * X = b
			EqSystem_solve(&EqStokes, &SolverStokes, &Grid, &Physics, &BCStokes, &NumStokes);


			// 										COMPUTE STOKES									//
			//																						//
			// =====================================================================================//








#if (LINEAR_VISCOUS)
			printf("/!\\ /!\\ LINEAR_VISCOUS==true, Non-linear iterations are ineffective/!\\ \n");
			Physics_get_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
			Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
			Physics_computeEta(&Physics, &Grid, &Numerics);
#else
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
				//printf("== Line search %i:  ", iLS);

				// X1 = X0 + a*(X-X0)
				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					EqStokes.x[iEq] = NonLin_x0[iEq] + Numerics.glob[iLS]*(NonLin_dx[iEq]);
				}

				// Update the stiffness matrix
				Physics_get_VxVy_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_get_P_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
#if (DARCY)
				Physics_computePhi(&Physics, &Grid, &Numerics, &BCStokes);
				Physics_computePerm(&Physics, &Grid, &Numerics, &BCStokes);
				Physics_computeRho(&Physics, &Grid);
#endif
				Physics_computeEta(&Physics, &Grid, &Numerics, &BCStokes);

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
#endif









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
				printf("/!\\ /!\\  Warning  /!\\ /!\\ : The residual is larger than the tolerance. The non linear iterations might be diverging. Wiping up the solution and starting the iteration again\n");
			}
			*/

			// 						Blowing up check: if the residual is too large					//
			// =====================================================================================//


			Numerics.itNonLin++;
		} // end of non-linear loop


#if (!LINEAR_VISCOUS)
		free(NonLin_x0);
		free(NonLin_dx);
#endif
#if (VISU)
		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);
#endif


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
			Physics_updateDt(&Physics, &Grid, &MatProps, &Numerics);
		}
#endif
		printf("dt = %.3e, dtAdv = %.3e\n",Physics.dt,Physics.dtAdv);
		// update stress on the particles
		// =============================
		Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
		Physics_interpStressesFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal, &MatProps);

#if (DARCY)
		Physics_interpPhiFromCellsToParticle	(&Grid, &Particles, &Physics);

		int iy, ix, iCell;
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
			Particles_Periodicize(&Particles, &Grid);
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
		Particles_injectOrDelete(&Particles, &Grid);

		// Update the Physics on the Cells
		// =================================
		printf("Physics: Interp from particles to cell\n");
		Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);



		// Update BC
		// =================================
		printf("BC: Update\n");
		BCStokes.counter = 0;
		BC_updateStokes_Vel(&BCStokes, &Grid, true);
#if (HEAT)
		BCThermal.counter = 0;
		BC_updateThermal(&BCThermal, &Grid, true);
#endif


		// 										ADVECTION AND INTERPOLATION 									//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//





		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);
		Physics.time += Physics.dt;

		Numerics.timeStep++;
#if VISU
		/*
		//int iy, ix, iCell;
		for (iy = 0; iy < Grid.nyEC; ++iy) {
			for (ix = 0; ix < Grid.nxEC; ++ix) {
				iCell = ix+iy*Grid.nxEC;
				if (iy==5) {
					printf("Physics->phi [iCell] = %.2e\n",Physics.phi [iCell]);
				}

			}

		}
		*/

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
	printf("Free Grid...\n");
	Grid_freeMemory(&Grid);
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

