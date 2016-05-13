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

	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                               INITIALIZATION                               //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	printf("\n\n\n\n\n\n");
	printf("               ============================\n"
		   "               ============================\n");
	printf("\n\n\n\n\n\nBeginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());

	//int C = 0;
	//int ix, iy;
	int i;
	// Declare structures
	// =================================
	// General
	Grid 		Grid;
	MatProps 	MatProps;
	Particles 	Particles;
	Physics 	Physics;
	Visu 		Visu;
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

	// Darcy
	Darcy Darcy;




	INIT_TIMER

	// Set model properties
	// =================================
	int nTimeSteps  = 1; //  negative value for infinite
	int nLineSearch = 3;
	int maxNonLinearIter = 10; // should always be greater than the number of line searches
	int minNonLinearIter = 5; // should always be greater than the number of line searches
	compute relativeTolerance = 3E-5; // relative tolerance to the one of this time step
	compute absoluteTolerance = 3E-5; // relative tolerance to the first one of the simulation
	compute maxCorrection = 1.0;

	Grid.nxC = 256;
	Grid.nyC = 128;

	Particles.nPCX = 5;
	Particles.nPCY = 5;

	//Grid.xmin = 0;
	//Grid.xmax = (compute) Grid.nxC;
	//Grid.ymin = 0;
	//Grid.ymax = (compute) Grid.nyC;
	Grid.xmin =   0*50E3;
	Grid.xmax =   6*50E3;
	Grid.ymin =   0*50E3;
	Grid.ymax =   1*50E3;

	MatProps.nPhase  = 4;

	//MatProps.rho0[0] = 1; 		MatProps.eta0[0] = 1.0;  		MatProps.n[0] = 1.0; 		MatProps.flowLaw[0] = PowerLawViscous;
	//MatProps.rho0[1] = 1;		MatProps.eta0[1] = 0.001; 		MatProps.n[1] = 1.0;		MatProps.flowLaw[1] = PowerLawViscous;

	MatProps.rho0[0] = 10; 			MatProps.eta0[0] = 1E17;  		MatProps.n[0] = 1.0; 		MatProps.flowLaw[0] = PowerLawViscous;
	MatProps.rho0[1] = 2700;		MatProps.eta0[1] = 1E23; 		MatProps.n[1] = 1.0;		MatProps.flowLaw[1] = PowerLawViscous;
	MatProps.rho0[2] = 2700;		MatProps.eta0[2] = 1E20; 		MatProps.n[2] = 1.0;		MatProps.flowLaw[2] = PowerLawViscous;
	MatProps.rho0[3] = 2700;		MatProps.eta0[3] = 1E23; 		MatProps.n[3] = 1.0;		MatProps.flowLaw[3] = PowerLawViscous;

	MatProps.alpha[0] = 1E-5;  	MatProps.beta [0] = 0.0;  		MatProps.k[0] = 1E-2; 			MatProps.G[0] = 1E11;
	MatProps.alpha[1] = 1E-5; 	MatProps.beta [1] = 0.0;  		MatProps.k[1] = 1E-2; 			MatProps.G[1] = 1E11;
	MatProps.alpha[2] = 1E-5; 	MatProps.beta [2] = 0.0;  		MatProps.k[2] = 1E-2; 			MatProps.G[2] = 1E11;
	MatProps.alpha[3] = 1E-5; 	MatProps.beta [3] = 0.0;  		MatProps.k[3] = 1E-2; 			MatProps.G[3] = 1E11;

	MatProps.cohesion[0] = 10000.0*1E6; 	MatProps.frictionAngle[0] = 30*PI/180; //air
	MatProps.cohesion[1] = 10.0*1E6;		MatProps.frictionAngle[1] = 30*PI/180; // green
	MatProps.cohesion[2] = 10.0*1E6;		MatProps.frictionAngle[2] = 5*PI/180; // orange
	MatProps.cohesion[3] = 100.0*1E6;		MatProps.frictionAngle[3] = 30*PI/180; // blue


	MatProps.SD[0] = 1E-3; // 1/m
	MatProps.SD[1] = 1E-3; // 1/m
	MatProps.SD[2] = 1E-3; // 1/m
	MatProps.SD[3] = 1E-3; // 1/m

	MatProps.kD[0] = 1.0E-7; // m/s
	MatProps.kD[1] = 1.0E-7; // m/s
	MatProps.kD[2] = 1.0E-7; // m/s
	MatProps.kD[3] = 1.0E-7; // m/s

	// /!\ for a yet unknwon reason cohesion <100E6 gives a dirty viscosity jump at the interface with the sticky air


	Physics.Cp = 1.0;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	BCStokes.SetupType = Sandbox;
	BCStokes.backStrainRate = -1.0E-14;//+0.00001;

	BCThermal.TT = 0.0;
	BCThermal.TB = 0.0;

	compute dtMax = 3600*24*365.25 * 1E6;
	compute time = 0;
	Physics.dt = 3600*24*365.25 * 100E6; // initial value is really high to set the temperature profile. Before the advection, dt is recomputed to satisfy CFL
	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);

	Physics.g[0] = -9.81*sin( 0*PI/180);
	Physics.g[1] = -9.81*cos( 0*PI/180);



	Darcy.hOcean = Grid.ymin + (Grid.ymax-Grid.ymin)*0.42;
	Darcy.rainFlux = 0.2/(3600.0*24.0*365.0);








	compute CFL_fac = 1.0; // 0.5 ensures stability
	Particles.noiseFactor = 0.0; // between 0 and 1

	Visu.type 			= WaterPressureHead; // Default
	Visu.typeParticles	= Phase; // Default
	Visu.showParticles  = true;
	Visu.shiftFac[0]    = 0.0;
	Visu.shiftFac[1] 	= -.51;
	Visu.shiftFac[2] 	= -.05;
	Visu.writeImages 	= true;
	Visu.transparency 	= false;
	//Visu.outputFolder 	= "../StokesFD/OutputTest/";
	strcpy(Visu.outputFolder, "../StokesFD_OutputTest/");

	Visu.retinaScale = 2;

	// Initialize some arrays for comparing with numerical and analytical solutions
	/*
	compute* NumericalSolution  = (compute*) malloc(nTimeSteps * sizeof(compute));
	compute* AnalyticalSolution = (compute*) malloc(nTimeSteps* sizeof(compute));
	compute* timeArray 			= (compute*) malloc(nTimeSteps* sizeof(compute));
	*/
	compute a[20];

	//EqSystem.penaltyMethod = false;
	//EqSystem.penaltyFac = 1000000;


	// Set characteristic quantities
	// =================================
	Char.length 		= (Grid.ymax-Grid.ymin)*0.5    ;//fmin(Grid.dx,Grid.dy);
	Char.density 		= 0.5*(MatProps.rho0[0]+MatProps.rho0[1]);
	Char.acceleration 	= fabs(Physics.g[1]);

	Char.viscosity  	= 0.5*(MatProps.eta0[2]+MatProps.eta0[1]);//pow( 10, (log10(MatProps.eta0[0])+log10(MatProps.eta0[1]))/2 );

	//Char.stress 		= 2.0*fabs(BCStokes.backStrainRate)*Char.viscosity;



	//Char.viscosity 		= Char.stress/(2.0*fabs(BCStokes.backStrainRate));
	//Char.viscosity 		= Char.stress/(2.0*fabs(BCStokes.backStrainRate));


	Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h
	//Char.stress 		= 2*Char.viscosity*fabs(BCStokes.backStrainRate); // i.e. rho*g*h



	Char.time 			= Char.viscosity/Char.stress;
	Char.velocity 		= Char.length/Char.time;
	//Char.acceleration 	= Char.length/Char.time/Char.time;
	Char.strainrate 	= 1.0/Char.time;
	Char.mass			= Char.density*Char.length*Char.length*Char.length;

	//Char.temperature 	= (BCThermal.TB+BCThermal.TT)*0.5;
	Char.temperature 	= (BCThermal.TB);
	if (Char.temperature == 0)
		Char.temperature = 1;

	MatProps.maxwellTime[0] = MatProps.eta0[0]/MatProps.G[0];
	MatProps.maxwellTime[1] = MatProps.eta0[1]/MatProps.G[1];
	MatProps.maxwellTime[2] = MatProps.eta0[1]/MatProps.G[2];
	MatProps.maxwellTime[3] = MatProps.eta0[1]/MatProps.G[3];
	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes, &BCThermal, &Darcy);

	printf("Eta0[1] = %.3e", MatProps.eta0[1]);

	Physics.etaMin = 1E-5;
	Physics.etaMax = 1E3;
	Physics.epsRef = fabs(BCStokes.backStrainRate);

	printf("max backStrainRate = %.3e\n",BCStokes.backStrainRate);
	if (Physics.epsRef == 0)
		Physics.epsRef = 1E0;




	// Determine dtmin and dt max based on maxwell times on materials

	compute dtmin = 1E100;
	compute dtmax = 0;
	for (i=0;i<MatProps.nPhase;++i) {
		if (MatProps.maxwellTime[i]<dtmin) {
			dtmin = MatProps.maxwellTime[i];
		}
		if (MatProps.maxwellTime[i]>dtmax) {
			dtmax = MatProps.maxwellTime[i];
		}
		printf("maxwell time = %.2e\n", MatProps.maxwellTime[i]);
	}
	dtmin = 1E-100;
	//dtmax = dtmax;

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

	Visu.ntri   	= 2;//Grid.nxC*Grid.nyC*2;
	Visu.ntrivert 	= Visu.ntri*3;
	Visu.nParticles = Particles.n+ (int) (Particles.n*0.1); // overallocate 5% of the number of particles
	Visu.particleMeshRes = 6;

	printf("xmin = %.3f, ymin = %.3f\n", Grid.xmin, Grid.ymin);


	coord xmax_ini = Grid.xmax;
	coord xmin_ini = Grid.xmin;
	coord ymax_ini = Grid.ymax;
	coord ymin_ini = Grid.ymin;


	printf("kD = %.2e, SD=%.2e\n", MatProps.kD[0]*Char.length/Char.time,MatProps.SD[0]/Char.length);








	// Other variables
	// =================================
	int iEq, iLS;

	// Allocate memory
	// =================================
	printf("Allocate memory\n");
	Memory_allocateMain(&Grid, &Particles, &Physics, &EqStokes, &NumStokes, &NumThermal);


	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes	(&BCStokes , &Grid, &EqStokes);
	BC_initThermal	(&BCThermal, &Grid, &EqThermal);


	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_allocateI(&EqStokes);
	Numbering_init(&BCStokes, &Grid, &EqStokes, &NumStokes);


	printf("Numbering: init Thermal\n");
	EqSystem_allocateI(&EqThermal);
	Numbering_init(&BCThermal, &Grid, &EqThermal, &NumThermal);


	// Initialize Particles' coordinates
	// =================================
	printf("Particles: Init Coord\n");
	Particles_initCoord(&Grid, &Particles);
	printf("Particles: Init Coord\n");
	Particles_updateLinkedList(&Grid, &Particles, &Physics); // in case a ridiculous amount of noise is put on the particle
	printf("Particles: Init Coord\n");


	// Initialize Particles' phase
	// =================================
	printf("Particles: Init Phase\n");
	Particles_initPhase(&Grid, &Particles);

	// Initialize Particles' passive
	// =================================
	printf("Particles: Init Passive\n");
	Particles_initPassive(&Grid, &Particles);



	// Get Physics from particles to cell and to nodes (important for Neumann conditions)
	// =================================
	printf("Physics: Interp from particles to cell\n");
	printf("kD = %.2e, SD=%.2e\n", MatProps.kD[0]*Char.length/Char.time,MatProps.SD[0]/Char.length);
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);
	Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);
	Physics_interpFromCellToNode(&Grid, Physics.G  , Physics.GShear  );


	// Allocate memory for the system of equations
	// =================================
	EqSystem_allocateMemory(&EqStokes );
	EqSystem_allocateMemory(&EqThermal);



	// Init Solver
	// =================================
	printf("EqStokes: Init Solver\n");
	EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes); // dummy assembly to give the EqSystem initSolvers
	//EqSystem_check(&EqStokes);
	EqSystem_initSolver (&EqStokes, &SolverStokes);


	printf("EqThermal: Init Solver\n");
	EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (&EqThermal, &SolverThermal);

	// Initial temperature profile
	EqSystem_solve(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);
	Physics_set_T_FromSolution(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);
	Physics_interpTempFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal);
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);


	// Initial Darcy profile
	Physics.dt = 1E4*3600*24*365/Char.time; // run Darcy for this period of time initially
	printf("Darcy: solve\n");
	Darcy_solve(&Darcy, &Grid, &Physics, &MatProps, &Particles);
	printf("Darcy: interp from cells to particles\n");
	Physics_interpPsiFromCellsToParticle(&Grid, &Particles, &Physics);

	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          	INIT VISUALIZATION                            //
	//                                                                            //
	//============================================================================//
	//============================================================================//

#if (VISU)
	// Init GLFW
	// =======================================
	GLFWwindow* window = NULL;
	Visu_initWindow(&window, &Visu);
	Visu_allocateMemory(&Visu, &Grid);
	Visu_init(&Visu, &Grid, &Particles);
	Visu_initOpenGL(&Visu, &Grid, window);
#endif














//======================================================================================================//
//======================================================================================================//
//                                                                      				      			//
//                          				TIME LOOP             		  					        	//
//                                                                      							    //
//======================================================================================================//
//======================================================================================================//
	int timeStep = 0;
	Physics.dt = dtmax*1000;// pow(10,(log10(dtmin)+log10(dtmax))/2);
	Physics.time = 0;
	Physics.itNonLin = -1;
	Physics.glob = 1.0;
	//printf("dt ini = %.2e, meandt = %.2e\n", Physics.dt, (dtmax));
	while(timeStep!=nTimeSteps) {

		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                          COMPUTE AND UPDATE STOKES                         //
		//                                                                            //
		//============================================================================//
		//============================================================================//
		printf("\n\n\n          ========  Time step %i  ========   \n"
				"       =====================================\n\n",timeStep);

		// Get Physics from particles to cell and to nodes
		// =================================
		TIC


		Physics_computeEta(&Physics, &Grid);

		printf("Physics: Interp from cell to node\n");
		Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);
		Physics_interpFromCellToNode(&Grid, Physics.G  , Physics.GShear  );
		TOC
		printf("Physics: Interp took %.2fs\n",toc);

		// Update BC
		// =================================
		printf("BC: Update\n");
		BC_updateStokes(&BCStokes, &Grid);
		BC_updateThermal(&BCThermal, &Grid);








		// Solve the heat conservation
		// =================================
		EqSystem_solve(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);
		Physics_set_T_FromSolution(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);
		// update temperature on markers
		Physics_interpTempFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal);


		// Assemble the systems of equations
		// =================================
		printf("EqStokes: Assemble\n");
		EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);
		EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal);




		// ============================================================================
		// 							Non-linear interation
		// ============================================================================
		Physics.itNonLin = 0;
		EqStokes.normResidual = 1.0;
		compute normRes0 = 1.0;
		compute normResRef = 1.0;
		compute* NonLin_x0 = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		compute* NonLin_dx = (compute*) malloc(EqStokes.nEq * sizeof(compute));




		double timeStepTic = glfwGetTime();
		while(( (EqStokes.normResidual/normRes0 > relativeTolerance && EqStokes.normResidual/normResRef > absoluteTolerance ) && Physics.itNonLin!=maxNonLinearIter ) || Physics.itNonLin<minNonLinearIter) {

			printf("\n\n  ==== Non linear iteration %i ==== \n\n",Physics.itNonLin);

			// Solve: A(X0) * X = b
			// ====================
			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_x0[iEq] = EqStokes.x[iEq];
			}
			printf("EqStokes: Solve\n");


			if (fabs(Physics.maxV)<1E-6)
				Physics.maxV = 1E-6;
			Physics.dt = CFL_fac*fmin(Grid.dx,Grid.dy)/(Physics.maxV); // note: the min(dx,dy) is the char length, so = 1
			printf("maxV = %.3em, Physics.dt = %.3e, Physics.dt(SCALED)= %.3e yr, dtmin = %.2e, dtmax = %.2e, dtMax = %.2e\n",fabs(Physics.maxV), Physics.dt, Physics.dt*Char.time/3600/24/365, dtmin, dtmax, dtMax);


			if (Physics.dt<dtmin) {
				Physics.dt = dtmin;
			} else if (Physics.dt>dtmax) {
			//	Physics.dt = dtmax;
			}

			EqSystem_solve(&EqStokes, &SolverStokes, &Grid, &Physics, &BCStokes, &NumStokes);




			// ======================================
			// 				Line search
			// ======================================

			a[nLineSearch] = 1.0/nLineSearch; // this is the best value
			compute minRes = 1E100;

			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes.x[iEq] - NonLin_x0[iEq];
			}

			for (iLS= 0; iLS < nLineSearch+1; ++iLS) {
				printf("== Line search %i:  ", iLS);

				//compute a, the globalization parameter;
				if (iLS!=nLineSearch)
					a[iLS] = maxCorrection - maxCorrection/(nLineSearch) * (iLS);

				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					// X1 = X0 + a*(X-X0)
					EqStokes.x[iEq] = NonLin_x0[iEq] + a[iLS]*(NonLin_dx[iEq]);
				}

				Physics.glob = a[iLS];

				// Update the stiffness matrix
				Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);



				Physics_computeEta(&Physics, &Grid);



				Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);
				Physics_interpFromCellToNode(&Grid, Physics.G  , Physics.GShear  );
				TIC
				EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);

				TOC
				printf("Stokes Assembly: %.2fs\n", toc);


				// compute the norm of the  residual:
				// F = b - A(X1) * X1
				EqSystem_computeNormResidual(&EqStokes);
				printf("a = %.2f, |F| / |b|: %.3e, minRes = %.3e, best a = %.2f\n", a[iLS], EqStokes.normResidual/normResRef, minRes, a[nLineSearch]);
				if (EqStokes.normResidual<minRes) {
					a[nLineSearch] = a[iLS];
					minRes = EqStokes.normResidual;
					if (iLS==nLineSearch-1) // if the last one is the best one then don't recompute, i.e. next one would be the same
						break;
				}

				if (Physics.itNonLin==0) {
					normRes0 = EqStokes.normResidual;
					if (timeStep==0)
						normResRef = EqStokes.normResidual;

					if (maxNonLinearIter==1) {
						break;
					}
				}

			}
			// ======================================
			// 		   End of Line search
			// ======================================

			// Blowing up check: if the residual is too large
			// wipe up the solution vector and start the iteration again with 0 everywhere initial guess
			if (minRes>1000.0) {
				for (i=0; i<EqStokes.nEq; i++) {
					EqStokes.x[i] = 0;
				}
				Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Physics.itNonLin = 0;
				printf("/! /!  Warning  /! /! : The residual is larger than the tolerance. The non linear iterations might be diverging. Wiping up the solution and starting the iteration again\n");
			}



			Physics.itNonLin++;
		} // end of non-linear loop

		free(NonLin_x0);
		free(NonLin_dx);

		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);



		// ======================================
		// 		   Solve Darcy
		// ======================================
		printf("Darcy: solve\n");
		Darcy_solve(&Darcy, &Grid, &Physics, &MatProps, &Particles);
		printf("Darcy: interp from cells to particles\n");
		Physics_interpPsiFromCellsToParticle(&Grid, &Particles, &Physics);


		// ============================================================================
		// 						End of Non-linear interation
		// ============================================================================









		// ============================================================================
		// 									Advect
		// ============================================================================
		// Update dt



		/*
		if (fabs(Physics.maxV)<1E-6)
			Physics.maxV = 1E-6;
		Physics.dt = CFL_fac*fmin(Grid.dx,Grid.dy)/(Physics.maxV); // note: the min(dx,dy) is the char length, so = 1
		printf("maxV = %.3em Physics.dt = %.3e, dtmin = %.2e, dtmax = %.2e, dtMax = %.2e\n",fabs(Physics.maxV), Physics.dt, dtmin, dtmax, dtMax);


		if (Physics.dt<dtmin) {
			Physics.dt = dtmin;
		} else if (Physics.dt>dtmax) {
		//	Physics.dt = dtmax;
		}
		*/





		Physics.time += Physics.dt;
		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);


		// update stress on the particles
		Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);

		/*
		printf("DSigma_xx \n");
		C = 0;
		for (iy = 0; iy < Grid.nyEC; ++iy) {
			for (ix = 0; ix < Grid.nxEC; ++ix) {
				printf("%.3e  ", Physics.Dsigma_xx_0[C]);
				C++;
			}
			printf("\n");
		}
		printf("DSigma_xy \n");
		C = 0;
		for (iy = 0; iy < Grid.nyS; ++iy) {
			for (ix = 0; ix < Grid.nxS; ++ix) {
				printf("%.3e  ", Physics.Dsigma_xy_0[C]);
				C++;
			}
			printf("\n");
		}
		*/

		Physics_interpStressesFromCellsToParticle(&Grid, &Particles, &Physics, &BCStokes,  &BCThermal, &NumThermal);




		printf("Particles: Advect\n");
		Particles_advect(&Particles, &Grid, &Physics);
		//printf("Grid: Update pure shear\n");
		switch (BCStokes.SetupType) {
		case PureShear:
			Grid_updatePureShear(&Grid, &BCStokes, Physics.dt);
			Particles_teleportInsideTheDomain(&Grid, &Particles, &Physics);
			break;
		case SimpleShearPeriodic:
			Particles_Periodicize(&Grid, &Particles, &BCStokes);
			break;
		case FixedLeftWall:
			break;
		case Sandbox:
			Grid_updatePureShear(&Grid, &BCStokes, Physics.dt);
			Particles_teleportInsideTheDomain(&Grid, &Particles, &Physics);
			//Particles_deleteIfOutsideTheDomain(&Grid, &Particles);
			break;
		default:
			break;
		}

		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(&Grid, &Particles, &Physics);
		printf("Physics: Interp from particles to cell\n");
		Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes, &NumThermal, &BCThermal);
		// ============================================================================
		// 								End of	Advect
		// ============================================================================








		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                     END OF COMPUTE AND UPDATE STOKES                       //
		//                                                                            //
		//============================================================================//
		//============================================================================//






		printf("Physics EpsRef = %.3e\n", Physics.epsRef);








#if VISU
		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                                 VISUALIZATION                              //
		//                                                                            //
		//============================================================================//
		//============================================================================//
		GLfloat shiftIni[3];
		do  {
				glfwPollEvents();
				///printf("A-3\n");
				Visu_checkInput(&Visu, window);
				//printf("A-4\n");
				glClearColor(0, 0, 0, 1); // black
				//glClear(GL_COLOR_BUFFER_BIT);
				//printf("A-1\n");
				glEnable(GL_DEPTH_TEST);
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
				glStencilMask(0xFF);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

				glStencilMask(0x00);



				//glClear(GL_COLOR_BUFFER_BIT);

				if (Visu.initPassivePart) {
					Particles_initPassive(&Grid, &Particles);
					Visu.initPassivePart = false;
				}

				shiftIni[0] = Visu.shift[0];
				shiftIni[1] = Visu.shift[1];
				shiftIni[2] = Visu.shift[2];

				Visu.shift[0] -= (xmax_ini-xmin_ini)*Visu.shiftFac[0]*Visu.scale;
				Visu.shift[1] += (ymax_ini-ymin_ini)*Visu.shiftFac[1]*Visu.scale;
				Visu.shift[2] +=                 1.0*Visu.shiftFac[2];

				//============================================================================
				// 								PLOT PARTICLE

				if (Visu.showParticles) {
					glEnable(GL_STENCIL_TEST);
					// Draw the box in black and use it to set the stencil
					// Particles fragment will be drawn only where the stencil is 1 (i.e. inside the box)
					glStencilFunc(GL_ALWAYS, 1, 0xFF); // All fragments should update the stencil buffer
					glStencilMask(0xFF); // Enable writing to the stencil buffer

					glBindVertexArray(Visu.VAO);
					glUseProgram(Visu.ParticleBackgroundShaderProgram);
					glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu.EBO);
					Visu_updateUniforms(&Visu, window);
						//Visu_update(&Visu, window, &Grid, &Physics, &BCStokes, &Char);
						glDrawElements(GL_TRIANGLES, Visu.ntrivert, GL_UNSIGNED_INT, 0);

					glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
					glUseProgram(0);
					glBindVertexArray(0);


					glStencilFunc(GL_EQUAL, 1, 0xFF);
					glStencilMask(0x00);  // disable writing to the buffer
					//glDisable(GL_DEPTH_TEST);



					glBindVertexArray(Visu.VAO_part);
					glUseProgram(Visu.ParticleShaderProgram);


					glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_part);

					//glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_partMesh);
					// update the buffer containing the particles
						Visu_particles(&Visu, &Particles, &Grid);
						Visu_updateUniforms(&Visu, window);
						glBufferSubData(GL_ARRAY_BUFFER, 0, 4*Particles.n*sizeof(GLfloat), Visu.particles);
						glBindBuffer(GL_ARRAY_BUFFER, 0);

						//glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_partMesh);
						//glBindBuffer(GL_ARRAY_BUFFER, 0);
						glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, Visu.particleMeshRes+2, Particles.n);
						//printf("Visu.particleMeshRes= %i\n",Visu.particleMeshRes);
						//glDrawArraysInstanced(GL_TRIANGLES, 0, 3, Particles.n);
						//
					glUseProgram(0);
					glBindVertexArray(0);
					glDisable(GL_STENCIL_TEST);
				}
				// 								PLOT PARTICLE
				//============================================================================



				Visu.shift[0] += 2*(xmax_ini-xmin_ini)*Visu.shiftFac[0]*Visu.scale;
				Visu.shift[1] -= 2*(ymax_ini-ymin_ini)*Visu.shiftFac[1]*Visu.scale;
				Visu.shift[2] -=                   2.0*Visu.shiftFac[2];

				//============================================================================
				// 								PLOT GRID DATA


				// ****** Bind shader textures, arrays and buffers
				glBindVertexArray(Visu.VAO);
				glUseProgram(Visu.ShaderProgram);
				glBindTexture(GL_TEXTURE_2D, Visu.TEX);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu.EBO);

					// 1. Update data
					if (BCStokes.SetupType==PureShear || BCStokes.SetupType==Sandbox) {
						Visu_updateVertices(&Visu, &Grid);
						glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO);
								glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu.vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}
					Visu_update(&Visu, window, &Grid, &Physics, &BCStokes, &Char);
					Visu_alphaValue(&Visu, &Grid, &Particles);
					// update the content of Visu.U
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RG, Grid.nxS, Grid.nyS, 0, GL_RG, GL_FLOAT, Visu.U);	// load the updated Visu.U in the texture
					// 2. Draw
					glDrawElements(GL_TRIANGLES, Visu.ntrivert, GL_UNSIGNED_INT, 0);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
				glBindTexture(GL_TEXTURE_2D, 0);
				glUseProgram(0);
				glBindVertexArray(0);



				// ****** Unbind textures, arrays and buffers


				// 								PLOT GRID DATA
				//============================================================================


				Visu.shift[0] = shiftIni[0];
				Visu.shift[1] = shiftIni[1];
				Visu.shift[2] = shiftIni[2];



				//============================================================================
				// 							  SAVE TO IMAGE FILE

				if (Visu.writeImages) {
					FILE *fptr;
					char fname[1024];
					char ftitle[1024];
					sprintf(fname,"%sFrame_%05i.png",Visu.outputFolder,timeStep);
					sprintf(ftitle,"time_%5.5e.png",time);
					//sprintf(fname,"Frame_%04i.raw",timeStep);
					if ((fptr = fopen(fname,"w")) == NULL) {
						fprintf(stderr,"Failed to open the file for window dump\n");
						exit(0);
					}

					glPixelStorei(GL_PACK_ALIGNMENT,1);
					glReadBuffer(GL_BACK);
					glReadPixels(0,0,Visu.retinaScale*WIDTH,Visu.retinaScale*HEIGHT,GL_RGB,GL_UNSIGNED_BYTE,Visu.imageBuffer);

					//fwrite(Visu.imageBuffer,WIDTH*HEIGHT*3,1,fptr);
					int result = writePNGImage(fname, Visu.retinaScale*WIDTH, Visu.retinaScale*HEIGHT, Visu.imageBuffer, ftitle);
					if (result!=0) {
						printf("error: couldn't write png file\n");
						exit(0);
					}
					fclose(fptr);


				}


				// 							  SAVE TO IMAGE FILE
				//============================================================================


				glfwSwapBuffers(window);

				if (glfwWindowShouldClose(window))
					break;
				if (timeStep==nTimeSteps-1)
					Visu.paused = true;
			} while (Visu.paused);

		if (glfwWindowShouldClose(window))
			break;

		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                          END OF VISUALIZATION                              //
		//                                                                            //
		//============================================================================//
		//============================================================================//

#endif
		/*
		int I = 2+2*Grid.nxEC;
		NumericalSolution[timeStep] = Physics.sigma_xx_0[I];
		AnalyticalSolution[timeStep] = 2*BCStokes.backStrainRate*Physics.eta[I]* (1-exp(-time*Physics.G[I]/Physics.eta[I]));
		timeArray[timeStep] = time;

		printf("time = %.2e\n", time);
		printf("epsref = %.2e, Physics.eta[I] = %.2e, Physics.G[I] = %.2e, Physics.sigma_xx_0[I] = %.2e\n", Physics.epsRef, Physics.eta[I], Physics.G[I], Physics.sigma_xx_0[I]);
		printf("Ana = %.2e\n", AnalyticalSolution[timeStep]);
		*/
		timeStep++;
	}







//======================================================================================================//
//======================================================================================================//
//                                                                      				      			//
//                          				END OF TIME LOOP           	  					        	//
//                                                                      							    //
//======================================================================================================//
//======================================================================================================//
	printf("Simulation successfully completed\n");

	/*

	printf("== Numerical Solution ==\n");
	for (i = 0; i < nTimeSteps; ++i) {
		printf("%.2e  ",NumericalSolution[i]);
	}
	printf("\n");
	printf("== Analytical Solution ==\n");
	for (i = 0; i < nTimeSteps; ++i) {
		printf("%.2e  ",AnalyticalSolution[i]);
	}
	printf("\n");
	printf("== TimeArray ==\n");
	for (i = 0; i < nTimeSteps; ++i) {
		printf("%.2e  ",timeArray[i]);
	}
	printf("\n");

	*/

	//printListd(NumericalSolution, nTimeSteps);


	//printListd(AnalyticalSolution, nTimeSteps);


	//printListd(timeArray, nTimeSteps);


	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	// Free memory
	printf("Start free\n");
	Memory_freeMain(&Particles, &Physics, &NumStokes, &NumThermal, &BCStokes, &Grid);
	EqSystem_freeMemory(&EqStokes, &SolverStokes);

	/*
	free(NumericalSolution);
	free(AnalyticalSolution);
	free(timeArray);
	*/

#if VISU
	// Quit glfw
	glfwDestroyWindow(window);
	glfwTerminate();
	Visu_freeMemory(&Visu);
#endif

	printf("SUCCESS!\n\n\n");

	return EXIT_SUCCESS;
}

