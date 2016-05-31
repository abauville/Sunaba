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

	// Numerics
	Numerics Numerics;


	INIT_TIMER



	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                       		 USER INPUT      		   	              	  //
	//                                                                            //
	//============================================================================//
	//============================================================================//


	// Set model properties
	// =================================
	Numerics.nTimeSteps  = 1; //  negative value for infinite
	Numerics.nLineSearch = 4;
	Numerics.maxNonLinearIter = 10; // should always be greater than the number of line searches
	Numerics.minNonLinearIter = 5; // should always be greater than the number of line searches
	Numerics.relativeTolerance = 3E-5; // relative tolerance to the one of this time step
	Numerics.absoluteTolerance = 3E-5; // relative tolerance to the first one of the simulation
	Numerics.maxCorrection = 1.0;


	Grid.nxC = 256;
	Grid.nyC = 128;

	Particles.nPCX = 4;	Particles.nPCY = 4;

	//Grid.xmin = 0;
	//Grid.xmax = (compute) Grid.nxC;
	//Grid.ymin = 0;
	//Grid.ymax = (compute) Grid.nyC;
	Grid.xmin =  -6*50E3;
	Grid.xmax =   0*50E3;
	Grid.ymin =   0*50E3;
	Grid.ymax =   1*50E3;

	MatProps.nPhase  = 4;

	MatProps.rho0[0] = 10; 			MatProps.eta0[0] = 1E18;  		MatProps.n[0] = 1.0;
	MatProps.rho0[1] = 1000; 		MatProps.eta0[1] = 1E18;  		MatProps.n[1] = 1.0;
	MatProps.rho0[2] = 2700;		MatProps.eta0[2] = 1E23; 		MatProps.n[2] = 1.0;
	MatProps.rho0[3] = 2700;		MatProps.eta0[3] = 1E20; 		MatProps.n[3] = 1.0;
	MatProps.rho0[4] = 2700;		MatProps.eta0[4] = 1E23; 		MatProps.n[4] = 1.0;

	MatProps.alpha[0] = 1E-5;  	MatProps.beta [0] = 0.0;  		MatProps.k[0] = 1E-2; 			MatProps.G[0] = 1E11;
	MatProps.alpha[1] = 1E-5;  	MatProps.beta [1] = 0.0;  		MatProps.k[1] = 1E-2; 			MatProps.G[1] = 1E11;
	MatProps.alpha[2] = 1E-5; 	MatProps.beta [2] = 0.0;  		MatProps.k[2] = 1E-2; 			MatProps.G[2] = 1E11;
	MatProps.alpha[3] = 1E-5; 	MatProps.beta [3] = 0.0;  		MatProps.k[3] = 1E-2; 			MatProps.G[3] = 1E11;
	MatProps.alpha[4] = 1E-5; 	MatProps.beta [4] = 0.0;  		MatProps.k[4] = 1E-2; 			MatProps.G[4] = 1E11;

	MatProps.cohesion[0] = 10000.0*1E6; 	MatProps.frictionAngle[0] = 30*PI/180; //air
	MatProps.cohesion[1] = 10000.0*1E6; 	MatProps.frictionAngle[1] = 30*PI/180; //air
	MatProps.cohesion[2] = 10.0*1E6;		MatProps.frictionAngle[2] = 30*PI/180; // green
	MatProps.cohesion[3] = 10.0*1E6;		MatProps.frictionAngle[3] = 5*PI/180;  // orange
	MatProps.cohesion[4] = 100.0*1E6;		MatProps.frictionAngle[4] = 30*PI/180; // blue


	MatProps.SD[0] = 1E-3; // 1/m
	MatProps.SD[1] = 1E-3; // 1/m
	MatProps.SD[2] = 1E-3; // 1/m
	MatProps.SD[3] = 1E-3; // 1/m
	MatProps.SD[4] = 1E-3; // 1/m

	MatProps.kD[0] = 1.0E-7; // m/s
	MatProps.kD[1] = 1.0E-7; // m/s
	MatProps.kD[2] = 1.0E-7; // m/s
	MatProps.kD[3] = 1.0E-7; // m/s
	MatProps.kD[4] = 1.0E-7; // m/s



	Physics.Cp = 1.0;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	BCStokes.SetupType = Sandbox;
	BCStokes.backStrainRate = -1.0E-14;//+0.00001;

	BCThermal.TT = 0.0;
	BCThermal.TB = 0.0;

	compute dtMax = 3600*24*365.25 * 1E6;
	Physics.dt = 3600*24*365.25 * 100E6; // initial value is really high to set the temperature profile. Before the advection, dt is recomputed to satisfy CFL
	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);

	Physics.g[0] = -9.81*sin( 0*PI/180);
	Physics.g[1] = -9.81*cos( 0*PI/180);

	compute CFL_fac = 0.4; // 0.5 ensures stability
	Particles.noiseFactor = 0.3; // between 0 and 1

	Visu.type 			= StrainRate; // Default
	Visu.typeParticles	= Phase; // Default
	Visu.showParticles  = true;
	Visu.shiftFac[0]    = 0.00;
	Visu.shiftFac[1] 	= 0.00;
	Visu.shiftFac[2] 	= +.05;
	Visu.writeImages 	= true;
	Visu.transparency 	= true;
	Visu.alphaOnValue 	= true;
	Visu.showGlyphs 	= true;
	Visu.glyphType		= StokesVelocity;
	Visu.glyphMeshType	= ThinArrow;
	Visu.glyphScale		= 0.05;
	Visu.glyphSamplingRateX  = 3; // sample every Visu.glyphSampling grid points
	Visu.glyphSamplingRateY  = 6; // sample every Visu.glyphSampling grid points

	//Visu.outputFolder 	= "../StokesFD/OutputTest/";
	strcpy(Visu.outputFolder, "../StokesFD_OutputTest/");

	Visu.retinaScale = 2;


	compute a[20];



	// Set characteristic quantities
	// =================================
	Char.length 		= (Grid.ymax-Grid.ymin)*0.5    ;//fmin(Grid.dx,Grid.dy);
	Char.density 		= 0.5*(MatProps.rho0[0]+MatProps.rho0[1]);
	Char.acceleration 	= fabs(Physics.g[1]);

	Char.viscosity  	= 0.5*(MatProps.eta0[2]+MatProps.eta0[1]);//pow( 10, (log10(MatProps.eta0[0])+log10(MatProps.eta0[1]))/2 );


	Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h



	Char.time 			= Char.viscosity/Char.stress;
	Char.velocity 		= Char.length/Char.time;
	Char.strainrate 	= 1.0/Char.time;
	Char.mass			= Char.density*Char.length*Char.length*Char.length;

	Char.temperature 	= (BCThermal.TB);
	if (Char.temperature == 0)
		Char.temperature = 1;

	MatProps.maxwellTime[0] = MatProps.eta0[0]/MatProps.G[0];
	MatProps.maxwellTime[1] = MatProps.eta0[1]/MatProps.G[1];
	MatProps.maxwellTime[2] = MatProps.eta0[1]/MatProps.G[2];
	MatProps.maxwellTime[3] = MatProps.eta0[1]/MatProps.G[3];
	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes, &BCThermal);

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


	Grid.xmax_ini = Grid.xmax;
	Grid.xmin_ini = Grid.xmin;
	Grid.ymax_ini = Grid.ymax;
	Grid.ymin_ini = Grid.ymin;


	printf("kD = %.2e, SD=%.2e\n", MatProps.kD[0]*Char.length/Char.time,MatProps.SD[0]/Char.length);






	//======================================================================================================//
	//======================================================================================================//
	//                                                                      				      			//
	//                          				INITIALIZATION         		  					        	//
	//                                                                      							    //
	//======================================================================================================//

	// Other variables
	// =================================
	int iEq, iLS;

	// Allocate memory
	// =================================
	printf("Allocate memory\n");
	Physics_allocateMemory(&Physics, &Grid);


	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes	(&BCStokes , &Grid, &EqStokes);
	BC_initThermal	(&BCThermal, &Grid, &EqThermal);


	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_allocateI(&EqStokes);
	Numbering_allocateMemory(&NumStokes, &EqStokes, &Grid);
	Numbering_init(&BCStokes, &Grid, &EqStokes, &NumStokes);


	printf("Numbering: init Thermal\n");
	EqSystem_allocateI(&EqThermal);
	Numbering_allocateMemory(&NumThermal, &EqThermal, &Grid);
	Numbering_init(&BCThermal, &Grid, &EqThermal, &NumThermal);


	// Initialize Particles' coordinates
	// =================================
	Particles_allocateMemory(&Particles, &Grid);
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

	Visu_allocateMemory(&Visu, &Grid);
	Visu_init(&Visu, &Grid, &Particles);
#endif






	//======================================================================================================//
	//======================================================================================================//
	//                                                                      				      			//
	//                          				TIME LOOP             		  					        	//
	//                                                                      							    //
	//======================================================================================================//
	//======================================================================================================//
	Numerics.timeStep = 0;
	Physics.dt = dtmax*1000;// pow(10,(log10(dtmin)+log10(dtmax))/2);
	Physics.time = 0;
	Numerics.itNonLin = -1;
	Numerics.glob = 1.0;
	//printf("dt ini = %.2e, meandt = %.2e\n", Physics.dt, (dtmax));


	while(Numerics.timeStep!=Numerics.nTimeSteps) {

		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                          COMPUTE AND UPDATE STOKES                         //
		//                                                                            //
		//============================================================================//
		//============================================================================//
		printf("\n\n\n          ========  Time step %i  ========   \n"
				"       ===================================== \n\n",Numerics.timeStep);

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
		Numerics.itNonLin = 0;
		EqStokes.normResidual = 1.0;
		compute normRes0 = 1.0;
		compute normResRef = 1.0;
		compute* NonLin_x0 = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		compute* NonLin_dx = (compute*) malloc(EqStokes.nEq * sizeof(compute));




		double timeStepTic = glfwGetTime();
		while(( (EqStokes.normResidual/normRes0 > Numerics.relativeTolerance && EqStokes.normResidual/normResRef > Numerics.absoluteTolerance ) && Numerics.itNonLin!=Numerics.maxNonLinearIter ) || Numerics.itNonLin<Numerics.minNonLinearIter) {

			printf("\n\n  ==== Non linear iteration %i ==== \n\n",Numerics.itNonLin);

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

			a[Numerics.nLineSearch] = 1.0/Numerics.nLineSearch; // this is the best value
			compute minRes = 1E100;

			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes.x[iEq] - NonLin_x0[iEq];
			}

			for (iLS= 0; iLS < Numerics.nLineSearch+1; ++iLS) {
				printf("== Line search %i:  ", iLS);

				//compute a, the globalization parameter;
				if (iLS!=Numerics.nLineSearch)
					a[iLS] = Numerics.maxCorrection - Numerics.maxCorrection/(Numerics.nLineSearch) * (iLS);

				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					// X1 = X0 + a*(X-X0)
					EqStokes.x[iEq] = NonLin_x0[iEq] + a[iLS]*(NonLin_dx[iEq]);
				}

				Numerics.glob = a[iLS];

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
				printf("a = %.2f, |F| / |b|: %.3e, minRes = %.3e, best a = %.2f\n", a[iLS], EqStokes.normResidual/normResRef, minRes, a[Numerics.nLineSearch]);
				if (EqStokes.normResidual<minRes) {
					a[Numerics.nLineSearch] = a[iLS];
					minRes = EqStokes.normResidual;
					if (iLS==Numerics.nLineSearch-1) // if the last one is the best one then don't recompute, i.e. next one would be the same
						break;
				}

				if (Numerics.itNonLin==0) {
					normRes0 = EqStokes.normResidual;
					if (Numerics.timeStep==0)
						normResRef = EqStokes.normResidual;

					if (Numerics.maxNonLinearIter==1) {
						break;
					}
				}

			}
			// ======================================
			// 		   End of Line search
			// ======================================

			// Blowing up check: if the residual is too large
			// wipe up the solution vector and start the iteration again with 0 everywhere initial guess

			if (Numerics.timeStep>1 && minRes>10.0) {
				for (i=0; i<EqStokes.nEq; i++) {
					EqStokes.x[i] = 0;
				}
				Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);
				Numerics.itNonLin = 0;
				printf("/! /!  Warning  /! /! : The residual is larger than the tolerance. The non linear iterations might be diverging. Wiping up the solution and starting the iteration again\n");
			}



			Numerics.itNonLin++;
		} // end of non-linear loop

		free(NonLin_x0);
		free(NonLin_dx);

		double timeStepToc = glfwGetTime();
		toc = timeStepToc-timeStepTic;
		printf("the timestep took: %.2f\n",toc);








		// ============================================================================
		// 						End of Non-linear interation
		// ============================================================================









		// ============================================================================
		// 									Advect
		// ============================================================================
		// Update dt
		Physics.time += Physics.dt;
		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);


		// update stress on the particles
		Physics_computeStressChanges  (&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);


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


		// Detect Faults
		Physics_changePhaseOfFaults(&Physics, &Grid, &MatProps, &Particles);

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
		Visu_main(&Visu, &Grid, &Physics, &Particles, &Numerics, &BCStokes, &Char, &MatProps);

		if (glfwWindowShouldClose(Visu.window))
			break;

#endif

		Numerics.timeStep++;
	}







	//======================================================================================================//
	//======================================================================================================//
	//                                                                      				      			//
	//                          				END OF TIME LOOP           	  					        	//
	//                                                                      							    //
	//======================================================================================================//
	//======================================================================================================//
	printf("Simulation successfully completed\n");


	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	// Free memory
	printf("Free Physics...\n");
	Physics_freeMemory(&Physics);
	printf("Free NumStokes...\n");
	Numbering_freeMemory(&NumStokes);
	printf("Free NumThermal...\n");
	Numbering_freeMemory(&NumThermal);
	printf("Free EqStokes...\n");
	EqSystem_freeMemory(&EqStokes, &SolverStokes);
	printf("Free EqThermal...\n");
	EqSystem_freeMemory(&EqThermal,&SolverThermal);
	printf("Free BCStokes...\n");
	BC_freeMemory(&BCStokes);
	printf("Free BCThermal...\n");
	BC_freeMemory(&BCThermal);
	printf("Free Particles...\n");
	Particles_freeMemory(&Particles, &Grid);
	printf("Memory freed successfully\n");


#if VISU
	// Quit glfw
	glfwDestroyWindow(Visu.window);
	glfwTerminate();
	Visu_freeMemory(&Visu);
#endif

	printf("SUCCESS!\n\n\n");

	return EXIT_SUCCESS;
}

