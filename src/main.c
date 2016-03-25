/*
 ============================================================================
 Name        : main.c
 Author      : Arthur Bauville
 Version     :
 Copyright   : 
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include "stokes.h"



int main(void) {

	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                              INIT STOKES PROGRAM                           //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	printf("\n\n\n\n\n\n");
	printf("               ============================\n"
		   "               ============================\n");
	printf("\n\n\n\n\n\nBeginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());
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


	INIT_TIMER

	// Set model properties
	// =================================
	int nTimeSteps  = 1; //  negative value for infinite
	int nLineSearch = 2;
	int maxNonLinearIter = 2;
	compute nonLinTolerance = 5E-4;

	Grid.nxC = 3;
	Grid.nyC = 4;

	Particles.nPCX = 3;
	Particles.nPCY = 4;

	//Grid.xmin = 0;
	//Grid.xmax = (compute) Grid.nxC;
	//Grid.ymin = 0;
	//Grid.ymax = (compute) Grid.nyC;
	Grid.xmin = -1.0;
	Grid.xmax =  1.0;
	Grid.ymin = -1.0;
	Grid.ymax =  1.0;

	MatProps.nPhase  = 2;
	MatProps.rho0[0] = 1; 		MatProps.eta0[0] = 1.0;  		MatProps.n[0] = 3.0; 		MatProps.flowLaw[0] = LinearViscous;
	MatProps.rho0[1] = 2.0;		MatProps.eta0[1] = 0.001; 		MatProps.n[1] = 3.0;		MatProps.flowLaw[1] = LinearViscous;



	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	BCStokes.SetupType = SimpleShearPeriodic;
	BCStokes.backStrainRate =  0.0;


	Physics.dt = 1.0;
	//Physics.epsRef = 1.0;//abs(BCStokes.backStrainRate);

	Physics.g[0] = 0.0;
	Physics.g[1] = -9.81;

	compute CFL_fac = 1.0; // 0.5 ensures stability
	Particles.noiseFactor = 1.0; // between 0 and 1

	//EqSystem.penaltyMethod = false;
	//EqSystem.penaltyFac = 1000000;


	// Set characteristic quantities
	// =================================
	Char.length 		= fmin(Grid.dx,Grid.dy);
	Char.density 		= MatProps.rho0[0];
	Char.acceleration 	= fabs(Physics.g[1]);
	Char.viscosity  	= MatProps.eta0[0];//pow( 10, (log10(MatProps.eta0[0])+log10(MatProps.eta0[0]))/2 );

	Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h
	//Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h

	Char.time 			= Char.viscosity/Char.stress;
	Char.velocity 		= Char.length/Char.time;
	Char.strainrate 	= 1.0/Char.time;

	Visu.type 			= Temperature; // Default
	Visu.showParticles  = false;

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BCStokes);
	Physics.etaMin = 1E-4;
	Physics.etaMax = 1E4;
	Physics.epsRef = abs(BCStokes.backStrainRate);
	if (Physics.epsRef == 0)
		Physics.epsRef = 1.0;

	// Init grid and particles
	// =================================
	Grid.nCTot  = Grid.nxC*Grid.nyC;


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
	Particles_updateLinkedList(&Grid, &Particles); // in case a ridiculous amount of noise is put on the particle
	printf("Particles: Init Coord\n");

	// Initialize Particles' phase
	// =================================
	printf("Particles: Init Phase\n");
	Particles_initPhase(&Grid, &Particles);

	// Initialize Particles' passive
	// =================================
	printf("Particles: Init Passive\n");
	Particles_initPassive(&Grid, &Particles);


	// Initialize the Physics stored on particles
	// =================================
	printf("Particles: Init Passive\n");
	Particles_initPhysics(&Grid, &Particles, &BCThermal);


	// Get Physics from particles to cell and to nodes (important for Neumann conditions)
	// =================================
	printf("Physics: Interp from particles to cell\n");
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes);
	printf("Physics: Interp from cell to node\n");
	Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BCStokes.SetupType);



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
	printf("EqThermal: Assembled successfully\n");
	EqSystem_check(&EqThermal);
	EqSystem_initSolver (&EqThermal, &SolverThermal);












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
	printf("Init Window\n");
	GLFWwindow* window = NULL;
	printf("Done\n");
	Visu.nParticles = Particles.n+ (int) (Particles.n*0.05); // overallocate 5% of the number of particles
	Visu.particleMeshRes = 6;
	Visu_initWindow(&window, &Visu);

	Visu_allocateMemory(&Visu, &Grid);
	Visu_init(&Visu, &Grid, &Particles);


	///Init shader
	// =======================================
	Visu.VertexShaderFile 			= "src/shader.vs";
	Visu.FragmentShaderFile 		= "src/shader.fs";
	Visu.ParticleVertexShaderFile 	= "src/particleShader.vs";
	Visu.ParticleGeometryShaderFile = "src/particleShader.gs";
	Visu.ParticleFragmentShaderFile = "src/particleShader.fs";


	Visu.ShaderProgram = 0;
	// Generate reference to objects (indexes that act as pointers to graphic memory)
	// =======================================
	Visu.VAO = 0; // Reference to the Vertex   array object
	Visu.VBO = 0; // Reference to the Vertex   buffer object
	Visu.EBO = 0; // Reference to the Element  buffer object
	Visu.TEX = 0; // Reference to the Element  buffer object
	Visu_initOpenGL(&Visu, &Grid);
#endif




	Particles_initPhase(&Grid, &Particles);










//======================================================================================================//
//======================================================================================================//
//                                                                      				      			//
//                          				TIME LOOP             		  					        	//
//                                                                      							    //
//======================================================================================================//
//======================================================================================================//
	int timeStep = 0;
	int itNonLin;
	while(timeStep!=nTimeSteps) {

		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                          COMPUTE AND UPDATE STOKES                         //
		//                                                                            //
		//============================================================================//
		//============================================================================//
		printf("\n\n\n          ========= Time step %i =========   \n"
				     "       =====================================\n\n",timeStep);


		// Get Physics from particles to cell and to nodes
		// =================================
		TIC
		printf("Physics: Interp from particles to cell\n");
		Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes);


		printf("Physics: Interp from cell to node\n");
		Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BCStokes.SetupType);
		TOC
		printf("Physics: Interp took %.2fs\n",toc);

		// Update BC
		// =================================
		printf("BC: Update\n");
		BC_updateStokes(&BCStokes, &Grid);
		BC_updateThermal(&BCStokes, &Grid);

		// Assemble the systems of equations
		// =================================
		printf("EqStokes: Assemble\n");
		EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);
		//EqSystem_check(&EqStokes);
		EqSystem_assemble(&EqThermal, &Grid, &BCThermal, &Physics, &NumThermal);



		// Solve the heat conservation
		// =================================
		EqSystem_solve(&EqThermal, &SolverThermal, &Grid, &Physics, &BCThermal, &NumThermal);
		Physics_set_T_FromSolution(&Physics, &Grid, &BCThermal, &NumThermal, &EqThermal);


		// ============================================================================
		// 							Non-linear interation
		// ============================================================================
		itNonLin = 0;
		EqStokes.normResidual = 1.0;
		compute* NonLin_x0 = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		compute* NonLin_dx = (compute*) malloc(EqStokes.nEq * sizeof(compute));
		while(EqStokes.normResidual > nonLinTolerance && itNonLin!=maxNonLinearIter) {

			printf("\n\n  ==== Non linear iteration %i ==== \n\n",itNonLin);

			// Solve: A(X0) * X = b
			// ====================
			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_x0[iEq] = EqStokes.x[iEq];
			}
			printf("EqStokes: Solve\n");
			EqSystem_solve(&EqStokes, &SolverStokes, &Grid, &Physics, &BCStokes, &NumStokes);




			// ======================================
			// 				Line search
			// ======================================
			compute a[20];
			a[nLineSearch] = 1.0/nLineSearch;; // this is the best value
			compute minRes = 1.0;

			for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes.x[iEq] - NonLin_x0[iEq];
			}

			for (iLS= 0; iLS < nLineSearch+1; ++iLS) {
				printf("== Line search %i:  ", iLS);

				//compute a, the globalization parameter;
				if (iLS!=nLineSearch)
					a[iLS] = 1.0 - 1.0/(nLineSearch) * (iLS);

				for (iEq = 0; iEq < EqStokes.nEq; ++iEq) {
					// X1 = X0 + a*(X-X0)
					EqStokes.x[iEq] = NonLin_x0[iEq] + a[iLS]*(NonLin_dx[iEq]);
				}


				// Update the stiffness matrix
				Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BCStokes, &NumStokes, &EqStokes);

				Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BCStokes);
				Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BCStokes.SetupType);

				EqSystem_assemble(&EqStokes, &Grid, &BCStokes, &Physics, &NumStokes);



				// compute the norm of the  residual:
				// F = b - A(X1) * X1
				EqSystem_computeNormResidual(&EqStokes);


				printf("a = %.2f, |F| / |b|: %.3e, minRes = %.3e, best a = %.2f\n", a[iLS], EqStokes.normResidual, minRes, a[nLineSearch]);
				if (EqStokes.normResidual<minRes) {
					a[nLineSearch] = a[iLS];
					minRes = EqStokes.normResidual;
					if (iLS==nLineSearch-1) // if the last one is the best one then don't recompute, i.e. next one would be the same
						break;
				}

				if (itNonLin==0)
					break;

			}
			// ======================================
			// 		   End of Line search
			// ======================================


			itNonLin++;
		} // end of non-linear loop

		free(NonLin_x0);
		free(NonLin_dx);

		// ============================================================================
		// 						End of Non-linear interation
		// ============================================================================












		// ============================================================================
		// 									Advect
		// ============================================================================
		// Update dt

		if (fabs(Physics.maxV)<1E-6)
			Physics.maxV = 1E-6;

			Physics.dt = CFL_fac*fmin(Grid.dx,Grid.dy)/(Physics.maxV); // note: the min(dx,dy) is the char length, so = 1
		if (Physics.dt>0.5) {
			Physics.dt = 0.5;
		}
		printf("maxV = %.3em Physics.dt = %.3e\n",fabs(Physics.maxV), Physics.dt);

		printf("Particles: Advect\n");
		Particles_advect(&Particles, &Grid, &Physics);
		//printf("Grid: Update pure shear\n");
		switch (BCStokes.SetupType) {
		case PureShear:
			Grid_updatePureShear(&Grid, &BCStokes, Physics.dt);
			break;
		case SimpleShearPeriodic:
			Particles_Periodicize(&Grid, &Particles, &BCStokes);
			break;
		default:
			break;
		}
		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(&Grid, &Particles);
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
















#if VISU
		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                                 VISUALIZATION                              //
		//                                                                            //
		//============================================================================//
		//============================================================================//
		do  {
			//printf("A-2\n");
				glfwPollEvents();
				///printf("A-3\n");
				Visu_checkInput(&Visu, window);
				//printf("A-4\n");
				glClearColor(0, 0, 0, 1); // black
				//glClear(GL_COLOR_BUFFER_BIT);
				//printf("A-1\n");
				glEnable(GL_DEPTH_TEST);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				//glClear(GL_COLOR_BUFFER_BIT);

				if (Visu.initPassivePart) {
					Particles_initPassive(&Grid, &Particles);
					Visu.initPassivePart = false;
				}


				//============================================================================
				// 								PLOT GRID DATA

				// ****** Bind shader textures, arrays and buffers
				glBindVertexArray(Visu.VAO);
				glUseProgram(Visu.ShaderProgram);
				glBindTexture(GL_TEXTURE_2D, Visu.TEX);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu.EBO);

					// 1. Update data
					if (BCStokes.SetupType==PureShear) {
						Visu_updateVertices(&Visu, &Grid);
						glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO);
								glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu.vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}
					Visu_update(&Visu, window, &Grid, &Physics, &BCStokes, &Char);
					// update the content of Visu.U
					glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, Grid.nxS, Grid.nyS, 0, GL_RED, GL_FLOAT, Visu.U);	// load the updated Visu.U in the texture

					// 2. Draw
					glDrawElements(GL_TRIANGLES, Visu.ntrivert, GL_UNSIGNED_INT, 0);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
				glBindTexture(GL_TEXTURE_2D, 0);
				glUseProgram(0);
				glBindVertexArray(0);
				// ****** Unbind textures, arrays and buffers


				// 								PLOT GRID DATA
				//============================================================================



				//============================================================================
				// 								PLOT PARTICLE
				if (Visu.showParticles) {
					//printf("A\n");
					glBindVertexArray(Visu.VAO_part);
					glUseProgram(Visu.ParticleShaderProgram);


					glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_part);

					//glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_partMesh);
					// update the buffer containing the particles
						Visu_particles(&Visu, &Particles, &Grid);
						//Visu_particleMesh(&Visu);
						Visu_updateUniforms(&Visu, window);
						//printf("B\n");
						//glBufferData(GL_ARRAY_BUFFER, 4*Visu.nParticles*sizeof(GLfloat), NULL, GL_STREAM_DRAW);
						glBufferSubData(GL_ARRAY_BUFFER, 0, 4*Particles.n*sizeof(GLfloat), Visu.particles);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
		//				glDrawArrays(GL_POINTS, 0, Particles.n);
						glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO_partMesh);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, Visu.particleMeshRes+2, Visu.nParticles);
						//printf("Visu.particleMeshRes= %i\n",Visu.particleMeshRes);
						//glDrawArraysInstanced(GL_TRIANGLES, 0, 3, Particles.n);
						//
					glUseProgram(0);
					glBindVertexArray(0);
				}
				// 								PLOT PARTICLE
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
#if VISU
	// Quit glfw
	glfwDestroyWindow(window);
	glfwTerminate();
	Visu_freeMemory(&Visu);
#endif

	printf("SUCCESS!\n\n\n");

	return EXIT_SUCCESS;
}

