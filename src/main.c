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
	printf("Beginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());
	// Declare structures
	// =================================
	Grid 		Grid;
	MatProps 	MatProps;
	Particles 	Particles;
	Physics 	Physics;
	Visu 		Visu;
	EqSystem 	EqSystem;
	BC 			BC;
	Numbering 	Numbering;
	Char 		Char;
	Solver 		Solver;

	INIT_TIMER

	// Set model properties
	// =================================
	int nTimeSteps = 1; //  negative value for infinite
	int nLineSearch = 2;
	int maxNonLinearIter = 5;
	compute nonLinTolerance = 5E-3;

	Grid.nxC = 3;
	Grid.nyC = 4;

	Particles.nPCX = 3;
	Particles.nPCY = 3;

	//Grid.xmin = 0;
	//Grid.xmax = (compute) Grid.nxC;
	//Grid.ymin = 0;
	//Grid.ymax = (compute) Grid.nyC;
	Grid.xmin = -1.0;
	Grid.xmax =  1.0;
	Grid.ymin = -1.0;
	Grid.ymax =  1.0;

	MatProps.nPhase  = 2;
	MatProps.rho0[0] = 0.9; 		MatProps.eta0[0] = 1.0;  		MatProps.n[0] = 3.0;
	MatProps.rho0[1] = 1;		MatProps.eta0[1] = 0.9; 		MatProps.n[1] = 3.0;

	MatProps.flowLaw[0] = LinearViscous;
	MatProps.flowLaw[1] = LinearViscous;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	BC.SetupType = 0;
	BC.backStrainRate = -0;
	BC.VxB =  1.0;	BC.VyB = 0.0;
	BC.VxT = -1.0;	BC.VyT = 0.0;

	//Physics.dt = Grid.dx/2*BC.VxL;
	Physics.dt = Grid.dx/(2*abs(BC.VxT-BC.VxB));
	Physics.epsRef = 1.0;//abs(BC.backStrainRate);

	Physics.g[0] = 0;
	Physics.g[1] = -9.81;

	compute CFL_fac = 2.0; // 0.5 ensures stability


	EqSystem.penaltyMethod = false;
	EqSystem.penaltyFac = 1000000;


	// Set characteristic quantities
	// =================================
	Char.length 		= fmin(Grid.dx,Grid.dy);
	Char.density 		= MatProps.rho0[0];
	Char.acceleration 	= Physics.g[1];
	Char.viscosity  	= MatProps.eta0[0];//pow( 10, (log10(MatProps.eta0[0])+log10(MatProps.eta0[0]))/2 );

	Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h
	//Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h

	Char.time 			= Char.viscosity/Char.stress;
	Char.velocity 		= Char.length/Char.time;
	Char.strainrate 	= 1.0/Char.time;

	Visu.type 			= StrainRate; // Default

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BC);
	Physics.etaMin = 1E-4;
	Physics.etaMax = 1E4;

	// Init grid and particles
	// =================================
	Grid.nCTot  = Grid.nxC*Grid.nyC;


	Grid.nxVx 	= Grid.nxC+1; 		Grid.nyVx	= Grid.nyC+2;
	Grid.nxVy 	= Grid.nxC+2;		Grid.nyVy	= Grid.nyC+1;
	Grid.nxS 	= Grid.nxC+1;		Grid.nyS	= Grid.nyC+1;
	Grid.nSTot  = Grid.nxS*Grid.nyS;

	Grid.nVxTot = Grid.nxVx*Grid.nyVx;
	Grid.nVyTot = Grid.nxVy*Grid.nyVy;


	EqSystem.nEqIni  = Grid.nxVx*Grid.nyVx + Grid.nxVy*Grid.nyVy + Grid.nxC*Grid.nyC;

	Particles.nPC 	= Particles.nPCX * Particles.nPCY;
	Particles.n 	= Grid.nCTot*Particles.nPC;

	Visu.ntri   	= Grid.nxC*Grid.nyC*2;
	Visu.ntrivert 	= Visu.ntri*3;


	// Other variables
	// =================================
	int iEq, iLS;

	// Allocate memory
	// =================================
	printf("Allocate memory\n");
	Memory_allocateMain(&Grid, &Particles, &Physics, &EqSystem, &Numbering);



	// Initialize Particles' coordinates
	// =================================
	printf("Particles: Init Coord\n");
	Particles_initCoord(&Grid, &Particles);



	// Initialize Particles' phase
	// =================================
	printf("Particles: Init Phase\n");
	Particles_initPhase(&Grid, &Particles);



	// Get Physics from particles to cell and to nodes (important for Neumann conditions)
	// =================================
	printf("Physics: Interp from particles to cell\n");
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BC);
	printf("Physics: Interp from cell to node\n");
	Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BC.SetupType);


	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_set(&BC, &Grid, &EqSystem, &Physics);



	// Initialize Numbering maps without dirichlet and EqSystem->I
	// =================================
	EqSystem_allocateI(&EqSystem);
	Numbering_initMapAndSparseTripletIJ(&BC, &Grid, &EqSystem, &Numbering);



	// Allocate memory for the system of equations
	// =================================
	EqSystem_allocateMemory(&EqSystem);


	// Init Solver
	// =================================
	printf("EqSystem: Init Solver\n");
	EqSystem_assemble(&EqSystem, &Grid, &BC, &Physics, &Numbering); // dummy assembly to give the EqSystem initSolvers
	//EqSystem_check(&EqSystem);
	EqSystem_initSolver (&EqSystem, &Solver);


	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          	INIT VISUALIZATION                            //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	if (VISU) {
		// Init GLFW
		// =======================================
		printf("Init Window\n");
			GLFWwindow* window = NULL;
			printf("Done\n");

		Visu_initWindow(&window);

		Visu_allocateMemory(&Visu, &Grid);
		Visu_init(&Visu, &Grid);


		///Init shader
		// =======================================
		Visu.VertexShaderFile = "src/shader.vs";
		Visu.FragmentShaderFile = "src/shader.fs";
		Visu.ShaderProgram = 0;
		// Generate reference to objects (indexes that act as pointers to graphic memory)
		// =======================================
		Visu.VAO = 0; // Reference to the Vertex   array object
		Visu.VBO = 0; // Reference to the Vertex   buffer object
		Visu.CBO = 0; // Reference to the Color    buffer object
		Visu.EBO = 0; // Reference to the Element  buffer object

		Visu_initOpenGL(&Visu, &Grid);
	//}





//	if (VISU) {

		int timeStep = 0;
		int itNonLin;
		while(!glfwWindowShouldClose(window) && timeStep!=nTimeSteps){

			//nLineSearch++;
			//maxNonLinearIter += 1;

			//============================================================================//
			//============================================================================//
			//                                                                            //
			//                          COMPUTE AND UPDATE STOKES                         //
			//                                                                            //
			//============================================================================//
			//============================================================================//
			printf("\n\n==== Time step %i ====\n",timeStep);
			int i;
			for (i = 0; i < Grid.nVxTot; ++i) {
				Physics.Vx[i] = 0;
			}
			for (i = 0; i < Grid.nVyTot; ++i) {
				Physics.Vy[i] = 0;
			}
			for (iEq = 0; iEq < EqSystem.nEq; ++iEq) {
				EqSystem.x[iEq] = 0;
			}
			itNonLin = 0;
			EqSystem.normResidual = 1.0;
			while(EqSystem.normResidual > nonLinTolerance && itNonLin!=maxNonLinearIter) {
				printf("==== Non linear iteration %i \n",itNonLin);



				// Get Physics from particles to cell and to nodes
				// =================================
				TIC
				printf("Physics: Interp from particles to cell\n");
				Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BC);
				printf("Physics: Interp from cell to node\n");
				Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BC.SetupType);
				TOC
				printf("Physics: Interp took %.2fs\n",toc);
				// Update BC
				// =================================
				printf("BC: Update\n");
				BC_updateDir(&BC, &Grid);
				BC_updateNeuCoeff(&BC, &Grid, &Physics);

				// Assemble the system of equations
				// =================================
				printf("EqSystem: Assemble\n");
				EqSystem_assemble(&EqSystem, &Grid, &BC, &Physics, &Numbering);


				//EqSystem_check(&EqSystem);

				// Solve
				// =================================
				for (iEq = 0; iEq < EqSystem.nEq; ++iEq) {
					EqSystem.x0[iEq] = EqSystem.x[iEq];
				}
				printf("EqSystem: Solve\n");
				//EqSystem_solve(&EqSystem, &Solver, &Grid, &Physics, &BC, &Numbering);

				compute a[20];

				a[0] = 1.0/nLineSearch;
				a[nLineSearch+1] = 1.0/nLineSearch;; // this is the best value
				compute minRes = 1.0;

				for (iEq = 0; iEq < EqSystem.nEq; ++iEq) {
					EqSystem.dx[iEq] = EqSystem.x[iEq] - EqSystem.x0[iEq];
				}

				for (iLS= 0; iLS < nLineSearch+1; ++iLS) {
					printf("== Line search %i\n", iLS);


					// In the previous solve step the solution to the following system of equation was computed:
					// A(X0) * X = b

					// We now update the solution in the following manner:
					// X1 = X0 + a*(X-X0)
					// where a is a globalization factor

					//compute a;
					if (iLS!=nLineSearch)
						a[iLS] = 1.0 - 1.0/nLineSearch * (iLS);
					for (iEq = 0; iEq < EqSystem.nEq; ++iEq) {
						EqSystem.x[iEq] = EqSystem.x0[iEq] + a[iLS]*(EqSystem.dx[iEq]);
					}
					// Compute A(X1), /!\ might not be valid for penalty method in the current state

					Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BC, &Numbering, &EqSystem);

					Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps, &BC);
					Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BC.SetupType);

					BC_updateNeuCoeff(&BC, &Grid, &Physics);

					EqSystem_assemble(&EqSystem, &Grid, &BC, &Physics, &Numbering);

					// compute the norm of the  residual:
					// F = b - A(X1) * X1
					EqSystem_computeNormResidual(&EqSystem);

					// compute the norm of b
					compute norm_b = 0;
					for (iEq = 0; iEq < EqSystem.nEq; ++iEq) {
						norm_b += EqSystem.b[iEq]*EqSystem.b[iEq];
					}
					norm_b = sqrt(norm_b);

					EqSystem.normResidual /= norm_b; // Normalize the norm of the residual by the norm of the right hand side

					//printf("x[0]=%.3e, x[25]=%.3e\n", EqSystem.x[0], EqSystem.x[25]);
					//printf("max dx=%.3e\n", max(EqSystem.dx, EqSystem.nEq));

					printf("a = %.2f, |F| / |b|: %.3e, minRes = %.3e, best a = %.2f\n", a[iLS], EqSystem.normResidual, minRes, a[nLineSearch]);
					if (EqSystem.normResidual<minRes) {
						a[nLineSearch] = a[iLS];
						minRes = EqSystem.normResidual;
						if (iLS==nLineSearch-1) // if the last one is the best one then don't recompute, i.e. next one would be the same
							break;
					}
					if (itNonLin==0) {
						break;
					}


				}

				itNonLin++;
			} // end of non-linear loop






			// Update dt
			// =================================
			//printf("maxV = %.2f\n",Physics.maxV);
			Physics.dt = CFL_fac*fmin(Grid.dx,Grid.dy)/(Physics.maxV); // note: the min(dx,dy) is the char length, so = 1

			// Advect particles and update grid
			// =================================
			printf("Particles: Advect\n");
			Particles_advect(&Particles, &Grid, &Physics);
			printf("Grid: Update pure shear\n");
			switch (BC.SetupType) {
			case 0:
				Grid_updatePureShear(&Grid, &BC, Physics.dt);
				break;
			case 1:
				Particles_Periodicize(&Grid, &Particles, &BC);
				break;
			default:
				break;
			}



			// Update the linked list of particles
			// =================================
			printf("Particles Update Linked List\n");
			Particles_updateLinkedList(&Grid, &Particles);



			//printf("minCoeff = %.2f, maxCoeff = %.2f\n",absmin(EqSystem.V,EqSystem.I[EqSystem.nRow-1]), absmax(EqSystem.V,EqSystem.I[EqSystem.nRow-1]));
			/*
			printf("===== RHS =====\n");
			int i;
			for (i=0; i<EqSystem.nEq; i++) {
				printf("RHS[%i] = %.2f\n", i, EqSystem.b[i]);
			}
			 */
			//printf("\n");
			//printf("minRHS = %.2f, maxRHS = %.2f\n",min(EqSystem.b,EqSystem.nEq), max(EqSystem.b,EqSystem.nEq));

			//printf("xmin = %.2f, xmax = %.2f, ymin = %.2f, ymax = %.2f\n",Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax);









			//============================================================================//
			//============================================================================//
			//                                                                            //
			//                                 VISUALIZATION                              //
			//                                                                            //
			//============================================================================//
			//============================================================================//


			// process pending events
			glfwPollEvents();
			// clear everything
			glClearColor(0, 0, 0, 1); // black
			glClear(GL_COLOR_BUFFER_BIT);

			// bind the program (the shaders)
			glUseProgram(Visu.ShaderProgram);

			// Update vertices
			glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO);
			Visu_updateVertices(&Visu, &Grid);

			glBufferSubData(GL_ARRAY_BUFFER, 0, 2*Grid.nxS*Grid.nyS*sizeof(GLfloat), Visu.vertices);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			// update nodal data1
			glBindBuffer(GL_ARRAY_BUFFER, Visu.CBO);
			Visu_update(&Visu, window, &Grid, &Physics, &BC, &Char);
			//printf("minU = %.2f, maxU = %.2f\n",minf(Visu.U,Grid.nSTot)/Visu.valueScale, maxf(Visu.U, Grid.nSTot)/Visu.valueScale);
			glBufferSubData(GL_ARRAY_BUFFER, 0, Grid.nxS*Grid.nyS*sizeof(GLfloat), Visu.U);
			glBindBuffer(GL_ARRAY_BUFFER, 0);


			// bind the VAO (the triangle)
			glBindVertexArray(Visu.VAO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu.EBO);

			// draw the VAO
			glDrawElements(GL_TRIANGLES, Visu.ntrivert, GL_UNSIGNED_INT, 0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

			// unbind the VAO
			glBindVertexArray(0);


			// unbind the program
			glUseProgram(0);

			//Swap windows
			glfwSwapBuffers(window);
			timeStep++;


		}

		// Quit glfw
		glfwDestroyWindow(window);
		glfwTerminate();



	}

	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	// Free memory
	Memory_freeMain(&Particles, &Physics, &Numbering);
	EqSystem_freeMemory(&EqSystem, &Solver);
	if (VISU) {
		Visu_freeMemory(&Visu);
	}

	printf("SUCCESS!\n");

	return EXIT_SUCCESS;
}

