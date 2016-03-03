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



	// Set model properties
	// =================================
	Grid.nxC = 16;
	Grid.nyC = 16;

	Particles.nPCX = 3;
	Particles.nPCY = 3;

	Grid.xmin = -1;
	Grid.xmax =  1;
	Grid.ymin = -1;
	Grid.ymax =  1;

	MatProps.nPhase  = 2;
	MatProps.rho0[0] = 1.0; 	MatProps.eta0[0] = 1.0;
	MatProps.rho0[1] = 1.0;		MatProps.eta0[1] = 0.001;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

	BC.SetupType = 0;
	BC.backStrainRate = 1.0;
	BC.VxB =  1.0;	BC.VyB = 0.0;
	BC.VxT = -1.0;	BC.VyT = 0.0;

	//Physics.dt = Grid.dx/2*BC.VxL;
	Physics.dt = Grid.dx/(2*abs(BC.VxT-BC.VxB));

	Physics.g[0] = 0;
	Physics.g[1] = 9.81;

	// Set characteristic quantities
	// =================================
	Char.length 		= fmin(Grid.dx,Grid.dy);
	Char.density 		= MatProps.rho0[0];
	Char.acceleration 	= Physics.g[1];
	Char.viscosity  	= 1.0;//pow( 10, (log10(MatProps.eta0[0])+log10(MatProps.eta0[0]))/2 );
	Char.stress 		= Char.density*Char.acceleration*Char.length; // i.e. rho*g*h

	Char.time 			= Char.viscosity/Char.stress;
	Char.velocity 		= Char.length/Char.time;
	Char.strainrate 	= 1.0/Char.time;

	Visu.type = Viscosity; // Default

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(&Char, &Grid, &Physics, &MatProps, &BC);

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
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps);
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






	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          	INIT VISUALIZATION                            //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	GLFWwindow* window = NULL;

	if (VISU) {
		// Init GLFW
		// =======================================

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
	}





	if (VISU) {

		int timeStep = 0;
		int nsteps = 1;
		int istep;
		//for (istep = 0; istep < nsteps; ++istep) {
		while(!glfwWindowShouldClose(window)){

			//============================================================================//
			//============================================================================//
			//                                                                            //
			//                          COMPUTE AND UPDATE STOKES                         //
			//                                                                            //
			//============================================================================//
			//============================================================================//
			printf("\n\n==== Time step %i ====\n",timeStep);

			// Get Physics from particles to cell and to nodes
			// =================================
			printf("Physics: Interp from particles to cell\n");
			Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps);
			printf("Physics: Interp from cell to node\n");
			Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear, BC.SetupType);

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



			printf("=== Check eta cell ===\n");
			int C = 0;
			int ix, iy;
			for (iy = 0; iy < Grid.nyC; ++iy) {
				for (ix = 0; ix < Grid.nxC; ++ix) {
					printf("%.3f  ", Physics.eta[C]);
					C++;
				}
				printf("\n");
			}

			printf("=== Check eta node ===\n");
			C = 0;
			//int ix, iy;
			for (iy = 0; iy < Grid.nyS; ++iy) {
				for (ix = 0; ix < Grid.nxS; ++ix) {
					printf("%.3f  ", Physics.etaShear[C]);
					C++;
				}
				printf("\n");
			}

			int iCell = 0;
			//int ix, iy;
			printf("=== Linked list heads ===\n");
			for (iy = 0; iy < Grid.nyC; ++iy) {
				for (ix = 0; ix < Grid.nxC; ++ix) {
					printf("%i  ", Particles.linkHead[iCell]);
					iCell++;
				}
				printf("\n");
			}


			// Solve
			// =================================
			printf("EqSystem: Solve\n");
			EqSystem_solve(&EqSystem);

			// Reconstruct Vx, Vy, P from the solution vector
			// =================================
			printf("Physics: SetVx Vy from Sol\n");
			Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BC, &Numbering, EqSystem.x);

			// Update dt
			// =================================
			printf("maxV = %.2f\n",Physics.maxV);
			Physics.dt = 0.9*fmin(Grid.dx,Grid.dy)/(Physics.maxV); // note: the min(dx,dy) is the char length, so = 1

			// Advect particles and update grid
			// =================================
			if (ADVECT) {
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
			}



			// Update the linked list of particles
			// =================================
			printf("Particles Update Linked List\n");
			Particles_updateLinkedList(&Grid, &Particles);



			printf("minCoeff = %.2f, maxCoeff = %.2f\n",absmin(EqSystem.V,EqSystem.nnz), absmax(EqSystem.V,EqSystem.nnz));


			printf("xmin = %.2f, xmax = %.2f, ymin = %.2f, ymax = %.2f\n",Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax);









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

			// Update U data
			glBindBuffer(GL_ARRAY_BUFFER, Visu.VBO);
				Visu_updateVertices(&Visu, &Grid);

				glBufferSubData(GL_ARRAY_BUFFER, 0, 2*Grid.nxS*Grid.nyS*sizeof(GLfloat), Visu.vertices);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			glBindBuffer(GL_ARRAY_BUFFER, Visu.CBO);
				Visu_update(&Visu, window, &Grid, &Physics, &BC, &Char);
				//printf("minEps = %.2f, maxEps = %.2f\n",minf(Visu.U,Grid.nSTot), maxf(Visu.U, Grid.nSTot));
				printf("minEps = %.2f, maxEps = %.2f\n",minf(Visu.U,Grid.nSTot)/BC.backStrainRate, maxf(Visu.U, Grid.nSTot)/BC.backStrainRate);
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
	EqSystem_freeMemory(&EqSystem);
	if (VISU) {
		Visu_freeMemory(&Visu);
	}

	printf("SUCCESS!");

	return EXIT_SUCCESS;
}

