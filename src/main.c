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



	// Set model properties
	// =================================
	Grid.nxC = 128;
	Grid.nyC = 128;

	Particles.nPCX = 2;
	Particles.nPCY = 2;

	Grid.xmin = -1;
	Grid.xmax =  1;
	Grid.ymin = -1;
	Grid.ymax =  1;

	MatProps.nPhase  = 2;
	MatProps.rho0[0] = 1.0; 	MatProps.eta0[0] = 1.0;
	MatProps.rho0[1] = 1.0;		MatProps.eta0[1] = 100.0;

	Physics.dt = 0.01;

	BC.VxL = -1.0*Grid.xmin; BC.VxR = -1.0*Grid.xmax;
	BC.VyB =  1.0*Grid.ymin; BC.VyT =  1.0*Grid.ymax;

	// Init grid and particles
	// =================================
	Grid.nCTot  = Grid.nxC*Grid.nyC;

	Grid.nxVx 	= Grid.nxC+1; 		Grid.nyVx	= Grid.nyC+2;
	Grid.nxVy 	= Grid.nxC+2;		Grid.nyVy	= Grid.nyC+1;
	Grid.nxS 	= Grid.nxC+1;		Grid.nyS	= Grid.nyC+1;
	//Grid.nxN 	= Grid.nxC;			Grid.nyN	= Grid.nyC;

	Grid.nVxTot = Grid.nxVx*Grid.nyVx;
	Grid.nVyTot = Grid.nxVy*Grid.nyVy;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;

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
	Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);


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
		for (istep = 0; istep < nsteps; ++istep) {
		//while(!glfwWindowShouldClose(window)){

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
			Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);

			printf("Finished Updating Linked List\n");

			// Assemble the system of equations
			// =================================
			EqSystem_assemble(&EqSystem, &Grid, &BC, &Physics, &Numbering);
			//EqSystem_check(&EqSystem);


			// Solve
			// =================================
			EqSystem_solve(&EqSystem);

			// Reconstruct Vx, Vy, P from the solution vector
			// =================================
			Physics_set_VxVyP_FromSolution(&Physics, &Grid, &BC, &Numbering, EqSystem.x);

			// Advect particles and update grid
			// =================================
			Particles_advect(&Particles, &Grid, &Physics);
			//Grid_updatePureShear(&Grid, &BC, Physics.dt);

			// update boundary conditions
			// =================================
			printf("BC: Update\n");
			BC.VxL = -1.0*Grid.xmin; BC.VxR = -1.0*Grid.xmax;
			BC.VyB =  1.0*Grid.ymin; BC.VyT =  1.0*Grid.ymax;
			BC_updateDir(&BC, &Grid);

			// Update the linked list
			// =================================
			printf("Particles Update Linked List\n");
			//Particles_updateLinkedList(&Grid, &Particles);




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
			glBindBuffer(GL_ARRAY_BUFFER, Visu.CBO);
			//Visu_updateCenterValue(&Visu, &Grid, Physics.eta);
			Visu_updateVertices(&Visu, &Grid);
			Visu_StrainRate(&Visu, &Grid, &Physics);
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
			printf("xmin = %.2f, xmax = %.2f, ymin = %.2f, ymax = %.2f",Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax);
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

