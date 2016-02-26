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
	Grid.nxC = 64;
	Grid.nyC = 64;

	Particles.nPCX = 2;
	Particles.nPCY = 2;

	Grid.xmin = -1;
	Grid.xmax =  1;
	Grid.ymin = -1;
	Grid.ymax =  1;

	MatProps.nPhase  = 2;
	MatProps.rho0[0] = 1; 	MatProps.eta0[0] = 100;
	MatProps.rho0[1] = 1;	MatProps.eta0[1] = 1;




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



	// Get Physics from particles to cell and to nodes
	// =================================
	printf("Physics: Interp from particles to cell\n");
	Physics_interpFromParticlesToCell(&Grid, &Particles, &Physics, &MatProps);
	printf("Physics: Interp from cell to node\n");
	Physics_interpFromCellToNode(&Grid, Physics.eta, Physics.etaShear);



	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_set(&BC, &Grid, &EqSystem, &Physics);;



	// Initialize Numbering maps without dirichlet and EqSystem->I
	// =================================
	EqSystem_allocateI(&EqSystem);
	Numbering_initMapAndSparseTripletIJ(&BC, &Grid, &EqSystem, &Numbering);



	// Allocate memory for the system of equations
	// =================================
	EqSystem_allocateMemory(&EqSystem);



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



	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          COMPUTE AND UPDATE STOKES                         //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	// Compute Physics variable on the base grid
	// based on the phase of particles
	// =================================
	printf("Particles: Get Physics from\n");
	//Particles_getPhysicsFrom(&Grid, &Particles, &Physics, &MatProps);

	printf("Particles Update Linked List\n");
	Particles_updateLinkedList(&Grid, &Particles);


	printf("Finished Updating Linked List\n");





	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                 VISUALIZATION                              //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	if (DEBUG) {
		printf("=== Check eta ===\n");
		int C = 0;
		int ix, iy;
		for (iy = 0; iy < Grid.nyC; ++iy) {
			for (ix = 0; ix < Grid.nxC; ++ix) {
				printf("%.2f  ", Physics.eta[C]);
				C++;
			}
			printf("\n");
		}
	}



	if (VISU) {
		//Visu_updateCenterValue(&Visu, &Grid, Physics.eta);

		Visu_StrainRate(&Visu, &Grid, &Physics);

		while(!glfwWindowShouldClose(window)){
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

