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
	Grid Grid;
	MatProps MatProps;
	Particles Particles;
	Physics Physics;
	Visu Visu;




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

	MatProps.nPhase = 2;
	MatProps.rho0[0] = 1; 	MatProps.eta0[0] = 0;
	MatProps.rho0[1] = 1;	MatProps.eta0[1] = 1;




	// Init grid and particles
	// =================================
	Grid.nCTot  = Grid.nxC*Grid.nyC;

	Grid.nxVx 	= Grid.nxC; 		Grid.nyVx	= Grid.nyC+1;
	Grid.nxVy 	= Grid.nxC+1;		Grid.nyVy	= Grid.nyC;
	Grid.nxS 	= Grid.nxC+1;		Grid.nyS	= Grid.nyC+1;
	Grid.nxN 	= Grid.nxC;			Grid.nyN	= Grid.nyC;

	Grid.nVxTot = Grid.nxVx*Grid.nyVx;
	Grid.nVyTot = Grid.nxVy*Grid.nyVy;

	Grid.dx = (Grid.xmax-Grid.xmin)/Grid.nxC;
	Grid.dy = (Grid.ymax-Grid.ymin)/Grid.nyC;


	Particles.nPC 	= Particles.nPCX * Particles.nPCY;
	Particles.n 	= Grid.nCTot*Particles.nPC;


	Visu.ntri   	= Grid.nxC*Grid.nyC*2;
	Visu.ntrivert 	= Visu.ntri*3;



	// Allocate memory
	// =================================
	printf("Allocate memory\n");
	allocateMemory(&Grid, &Particles,&Physics);



	// Initialize Particles' coordinates
	// =================================
	printf("Particles: Init Coord\n");
	Particles_initCoord(&Grid, &Particles);


	// Initialize Particles' phase
	// =================================
	Particles_initPhase(&Grid, &Particles);











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
	Particles_getPhysicsFrom(&Grid, &Particles, &Physics, &MatProps);

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


	Visu_updateCenterValue(&Visu, &Grid, Physics.eta);
	if (VISU) {
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
	freeMemory(&Particles,&Physics);
	Visu_freeMemory(&Visu);



	printf("SUCCESS!");
	return EXIT_SUCCESS;
}

