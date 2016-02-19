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


static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

int main(void) {

	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                              INIT STOKES PROGRAM                           //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	// Declare structures
	// =================================
	Grid Grid;
	MatProps MatProps;
	Particles Particles;
	Physics Physics;
	Visu Visu;




	// Set model properties
	// =================================
	Grid.nxC = 3;
	Grid.nyC = 4;

	Particles.nPCX = 2;
	Particles.nPCY = 2;

	Grid.xmin = 0;
	Grid.xmax = 3;
	Grid.ymin = 0;
	Grid.ymax = 4;

	MatProps.nPhase = 2;
	MatProps.rho0[0] = 1; 	MatProps.eta0[0] = 1;
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
	allocateMemory(&Grid, &Particles,&Physics);


	// Initialize Particle coordinates
	// =================================
	initParticlesCoord(&Grid, &Particles);


	// Compute Physics variable on the base grid
	// based on the phase of particles
	// =================================
	getPhysicsFromParticles2Grid(&Grid, &Particles, &Physics, &MatProps);








	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                           INIT WINDOW AND OPENGL                           //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	/// Init GLFW
	// =======================================
	GLFWwindow* window;
	glfwSetErrorCallback(error_callback);
	if (!glfwInit()){
		exit(EXIT_FAILURE);
	}
	#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	#endif
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);


	/// Create window
	// =======================================
	window = glfwCreateWindow(1024, 1024, "Simple example", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);

	/// Init Glew - Must be done after glut is initialized!
	// =======================================
	glewExperimental = GL_TRUE;
	GLenum res = glewInit();
	if (res != GLEW_OK) {
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
	}
	if(!GLEW_VERSION_3_2){
		fprintf(stderr, "OpenGL 3.2 API is not available.");
		return 1;
	}

	/// Test GL version
	// =======================================
	const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
	const GLubyte* version = glGetString (GL_VERSION); // version as a string
	const GLubyte* glslversion = glGetString (GL_SHADING_LANGUAGE_VERSION); // version as a string

	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);
	printf("GLSL version supported %s\n", glslversion);



	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          	INIT VISUALIZATION                            //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	allocateVisuMemory(&Visu);
	compute* U = Physics.rho;


	initVisualization(&Visu, &Grid);


    // Loop

		// Compute



		// PlotScalar = pointer to something

		// Plot

	// end loop













	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          COMPUTE AND UPDATE STOKES                         //
	//                                                                            //
	//============================================================================//
	//============================================================================//










	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                 VISUALIZATION                              //
	//                                                                            //
	//============================================================================//
	//============================================================================//









	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	// Free memory
	freeMemory(&Particles,&Physics);
	freeVisuMemory(&Visu);

	printf("SUCCESS!");
	return EXIT_SUCCESS;
}

