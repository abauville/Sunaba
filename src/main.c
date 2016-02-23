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
	Grid.nxC = 256;
	Grid.nyC = 256;

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

	Visu_allocateMemory(&Visu, &Grid);
	//compute* U = Physics.rho;


	Visu_init(&Visu, &Grid);
	Visu_plotCenterValue(&Visu, &Grid, Physics.eta);



	if (VISU) {

		///Init shader
		// =======================================
		const char* pVSFileName = "src/shader.vs";
		const char* pFSFileName = "src/shader.fs";
		GLuint ShaderProgram = 0;
		// Generate reference to objects (indexes that act as pointers to graphic memory)
		// =======================================
		GLuint VAO = 0; // Reference to the Vertex   array object
		GLuint VBO = 0; // Reference to the Vertex   buffer object
		GLuint CBO = 0; // Reference to the Color    buffer object
		GLuint EBO = 0; // Reference to the Element  buffer object

		// And assigned them to objects (stored in the graphic memory)
		// =======================================
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &CBO);
		glGenBuffers(1, &EBO);

		// Bind Vertex Array object
		// =======================================
		glBindVertexArray(VAO);
		// compile shaders
		// =======================================
		compileShaders(&ShaderProgram, pVSFileName, pFSFileName);

		glUseProgram(ShaderProgram);

		// Get IDs for the in attributes of the shader
		// =======================================
		GLint VertAttrib    = glGetAttribLocation(ShaderProgram,"in_Vertex");
		GLint SolAttrib     = glGetAttribLocation(ShaderProgram,"U");
		// Bind objects and associate with data tables
		// =======================================
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, Grid.nxS*Grid.nyS*sizeof(coord), Visu.vertices, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, CBO);
		glBufferData(GL_ARRAY_BUFFER, Grid.nxS*Grid.nyS*sizeof(GLfloat), Visu.U, GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, Visu.ntrivert*sizeof( GLuint ), Visu.elements, GL_STATIC_DRAW);

		// Connect Vertex data (stored in VBO) to the "in_Vertex" attribute of the shader
		// =======================================
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(VertAttrib, 2, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(VertAttrib);

		// Connect Color data (stored in CBO) to the "in_Color" attribute of the shader
		// =======================================
		glBindBuffer(GL_ARRAY_BUFFER, CBO);
		glVertexAttribPointer(SolAttrib, 1, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(SolAttrib);



		// Declare the scale as a uniform
		// =======================================
		GLfloat Scale;

		if ((Grid.xmax-Grid.xmin)>(Grid.ymax-Grid.ymin)){
			Scale = 2.0/(1.1*(Grid.xmax-Grid.xmin));
		}
		else {
			Scale = 2.0/(1.1*(Grid.ymax-Grid.ymin));
		}
		//Scale = 1.0;
		GLfloat Transform[] = {Scale,0.0f,0.0f,0.0f , 0.0f,Scale,0.0f,0.0f , 0.0f,0.0f,1.0f,0.0f , 0.0f,0.0f,0.0f,1.0f};
		GLuint transformLoc = glGetUniformLocation(ShaderProgram, "transform");
		glUniformMatrix4fv(transformLoc, 1, GL_FALSE, &Transform[0]);


		// unbind the Buffer object (VBO, CBO) and Vertex array object (VAO)
		// =======================================
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
		glUseProgram(0);












		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                          COMPUTE AND UPDATE STOKES                         //
		//                                                                            //
		//============================================================================//
		//============================================================================//

		//Visu_plotCenterValue(&Visu, &Grid, Physics.eta);




		//============================================================================//
		//============================================================================//
		//                                                                            //
		//                                 VISUALIZATION                              //
		//                                                                            //
		//============================================================================//
		//============================================================================//

		while(!glfwWindowShouldClose(window)){
			// process pending events
			glfwPollEvents();
			// clear everything
			glClearColor(0, 0, 0, 1); // black
			glClear(GL_COLOR_BUFFER_BIT);

			// bind the program (the shaders)
			glUseProgram(ShaderProgram);
			// bind the VAO (the triangle)
			glBindVertexArray(VAO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
			// draw the VAO
			glDrawElements(GL_TRIANGLES, Visu.ntrivert, GL_UNSIGNED_INT, 0);
			// unbind the VAO
			glBindVertexArray(0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

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

