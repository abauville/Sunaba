/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Visu_allocateMemory( Visu* Visu, Grid* Grid )
{
	Visu->elements      = (GLuint*)   malloc(Visu->ntrivert    		* sizeof( GLuint ));
	Visu->U             = (GLfloat*)  malloc(Grid->nxS*Grid->nyS    * sizeof( GLfloat ));
	Visu->vertices      = (GLfloat*)  malloc(Grid->nxS*Grid->nyS*2  * sizeof( GLfloat ));
}




void Visu_freeMemory( Visu* Visu )
{
	free(Visu->elements);
	free(Visu->U);
	glDeleteProgram(Visu->ShaderProgram);
	glDeleteVertexArrays(1, &Visu->VAO );
	glDeleteBuffers(1, &Visu->VBO);
	glDeleteBuffers(1, &Visu->CBO);
	glDeleteBuffers(1, &Visu->EBO);

}


void Visu_initWindow(GLFWwindow** window){

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
	*window = glfwCreateWindow(1024, 1024, "Simple example", NULL, NULL);
	if (!*window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(*window);
	glfwSetKeyCallback(*window, key_callback);

	/// Init Glew - Must be done after glut is initialized!
	// =======================================
	glewExperimental = GL_TRUE;
	GLenum res = glewInit();
	if (res != GLEW_OK) {
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		//return 1;
	}
	if(!GLEW_VERSION_3_2){
		fprintf(stderr, "OpenGL 3.2 API is not available.");
		//return 1;
	}

	/// Test GL version
	// =======================================
	const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
	const GLubyte* version = glGetString (GL_VERSION); // version as a string
	const GLubyte* glslversion = glGetString (GL_SHADING_LANGUAGE_VERSION); // version as a string

	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);
	printf("GLSL version supported %s\n", glslversion);
}







void Visu_initOpenGL(Visu* Visu, Grid* Grid) {


	// And assigned them to objects (stored in the graphic memory)
	// =======================================
	glGenVertexArrays(1, &Visu->VAO);
	glGenBuffers(1, &Visu->VBO);
	glGenBuffers(1, &Visu->CBO);
	glGenBuffers(1, &Visu->EBO);

	// Bind Vertex Array object
	// =======================================
	glBindVertexArray(Visu->VAO);
	// compile shaders
	// =======================================
	compileShaders(&Visu->ShaderProgram, Visu->VertexShaderFile, Visu->FragmentShaderFile);

	glUseProgram(Visu->ShaderProgram);

	// Get IDs for the in attributes of the shader
	// =======================================
	GLint VertAttrib    = glGetAttribLocation(Visu->ShaderProgram,"in_Vertex");
	GLint SolAttrib     = glGetAttribLocation(Visu->ShaderProgram,"U");
	// Bind objects and associate with data tables
	// =======================================
	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
	glBufferData(GL_ARRAY_BUFFER, Grid->nxS*Grid->nyS*sizeof(coord), Visu->vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, Visu->CBO);
	glBufferData(GL_ARRAY_BUFFER, Grid->nxS*Grid->nyS*sizeof(GLfloat), Visu->U, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, Visu->ntrivert*sizeof( GLuint ), Visu->elements, GL_STATIC_DRAW);

	// Connect Vertex data (stored in Visu->VBO) to the "in_Vertex" attribute of the shader
	// =======================================
	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
	glVertexAttribPointer(VertAttrib, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(VertAttrib);

	// Connect Color data (stored in Visu->CBO) to the "in_Color" attribute of the shader
	// =======================================
	glBindBuffer(GL_ARRAY_BUFFER, Visu->CBO);
	glVertexAttribPointer(SolAttrib, 1, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(SolAttrib);


	// Declare the Visu->scale as a uniform
	// =======================================


	if ((Grid->xmax-Grid->xmin)>(Grid->ymax-Grid->ymin)){
		Visu->scale = 2.0/(1.1*(Grid->xmax-Grid->xmin));
	}
	else {
		Visu->scale = 2.0/(1.1*(Grid->ymax-Grid->ymin));
	}

	GLfloat Transform[] = {Visu->scale,0.0f,0.0f,0.0f , 0.0f,Visu->scale,0.0f,0.0f , 0.0f,0.0f,1.0f,0.0f , 0.0f,0.0f,0.0f,1.0f};
	GLuint transformLoc = glGetUniformLocation(Visu->ShaderProgram, "transform");
	glUniformMatrix4fv(transformLoc, 1, GL_FALSE, &Transform[0]);


	// unbind the Buffer object (Visu->VBO, Visu->CBO) and Vertex array object (Visu->VAO)
	// =======================================
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	glUseProgram(0);




}


void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}


void Visu_init(Visu* Visu, Grid* Grid)
{

	// Create the element array
	// Fill elements, loop through cells
	int ix,iy,C;
	int nxS = Grid->nxS;
	int nyS = Grid->nyS;


	C = 0;
	for (iy=0;iy<nyS-1;iy++){
		for (ix=0;ix<nxS-1;ix++){
			// Triangle 1
			Visu->elements[C+0] = ix+iy*nxS;
			Visu->elements[C+1] = ix+1+iy*nxS;
			Visu->elements[C+2] = (iy+1)*nxS+ix;
			// Triangle 2
			Visu->elements[C+3] = ix+1+iy*nxS;
			Visu->elements[C+4] = (iy+1)*nxS+ix+1;
			Visu->elements[C+5] = (iy+1)*nxS+ix;
			C = C+6;
		}
	}

	C =0;
	for (iy = 0; iy < nyS; ++iy) {
		for (ix = 0; ix < nxS; ++ix) {
			Visu->vertices[C  ] = (Grid->xmin + ix*Grid->dx);
			Visu->vertices[C+1] = (Grid->ymin + iy*Grid->dy);
			C += 2;
		}
	}

}


void Visu_updateCenterValue(Visu* Visu, Grid* Grid, compute* Value)
{
	// UC is a scalar value defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	int iNW, iNE, iSW, iSE;
	// Value interpolated on the center nodes
	// ======================================
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			I = ix + iy*Grid->nxS;
			iNW = (ix-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			Visu->U[I] = (Value[iNW] + Value[iNE] + Value[iSW] + Value[iSE])/4;
		}
	}
	// Value extrapolated on the lower boundary
	// ======================================
	// o: centered value
	// x: value extrapolated (in 1D) fron the o nodes
	// X: valu interpolated between the two x
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1b   and 2b
	// |       |       |
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1a   and 1b
	// |       |       |
	// | - - - | - - - | -
	// |       |       |
	// | - x - X - x - |       nodes tempa and tempb
	//
	iy = 0;
	compute temp1, temp2;
	int i1a, i1b, i2a, i2b;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+ iy   *Grid->nxC;
		i2b =  ix   +(iy+1)*Grid->nxC;
		i2a =  ix   + iy   *Grid->nxC;


		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}
	// Value extrapolated on the upper boundary
	// ======================================
	//   x  X  x
	//  1a    2a
	//  1b    2b
	iy = Grid->nyS-1;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		i2b =  ix   +(iy-2)*Grid->nxC;
		i2a =  ix   +(iy-1) *Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Value extrapolated on the left boundary
	// ======================================
	//  x 1a   1b
	//  X
	//  x 2a   2b
	ix = 0;
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy  )*Grid->nxC;
		i1a =  ix   +(iy  )*Grid->nxC;
		i2b = (ix+1)+(iy-1)*Grid->nxC;
		i2a =  ix   +(iy-1)*Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Value extrapolated on the right boundary
	// ======================================
	//  1b   1a x
	//          X
	//  2b   2a x
	ix = Grid->nxS-1;
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy  )*Grid->nxC;
		i1a = (ix-1)+(iy  )*Grid->nxC;
		i2b = (ix-2)+(iy-1)*Grid->nxC;
		i2a = (ix-1)+(iy-1)*Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Lower left corner
	//          1b
	//      1a
	//   X
	ix = 0; iy = 0;
	I = ix + iy*Grid->nxS;
	i1b = (ix+1)+(iy+1)*Grid->nxC;
	i1a =  ix   +(iy  )*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Lower right corner
	//  1b
	//      1a
	//          X
	ix = Grid->nxS-1; iy = 0;
	I = ix + iy*Grid->nxS;
	i1b = (ix-2)+(iy+1)*Grid->nxC;
	i1a = (ix-1)+(iy  )*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Upper left corner
	//  X
	//      1a
	//          1b
	ix = 0; iy = Grid->nyS-1;
	I = ix + iy*Grid->nxS;
	i1b = (ix+1)+(iy-2)*Grid->nxC;
	i1a =  ix   +(iy-1)*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Upper right corner
	//          X
	//      1a
	//  1b
	ix = Grid->nxS-1; iy = Grid->nyS-1;
	I = ix + iy*Grid->nxS;
	i1b = (ix-2)+(iy-2)*Grid->nxC;
	i1a = (ix-1)+(iy-1)*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;


	if (DEBUG) {
		int C = 0;
		printf("=== Check Value ===\n");
		for (iy = 0; iy < Grid->nyC; ++iy) {
			for (ix = 0; ix < Grid->nxC; ++ix) {
				printf("%.2f  ", Value[C]);
				C++;
			}
			printf("\n");
		}


		C = 0;
		printf("=== Check U ===\n");
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				printf("%.2f  ", Visu->U[C]);
				C++;
			}
			printf("\n");
		}

	}


}








void Visu_StrainRate(Visu* Visu, Grid* Grid, Physics* Physics)
{

	compute* CenterEps = (compute*) malloc(Grid->nCTot * sizeof(compute));

	// Definition of second invariant: // E_II = sqrt( Eps_xx^2 + Eps_xy^2  );
	// Declarations
	// =========================
	int ix, iy, I, iNode, Ix, Iy;
	compute dVxdy, dVydx, dVxdx;
	// ix, iy modifiers
	int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
	int IyMod[4] = {0,0,1,1};
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			I = ix+iy*Grid->nxC;

			// Compute Eps_xy at the four nodes of the cell
			// 1. Sum contributions
			dVxdy = 0;
			dVydx = 0;
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = ix+IxMod[iNode];
				Iy = iy+IyMod[iNode];
				dVxdy += ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
						 - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;

				dVydx += ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
						 - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;
			}
			// 2. Average
			dVxdy /= 4;
			dVydx /= 4;

			dVxdx = (Physics->Vx[(ix+1) + (iy+1)*Grid->nxVx]
				   - Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx])/Grid->dx;

			CenterEps[I] = sqrt(  (0.5*(dVxdy+dVydx))*(0.5*(dVxdy+dVydx))    +    dVxdx*dVxdx  );
			//CenterEps[I] = sqrt(  dVxdx*dVxdx  );
		}
	}

	/*
	printf("=== Check Vx ===\n");
	int C = 0;
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			printf("%.2f  ", Physics->Vx[C]);
			C++;
		}
		printf("\n");
	}

	printf("=== Check CenterEps ===\n");
	C = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			printf("%.2f  ", CenterEps[C]);
			C++;
		}
		printf("\n");
	}
	*/

	Visu_updateCenterValue (Visu, Grid, CenterEps);
	free(CenterEps);
}
























