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
	*window = glfwCreateWindow(WIDTH, HEIGHT, "StokesFD", NULL, NULL);
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
	Visu->handCursor = glfwCreateStandardCursor(GLFW_HAND_CURSOR);

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


	// Declare the initial values of uniforms
	// =======================================
	if ((Grid->xmax-Grid->xmin)>(Grid->ymax-Grid->ymin)){
		Visu->scale = 2.0/(1.15*(Grid->xmax-Grid->xmin));
	}
	else {
		Visu->scale = 2.0/(1.15*(Grid->ymax-Grid->ymin));
	}

	GLint loc = glGetUniformLocation(Visu->ShaderProgram, "one_ov_log_of_10");
	glUniform1f(loc, 1.0/log(10));

	Visu->colorScale[0] = -0.5;
	Visu->colorScale[1] =  0.5;
	Visu->log10_on = 1;
	Visu->valueScale = 1.0;

	Visu->shift[0] = - ((Grid->xmax + Grid->xmin)/2.0)*Visu->scale;
	Visu->shift[1] = - ((Grid->ymax + Grid->ymin)/2.0)*Visu->scale;

	Visu->mouse1BeginDrag[0] = 0;
	Visu->mouse1BeginDrag[1] = 0;
	Visu->mouse2BeginDrag[0] = 0;
	Visu->mouse2BeginDrag[1] = 0;
	Visu->mouse1EndDrag[0] = 0;
	Visu->mouse1EndDrag[1] = 0;
	Visu->mouse2EndDrag[0] = 0;
	Visu->mouse2EndDrag[1] = 0;
	//	Visu_updateUniforms(Visu);

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

	Visu_updateVertices(Visu, Grid);

}


void Visu_updateVertices(Visu* Visu, Grid* Grid)
{
	int iy, ix, C;
	C =0;
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			Visu->vertices[C  ] = (Grid->xmin + ix*Grid->dx);
			Visu->vertices[C+1] = (Grid->ymin + iy*Grid->dy);
			C += 2;
		}
	}
}

void Visu_updateCenterValue(Visu* Visu, Grid* Grid, compute* CellValue, int BCType)
{
	// UC is a scalar CellValue defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	int iNW, iNE, iSW, iSE;
	// CellValue interpolated on the center nodes
	// ======================================
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			I = ix + iy*Grid->nxS;
			iNW = (ix-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			Visu->U[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}




	// CellValue extrapolated on the lower boundary
	// ======================================
	// o: centered CellValue
	// x: CellValue extrapolated (in 1D) fron the o nodes
	// X: value interpolated between the two x
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


		//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		//Visu->U[I] = (temp1+temp2)/2;
		Visu->U[I] = (CellValue[i1a]+CellValue[i2a])/2;
	}
	// CellValue extrapolated on the upper boundary
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
		//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		//Visu->U[I] = (temp1+temp2)/2;
		Visu->U[I] = (CellValue[i1a]+CellValue[i2a])/2;
	}


	if (BCType!=1) { // not periodic
		// CellValue extrapolated on the left boundary
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
			//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			//Visu->U[I] = (temp1+temp2)/2;
			Visu->U[I] = (CellValue[i1a]+CellValue[i2a])/2;
		}

		// CellValue extrapolated on the right boundary
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
			//temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
			//temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
			//Visu->U[I] = (temp1+temp2)/2;
			Visu->U[I] = (CellValue[i1a]+CellValue[i2a])/2;
		}

		// Lower left corner
		//          1b
		//      1a
		//   X
		ix = 0; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy+1)*Grid->nxC;
		i1a =  ix   +(iy  )*Grid->nxC;
		//Visu->U[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		Visu->U[I] = CellValue[i1a];

		// Lower right corner
		//  1b
		//      1a
		//          X
		ix = Grid->nxS-1; iy = 0;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+(iy  )*Grid->nxC;
		//Visu->U[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		Visu->U[I] = CellValue[i1a];

		// Upper left corner
		//  X
		//      1a
		//          1b
		ix = 0; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix+1)+(iy-2)*Grid->nxC;
		i1a =  ix   +(iy-1)*Grid->nxC;
		//Visu->U[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		Visu->U[I] = CellValue[i1a];

		// Upper right corner
		//          X
		//      1a
		//  1b
		ix = Grid->nxS-1; iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		//Visu->U[I] = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		Visu->U[I] = CellValue[i1a];



	}
	else { // if periodic boundaries

		for (iy = 1; iy < Grid->nyS-1; ++iy) {
			// Left and right boundary
			ix = 0;
			I = ix + iy*Grid->nxS;
			iNW = (ix+Grid->nxC-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix+Grid->nxC-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			Visu->U[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
			Visu->U[I+Grid->nxS-1] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}

		// Upper left and right corners
		// ======================================
		//   x  X  x
		//  1a    2a
		//  1b    2b
		iy = Grid->nyS-1;
		I = ix + iy*Grid->nxS;
		i1b = (Grid->nxC-1)+(iy-2)*Grid->nxC;
		i1a = (Grid->nxC-1)+(iy-1)*Grid->nxC;
		i2b =  0   +(iy-2)*Grid->nxC;
		i2a =  0   +(iy-1) *Grid->nxC;
		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
		Visu->U[I+Grid->nxS-1] = (temp1+temp2)/2;


		// Lower left and right corners
		// ======================================
		//  1b    2b
		//  1a    2a
		//   x  X  x
		iy = 0;
		compute temp1, temp2;
		int i1a, i1b, i2a, i2b;
		I = ix + iy*Grid->nxS;
		i1b = (Grid->nxC-1)+(iy+1)*Grid->nxC;
		i1a = (Grid->nxC-1)+ iy   *Grid->nxC;
		i2b =  0   +(iy+1)*Grid->nxC;
		i2a =  0   + iy   *Grid->nxC;
		temp1 = CellValue[i1a] - (CellValue[i1b] - CellValue[i1a])/2;
		temp2 = CellValue[i2a] - (CellValue[i2b] - CellValue[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
		Visu->U[I+Grid->nxS-1] = (temp1+temp2)/2;


	}


}








void Visu_strainRate(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC)
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

	Visu_updateCenterValue (Visu, Grid, CenterEps, BC->SetupType);
	free(CenterEps);
}


void Visu_velocity(Visu* Visu, Grid* Grid, Physics* Physics)
{
	int iy, ix;
	int I = 0;
    // Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
    for (iy=0; iy<Grid->nyS; iy++){
        for (ix=0; ix<Grid->nxS; ix++) {
        	I = ix+iy*Grid->nxS;
        	Visu->U[I]  = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2;
        	//Visu->U[I] = (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2;
        	//printf("%.2f  ",Visu->U[I]);
        }
        //printf("\n");
    }
}





void Visu_updateUniforms(Visu* Visu, GLFWwindow* window)
{
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	GLfloat ratio = (GLfloat)width/(GLfloat)height;
	//printf("ratio = %.2f, scale = %.2f\n\n\n\n",ratio, Visu->scale);

	GLfloat Transform[] = {Visu->scale,0.0f,0.0f,0.0f , 0.0f,Visu->scale*ratio,0.0f,0.0f , 0.0f,0.0f,1.0f,0.0f , Visu->shift[0],Visu->shift[1],0.0f,1.0f};
	GLuint loc = glGetUniformLocation(Visu->ShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);


	loc = glGetUniformLocation(Visu->ShaderProgram, "colorScale");
	glUniform2f(loc, Visu->colorScale[0], Visu->colorScale[1]);


	loc = glGetUniformLocation(Visu->ShaderProgram, "log10_on");
	glUniform1i(loc, Visu->log10_on);


	loc = glGetUniformLocation(Visu->ShaderProgram, "valueScale");
	glUniform1f(loc, Visu->valueScale);
}




void Visu_update(Visu* Visu, GLFWwindow* window, Grid* Grid, Physics* Physics, BC* BC, Char* Char)
{
	Visu_checkInput(Visu, window);

	switch (Visu->type) {
	case Viscosity:
		glfwSetWindowTitle(window, "Viscosity");
		Visu_updateCenterValue(Visu, Grid, Physics->eta, BC->SetupType);
		Visu->valueScale = 1.0;//Char->viscosity;
		Visu->colorScale[0] = -0.1;
		Visu->colorScale[1] =  0.1;
		Visu->log10_on = true;
		break;

	case StrainRate:
		glfwSetWindowTitle(window, "StrainRate");
		Visu_strainRate(Visu, Grid, Physics, BC);
		Visu->valueScale = Physics->epsRef;
		Visu->colorScale[0] = -3;
		Visu->colorScale[1] =  3;
		Visu->log10_on = true;
		break;

	case Velocity:
			glfwSetWindowTitle(window, "Velocity");
			Visu_velocity(Visu, Grid, Physics);
			Visu->valueScale = 1.0;//(Physics->epsRef*Grid->xmax);
			Visu->colorScale[0] = -1;
			Visu->colorScale[1] =  1;
			//Visu->scale 		= 2.0/(1.5*(Grid->xmax-Grid->xmin));
			Visu->log10_on = false;
			break;
	case Pressure:
			glfwSetWindowTitle(window, "Pressure");
			Visu_updateCenterValue(Visu, Grid, Physics->P, BC->SetupType);

			Visu->valueScale = 1.0;//Char->stress;
			Visu->colorScale[0] = -100;
			Visu->colorScale[1] =  100;
			Visu->log10_on = false;
			break;
	default:
		printf("Error: unknown Visu.type: %i",Visu->type);
	}

	Visu_updateUniforms(Visu, window);
}


void Visu_checkInput(Visu* Visu, GLFWwindow* window)
{
	// Check keyboard events
		if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
			Visu->type = Viscosity;
		}
		if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
			Visu->type = StrainRate;
		}
		if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
			Visu->type = Velocity;
		}
		if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) {
			Visu->type = Pressure;
		}



		// Check mouse events

		// Left click - shift visu
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);
			if (!Visu->mouse1Pressed) {
				Visu->mouse1BeginDrag[0] = xpos;
				Visu->mouse1BeginDrag[1] = ypos;
				glfwSetCursor(window,Visu->handCursor);
			}
			if (Visu->mouse1Pressed) {
				Visu->mouse1EndDrag[0] = xpos;
				Visu->mouse1EndDrag[1] = ypos;

				int width, height;
				glfwGetWindowSize(window, &width, &height);
				Visu->shift[0] += (Visu->mouse1EndDrag[0] - Visu->mouse1BeginDrag[0])/width*2.0;
				Visu->shift[1] -= (Visu->mouse1EndDrag[1] - Visu->mouse1BeginDrag[1])/height*2.0;

				Visu->mouse1BeginDrag[0] = xpos;
				Visu->mouse1BeginDrag[1] = ypos;

				//printf("xpos=%.1f, ypos%.1f\n\n\n", xpos,ypos);
			}
			Visu->mouse1Pressed = true;
		}
		else {

			Visu->mouse1Pressed = false;
			glfwSetCursor(window,NULL);
		}


		// Left click - zoom
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS){
			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);
			if (!Visu->mouse2Pressed) {
				Visu->mouse2BeginDrag[0] = xpos;
				Visu->mouse2BeginDrag[1] = ypos;
				glfwSetCursor(window,Visu->handCursor);
			}
			if (Visu->mouse2Pressed) {
				Visu->mouse2EndDrag[0] = xpos;
				Visu->mouse2EndDrag[1] = ypos;

				int width, height;
				glfwGetWindowSize(window, &width, &height);
				//Visu->shift[0] += (Visu->mouse2EndDrag[0] - Visu->mouse2BeginDrag[0])/width*2.0;
				Visu->scale *= 1+(Visu->mouse2EndDrag[1] - Visu->mouse2BeginDrag[1])/height;

				//Visu->shift[0] = (Visu->mouse2BeginDrag[0])/width*2.0 - 1.0;
				//Visu->shift[1] = (Visu->mouse2BeginDrag[1])/height*2.0;

				//Visu->mouse2BeginDrag[0] = xpos;
				//Visu->mouse2BeginDrag[1] = ypos;

				//printf("xpos=%.1f, ypos%.1f\n\n\n", xpos,ypos);
			}
			Visu->mouse2Pressed = true;
		}
		else {

			Visu->mouse2Pressed = false;
			glfwSetCursor(window,NULL);
		}




}







