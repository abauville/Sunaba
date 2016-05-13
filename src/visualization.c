/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Visu_allocateMemory( Visu* Visu, Grid* Grid )
{
	Visu->U             = (GLfloat*)  malloc(2*Grid->nxS*Grid->nyS    * sizeof( GLfloat ));
	//Visu->vertices      = (GLfloat*)  malloc(Grid->nxS*Grid->nyS*2  * sizeof( GLfloat ));
	Visu->elements      = (GLuint*)   malloc(Visu->ntrivert    		* sizeof( GLuint ));

	Visu->vertices      = (GLfloat*)  malloc(4 * 4 * sizeof( GLfloat )); // 4 corners only
	Visu->particles 	= (GLfloat*) malloc (Visu->nParticles*4*sizeof(GLfloat));
	printf("%i  \n", (Visu->particleMeshRes+1) *3);
	Visu->particleMesh 	= (GLfloat*) malloc ((Visu->particleMeshRes+2) *3*sizeof(GLfloat));
	//Visu->elements      = (GLuint*)   malloc(6  * sizeof( GLuint  )); // 2 triangles


	Visu->imageBuffer 	= (unsigned char*) malloc(Visu->retinaScale*Visu->retinaScale*3*WIDTH*HEIGHT*sizeof(unsigned char)); // does not consider image resizing

}




void Visu_freeMemory( Visu* Visu )
{
	free(Visu->elements);
	free(Visu->U);
	free(Visu->particles);
	free(Visu->particleMesh);
	free(Visu->imageBuffer);
	glDeleteProgram(Visu->ShaderProgram);
	glDeleteVertexArrays(1, &Visu->VAO );
	glDeleteVertexArrays(1, &Visu->VAO_part);
	glDeleteBuffers(1, &Visu->VBO);
	glDeleteBuffers(1, &Visu->VBO_part);
	glDeleteBuffers(1, &Visu->VBO_partMesh);
	glDeleteBuffers(1, &Visu->EBO);
	glDeleteTextures(1, &Visu->TEX);

}


void Visu_initWindow(GLFWwindow** window, Visu* Visu){

	glfwSetErrorCallback(error_callback);
	if (!glfwInit()){
		exit(EXIT_FAILURE);
	}
//#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
//#endif
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

	Visu->handCursor = glfwCreateStandardCursor(GLFW_HAND_CURSOR);
	Visu->paused 	 = false;

	/// Init Glew - Must be done after glut is initialized!
	// =======================================
//#ifdef __APPLE__
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
//#endif
	/// Test GL version
	// =======================================
	const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
	const GLubyte* version = glGetString (GL_VERSION); // version as a string
	const GLubyte* glslversion = glGetString (GL_SHADING_LANGUAGE_VERSION); // version as a string

	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);
	printf("GLSL version supported %s\n", glslversion);






}


void Visu_particles(Visu* Visu, Particles* Particles, Grid* Grid)
{
	if (Visu->nParticles<Particles->n)
	{
		Visu->nParticles = Particles->n + (int)(Particles->n*0.1);
		Visu->particles = (GLfloat*) realloc (Visu->particles, Visu->nParticles*4*sizeof(GLfloat));
		// Here I assume that the Visu->VBOPart is bound
		glBufferData(GL_ARRAY_BUFFER, 4*Visu->nParticles*sizeof(GLfloat), NULL, GL_STREAM_DRAW);
	}

	int C = 0;
	INIT_PARTICLE
	FOR_PARTICLES
		Visu->particles[C] = thisParticle->x;
		Visu->particles[C+1] = thisParticle->y;

		if (Visu->typeParticles == Phase) {
			Visu->particles[C+2] = thisParticle->phase;
		} else if (Visu->typeParticles == PartTemp) {
			Visu->particles[C+2] = thisParticle->T;
		} else if (Visu->typeParticles == PartSigma_xx) {
			Visu->particles[C+2] = thisParticle->sigma_xx_0;
		} else if (Visu->typeParticles == PartSigma_xy) {
			Visu->particles[C+2] = thisParticle->sigma_xy_0;
		}
		Visu->particles[C+3] = thisParticle->passive;

		C += 4;
	END_PARTICLES

}

void Visu_particleMesh(Visu* Visu)
{


	// Create the particle mesh, i.e. cone
	GLfloat radius = 2.0;// 1.0=dx or dy


	int i;
	Visu->particleMesh[ 0] = 0.0;
	Visu->particleMesh[ 1] = 0.0;
	Visu->particleMesh[ 2] = -0.1;
	/*
	Visu->particleMesh[ 3] = -radius;
	Visu->particleMesh[ 4] = -radius;
	Visu->particleMesh[ 5] =  0.1;

	Visu->particleMesh[ 6] =  radius;
	Visu->particleMesh[ 7] = -radius;
	Visu->particleMesh[ 8] =  0.1;

	Visu->particleMesh[ 9] = radius;
	Visu->particleMesh[10] = radius;
	Visu->particleMesh[11] =  0.1;

	Visu->particleMesh[12] = -radius;
	Visu->particleMesh[13] =  radius;
	Visu->particleMesh[14] =  0.1;

	Visu->particleMesh[15] = -radius;
	Visu->particleMesh[16] = -radius;
	Visu->particleMesh[17] =  0.1;
*/

	int C = 3;

	for (i=0;i<(Visu->particleMeshRes+1);i++) {
		Visu->particleMesh[C] 	= radius * cos((i*2*PI+PI)/Visu->particleMeshRes);
		Visu->particleMesh[C+1] = radius * sin((i*2*PI+PI)/Visu->particleMeshRes);
		Visu->particleMesh[C+2] = 0.1;
		//printf("%.3f   %.3f   %.3f\n",Visu->particleMesh[C],Visu->particleMesh[C+1],Visu->particleMesh[C+2]);
		C+=3;
	}

	//printf("C = %i", C);



}



void Visu_initOpenGL(Visu* Visu, Grid* Grid, GLFWwindow* window) {

	///Init shader
	// =======================================
	Visu->VertexShaderFile 				= "src/shader.vs";
	Visu->FragmentShaderFile 			= "src/shader.fs";
	Visu->ParticleVertexShaderFile 		= "src/particleShader.vs";
	Visu->ParticleGeometryShaderFile 	= "src/particleShader.gs";
	Visu->ParticleFragmentShaderFile 	= "src/particleShader.fs";
	Visu->ParticleBackgroundVertexShaderFile 	= "src/particleBackgroundShader.vs";
	Visu->ParticleBackgroundFragmentShaderFile 	= "src/particleBackgroundShader.fs";


	Visu->ShaderProgram = 0;
	Visu->ParticleShaderProgram = 0;
	Visu->ParticleBackgroundShaderProgram = 0;

	// Generate reference to objects (indexes that act as pointers to graphic memory)
	// =======================================
	Visu->VAO = 0; // Reference to the Vertex   array object
	Visu->VBO = 0; // Reference to the Vertex   buffer object
	Visu->EBO = 0; // Reference to the Element  buffer object
	Visu->TEX = 0; // Reference to the Element  buffer object



	// And assigned them to objects (stored in the graphic memory)
	// =======================================
	glGenVertexArrays(1, &Visu->VAO);
	glGenVertexArrays(1, &Visu->VAO_part);
	glGenBuffers(1, &Visu->VBO);
	glGenBuffers(1, &Visu->EBO);
	glGenTextures(1, &Visu->TEX);
	glGenBuffers(1, &Visu->VBO_part);
	glGenBuffers(1, &Visu->VBO_partMesh);

	// Bind Vertex Array object
	// =======================================
	glBindVertexArray(Visu->VAO);
	// compile shaders
	// =======================================
	const char* dum = NULL;
	compileShaders(&Visu->ShaderProgram, Visu->VertexShaderFile, Visu->FragmentShaderFile, dum, false);
	printf("Grid Shader succesfully compiled\n");
	glUseProgram(Visu->ShaderProgram);


	// Get IDs for the in attributes of the shader
	// =======================================
	GLint VertAttrib    	 = glGetAttribLocation(Visu->ShaderProgram,"in_Vertex");
	GLint TexCoordAttrib     = glGetAttribLocation(Visu->ShaderProgram,"in_TexCoord");

	// Bind objects and associate with data tables
	// =======================================
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, Visu->ntrivert*sizeof( GLuint ), Visu->elements, GL_STATIC_DRAW);



	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
		glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu->vertices, GL_STATIC_DRAW);
		glVertexAttribPointer(VertAttrib, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), 0);
		glEnableVertexAttribArray(VertAttrib);

		glVertexAttribPointer(TexCoordAttrib, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), (void*)(2*sizeof(GLfloat)));
		glEnableVertexAttribArray(TexCoordAttrib);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glBindTexture(GL_TEXTURE_2D, Visu->TEX);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RG, Grid->nxS, Grid->nyS, 0, GL_RG, GL_FLOAT, Visu->U);
	glBindTexture(GL_TEXTURE_2D, 0);










	// Declare the initial values of uniforms
	// =======================================
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	GLfloat ratio = (GLfloat)width/(GLfloat)height;

	if ((Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0])>(Grid->ymax-Grid->ymin)*(1+2*Visu->shiftFac[1])/ratio){
		Visu->scale = 2.0/(1.05*(Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0]));
	}
	else {
		Visu->scale = 2.0/(1.05*(Grid->ymax-Grid->ymin)*(1+2*Visu->shiftFac[1])*ratio);
	}

	GLint loc = glGetUniformLocation(Visu->ShaderProgram, "one_ov_log_of_10");
	glUniform1f(loc, 1.0/log(10));

	Visu->colorScale[0] = -0.5;
	Visu->colorScale[1] =  0.5;
	Visu->log10_on 		= 1;
	Visu->valueScale 	= 1.0;
	Visu->valueShift 	= 0.0;

	Visu->shift[0] = - ((Grid->xmax + Grid->xmin)/2.0)*Visu->scale;
	Visu->shift[1] = - ((Grid->ymax + Grid->ymin)/2.0)*Visu->scale;
	Visu->shift[2] = - 0.0;

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







	// =======================================
	// Particles part
	// =======================================
	glBindVertexArray(Visu->VAO_part);

	compileShaders(&Visu->ParticleShaderProgram, Visu->ParticleVertexShaderFile, Visu->ParticleFragmentShaderFile, Visu->ParticleGeometryShaderFile, false);
	glUseProgram(Visu->ParticleShaderProgram);

	GLint ParticleVertAttrib    	= glGetAttribLocation(Visu->ParticleShaderProgram,"PartVertex");
	GLint ParticleData    	 		= glGetAttribLocation(Visu->ParticleShaderProgram,"PartData");
	GLint ParticlePassiveData    	= glGetAttribLocation(Visu->ParticleShaderProgram,"PartPassiveData");
	GLint ParticleMeshVertex    	= glGetAttribLocation(Visu->ParticleShaderProgram,"PartMeshVertex");



	glBindVertexArray(0);




	// Create Mesh

	glBindVertexArray(Visu->VAO_part);
	Visu_particleMesh(Visu);


	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_partMesh);

		glBufferData(GL_ARRAY_BUFFER, 3*(Visu->particleMeshRes+2)*sizeof(GLfloat), Visu->particleMesh, GL_STATIC_DRAW);
		glVertexAttribPointer(ParticleMeshVertex , 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), 0);
		glEnableVertexAttribArray(ParticleMeshVertex );
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_part);
		glBufferData(GL_ARRAY_BUFFER, 4*Visu->nParticles*sizeof(GLfloat), NULL, GL_STREAM_DRAW);
		glVertexAttribPointer(ParticleVertAttrib , 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), 0);
		glEnableVertexAttribArray(ParticleVertAttrib );
		glVertexAttribPointer(ParticleData , 1, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), (void*)(2*sizeof(GLfloat)));
		glEnableVertexAttribArray(ParticleData );
		glVertexAttribPointer(ParticlePassiveData , 1, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));
		glEnableVertexAttribArray(ParticlePassiveData );
	glBindBuffer(GL_ARRAY_BUFFER, 0);




		glVertexAttribDivisor(ParticleMeshVertex , 0); // never changes
		glVertexAttribDivisor(ParticleVertAttrib, 1); // counter of +1 per instance
		glVertexAttribDivisor(ParticleData, 1);
		glVertexAttribDivisor(ParticlePassiveData, 1);

		glBindVertexArray(0);
	glUseProgram(0);



	// =======================================
	// Particles backgroundShader
	// =======================================
	// the VAO, VBO is the same as for the textured rectangle (that is used to plot the data on the grid)
	// only the shaders are different
	glBindVertexArray(Visu->VAO);
	compileShaders(&Visu->ParticleBackgroundShaderProgram, Visu->ParticleBackgroundVertexShaderFile, Visu->ParticleBackgroundFragmentShaderFile, dum, false);
	printf("particle background shader succesfully compiled\n");
	glUseProgram(Visu->ParticleBackgroundShaderProgram);
	// Get IDs for the in attributes of the shader
	// =======================================
	GLint PartBackVertAttrib    	 = glGetAttribLocation(Visu->ParticleBackgroundShaderProgram,"in_Vertex");

	// Bind objects and associate with data tables
	// =======================================
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, Visu->ntrivert*sizeof( GLuint ), Visu->elements, GL_STATIC_DRAW);



	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
		glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu->vertices, GL_STATIC_DRAW);
		glVertexAttribPointer(PartBackVertAttrib, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), 0);
		glEnableVertexAttribArray(PartBackVertAttrib);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);




	glUseProgram(0);
	glBindVertexArray(0);















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


void Visu_init(Visu* Visu, Grid* Grid, Particles* Particles)
{

	Visu->elements[0] = 0;
	Visu->elements[1] = 1;
	Visu->elements[2] = 2;
	Visu->elements[3] = 2;
	Visu->elements[4] = 3;
	Visu->elements[5] = 1;

	Visu_updateVertices(Visu, Grid);
	Visu_particles(Visu, Particles, Grid);
}


void Visu_updateVertices(Visu* Visu, Grid* Grid)
{
	/*
		int iy, ix, C;
		C =0;
		for (iy = 0; iy < Grid->nyS; ++iy) {
			for (ix = 0; ix < Grid->nxS; ++ix) {
				Visu->vertices[C  ] = (Grid->xmin + ix*Grid->dx);
				Visu->vertices[C+1] = (Grid->ymin + iy*Grid->dy);
				C += 2;
			}
		}
	*/
	// Coordinates of a simple rectangle for the texture;
	int ix, iy;
	int C = 0;
	for (iy = 0; iy < 2; ++iy) {
		for (ix = 0; ix < 2; ++ix) {
			Visu->vertices[C  ] = Grid->xmin + ix*(Grid->xmax-Grid->xmin);
			Visu->vertices[C+1] = Grid->ymin + iy*(Grid->ymax-Grid->ymin);
			Visu->vertices[C+2] = 1.0*ix;
			Visu->vertices[C+3] = 1.0*iy;

			C += 4;
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
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			I = 2* (ix + iy*Grid->nxS);
			iNW = (ix)+ (iy+1)   *Grid->nxEC;
			iNE = ix+1    + (iy+1)   *Grid->nxEC;
			iSW = (ix)+(iy)*Grid->nxEC;
			iSE = ix+1    +(iy)*Grid->nxEC;
			Visu->U[I] = (CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}



}








void Visu_strainRate(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC)
{

	compute* CenterEps = (compute*) malloc(Grid->nECTot * sizeof(compute));


	Physics_computeStrainRateInvariant(Physics, Grid, CenterEps);

	/*
	int iy, ix;
	int C = 0;
	printf("=== Cehck StrainRate invariant");
	for (iy=0;iy<Grid->nyEC;iy++) {
		for (ix=0;ix<Grid->nxEC;ix++) {
			printf("%.2e ", CenterEps[C]);
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
	compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	for (iy=0; iy<Grid->nyS; iy++){
		for (ix=0; ix<Grid->nxS; ix++) {
			I = 2*(ix+iy*Grid->nxS);
			//Visu->U[I]  = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2;
			//Visu->U[I] += (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2;


			A  = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2;
			B  = (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2;
			Visu->U[I] = sqrt(A*A + B*B);

		}
		//printf("\n");
	}
}


void Visu_stress(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC)
{

	int iy, ix;
	int I = 0;
	//compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	for (iy=0; iy<Grid->nyS; iy++){
		for (ix=0; ix<Grid->nxS; ix++) {
			I = 1*(ix+iy*Grid->nxS);


			Visu->U[2*I] = Physics->Dsigma_xy_0[I];

		}
		//printf("\n");
	}



	//compute* CenterEps = (compute*) malloc(Grid->nECTot * sizeof(compute));


	//Physics_computeStrainRateInvariant(Physics, Grid, CenterEps);





	/*
	int iy, ix;
	int C = 0;
	printf("=== Cehck StrainRate invariant");
	for (iy=0;iy<Grid->nyEC;iy++) {
		for (ix=0;ix<Grid->nxEC;ix++) {
			printf("%.2e ", CenterEps[C]);
			C++;
		}
		printf("\n");
	}
	*/



	//Visu_updateCenterValue (Visu, Grid, CenterEps, BC->SetupType);
	//free(CenterEps);
}


void Visu_alphaValue(Visu* Visu, Grid* Grid, Particles* Particles) {
	// Based on phase
	float C = 0;
	float alpha;
	INIT_PARTICLE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		thisParticle = Particles->linkHead[iNode];
		C = 0;
		alpha = 0;
		while (thisParticle != NULL) {
			if (thisParticle->phase!=0) {
				alpha ++;
			}
			C++;
			thisParticle = thisParticle->next;
		}

		Visu->U[2*iNode+1] = alpha/C;
		if (Visu->U[2*iNode+1]<0.99999999) {
			Visu->U[2*iNode+1] = 0;
		}
	}

}








void Visu_updateUniforms(Visu* Visu, GLFWwindow* window)
{
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	GLfloat ratio = (GLfloat)width/(GLfloat)height;
	//printf("ratio = %.2f, scale = %.2f\n\n\n\n",ratio, Visu->scale);

	GLfloat Transform[] = {Visu->scale,0.0f,0.0f,0.0f , 0.0f,Visu->scale*ratio,0.0f,0.0f , 0.0f,0.0f,1.0f,0.0f , Visu->shift[0],Visu->shift[1]*ratio,Visu->shift[2],1.0f};
	GLuint loc = glGetUniformLocation(Visu->ShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->ParticleBackgroundShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "size");
	//printf("scale: %.3f\n",Visu->scale);
	glUniform1f(loc, 1.0*Visu->scale);


	//printf("scale: %.3f\n",Visu->scale);
	GLfloat type;
	if (Visu->typeParticles == Phase ) {
		type = 0;
	} else {
		type = 1;
	}
	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "type");
	glUniform1f(loc, type);

	loc = glGetUniformLocation(Visu->ShaderProgram, "colorScale");
	glUniform2f(loc, Visu->colorScale[0], Visu->colorScale[1]);

	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "colorScale");
	glUniform2f(loc, Visu->partColorScale[0], Visu->partColorScale[1]);


	loc = glGetUniformLocation(Visu->ShaderProgram, "log10_on");
	glUniform1i(loc, Visu->log10_on);


	loc = glGetUniformLocation(Visu->ShaderProgram, "valueScale");
	glUniform1f(loc, Visu->valueScale);

	loc = glGetUniformLocation(Visu->ShaderProgram, "valueShift");
	glUniform1f(loc, Visu->valueShift);

	loc = glGetUniformLocation(Visu->ShaderProgram, "transparency");
	glUniform1i(loc, Visu->transparency);

}




void Visu_update(Visu* Visu, GLFWwindow* window, Grid* Grid, Physics* Physics, BC* BC, Char* Char)
{

	int i;
	switch (Visu->type) {
	case Viscosity:
		glfwSetWindowTitle(window, "Viscosity");
		Visu->valueScale = 1.0;//Char->viscosity;
		Visu->valueShift = 0;
		Visu_updateCenterValue(Visu, Grid, Physics->eta, BC->SetupType);


		Visu->colorScale[0] = -2;
		Visu->colorScale[1] =  2;
		Visu->log10_on = true;
		break;

	case StrainRate:
		glfwSetWindowTitle(window, "StrainRate");
		Visu->valueScale = Physics->epsRef;
		Visu->valueShift = 0;
		Visu_strainRate(Visu, Grid, Physics, BC);

		Visu->colorScale[0] = -0.5;
		Visu->colorScale[1] =  0.5;
		Visu->log10_on = true;
		break;
	case Stress:
		glfwSetWindowTitle(window, "Stress");
		Visu->valueScale = 1.0;
		Visu->valueShift = 0;
		Visu_stress(Visu, Grid, Physics, BC);

		Visu->colorScale[0] = -1;
		Visu->colorScale[1] =  1;
		Visu->log10_on = false;
		break;
	case Velocity:
		glfwSetWindowTitle(window, "Velocity");
		Visu_velocity(Visu, Grid, Physics);
		Visu->valueScale = 0.5*Physics->maxV;//(Physics->epsRef*Grid->xmax);
		Visu->valueShift = 0;
		Visu->colorScale[0] = -1.5;
		Visu->colorScale[1] =  1.5;
		Visu->log10_on = true;
		break;

	case Pressure:
		glfwSetWindowTitle(window, "Pressure");
		Visu_updateCenterValue(Visu, Grid, Physics->P, BC->SetupType);

		Visu->valueScale = 1.0;//Char->stress;
		Visu->valueShift = 0;
		Visu->colorScale[0] = -200;
		Visu->colorScale[1] =  200;
		Visu->log10_on = false;
		break;
	case Density:
		glfwSetWindowTitle(window, "Density");
		Visu_updateCenterValue(Visu, Grid, Physics->rho, BC->SetupType);
		Visu->valueScale = 1.0;
		Visu->valueShift = 0;
		Visu->colorScale[0] = -0.1;
		Visu->colorScale[1] =  0.1;
		Visu->log10_on = true;
		break;
	case Temperature:
			glfwSetWindowTitle(window, "Temperature");
			Visu_updateCenterValue(Visu, Grid, Physics->T, BC->SetupType); // Not optimal but good enough for the moment
			Visu->valueScale = 1.0;

			Visu->colorScale[0] = -0.5;
			Visu->colorScale[1] =  0.5;
			Visu->valueShift = 1*Visu->colorScale[0];
			Visu->log10_on = false;

			break;
	case WaterPressureHead:
			glfwSetWindowTitle(window, "Water pressure head");
			Visu_updateCenterValue(Visu, Grid, Physics->psi, BC->SetupType); // Not optimal but good enough for the moment
			Visu->valueScale = 1.0;

			Visu->colorScale[0] = -0.5;
			Visu->colorScale[1] =  0.5;
			Visu->valueShift = 1*Visu->colorScale[0];
			Visu->log10_on = false;

			break;

	case Blank:
			glfwSetWindowTitle(window, "Blank");
			for (i=0;i<Grid->nSTot;i++) {
				Visu->U[i] = 0;
			}
			Visu->valueScale = 1.0;
			Visu->colorScale[0] = -1;
			Visu->colorScale[1] =  1;
			Visu->valueShift = 0;
			Visu->log10_on = false;
			break;
	default:
		printf("Error: unknown Visu.type: %i",Visu->type);
	}

	switch (Visu->typeParticles) {
	case Phase:
		Visu->partColorScale[0] = -3;
		Visu->partColorScale[1] =  3;
		break;
	case PartTemp:
		Visu->partColorScale[0] = -1;
		Visu->partColorScale[1] =  1;
		break;
	case PartSigma_xx:
		Visu->partColorScale[0] = -100000;
		Visu->partColorScale[1] =  100000;
		break;
	case PartSigma_xy:
		Visu->partColorScale[0] = -100000;
		Visu->partColorScale[1] =  100000;
		break;
	default:
		printf("Error: unknown Visu.typeParticles: %i",Visu->typeParticles);
	}




	Visu_updateUniforms(Visu, window);
}


void Visu_checkInput(Visu* Visu, GLFWwindow* window)
{
	// Check keyboard events
	if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
		Visu->type = Viscosity;
	}
	else if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
		Visu->type = StrainRate;
	}
	else if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
		Visu->type = Stress;
	}
	else if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) {
		Visu->type = Pressure;
	}
	else if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS) {
		Visu->type = Density;
	}
	else if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS) {
		Visu->type = Temperature;
	}
	else if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS) {
		Visu->type = Velocity;
	}
	else if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS) {
		Visu->type = WaterPressureHead;
	}
	else if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS) {
		Visu->type = Blank;
	}



	else if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
		Visu->typeParticles = Phase;
	}
	else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		Visu->typeParticles = PartTemp;
	}
	else if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
		Visu->typeParticles = PartSigma_xx;
	}
	else if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
		Visu->typeParticles = PartSigma_xy;
	}





	else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		Visu->paused = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
		Visu->paused = false;
	}
	else if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
		Visu->showParticles = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) {
		Visu->showParticles = false;
	}
	else if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		Visu->initPassivePart = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
		Visu->transparency = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS) {
		Visu->transparency = false;
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
			Visu->shift[1] -= (Visu->mouse1EndDrag[1] - Visu->mouse1BeginDrag[1])/height*1.0;

			Visu->mouse1BeginDrag[0] = xpos;
			Visu->mouse1BeginDrag[1] = ypos;

		}
		Visu->mouse1Pressed = true;
	}

	// Righr click - zoom
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS){
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		int width, height;
		glfwGetWindowSize(window, &width, &height);

		if (!Visu->mouse2Pressed) {
			Visu->mouse2BeginDrag[0] = xpos/width*2-1;
			Visu->mouse2BeginDrag[1] = ypos/height*1-0.5;
			glfwSetCursor(window,Visu->handCursor);
		}

		if (Visu->mouse2Pressed) {
			Visu->mouse2EndDrag[0] = xpos/width*2-1;
			Visu->mouse2EndDrag[1] = ypos/height*1-0.5;


		double zoomFactor = 0.4;

		double scaleInc =  Visu->scale*zoomFactor*(Visu->mouse2EndDrag[1] - Visu->mouse2BeginDrag[1]);



		Visu->shift[0] += - 1*((Visu->mouse2BeginDrag[0]-Visu->shift[0])*scaleInc)/Visu->scale;
		Visu->shift[1] +=   1*((Visu->mouse2BeginDrag[1]+Visu->shift[1])*scaleInc)/Visu->scale;


		Visu->scale += scaleInc;
		}


		Visu->mouse2Pressed = true;

	}
	else {
		Visu->mouse1Pressed = false;
		Visu->mouse2Pressed = false;
		glfwSetCursor(window,NULL);
	}




}


/*
void Visu_SaveToImageFile(Visu* Visu) {

}
*/


