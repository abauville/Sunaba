/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"

#if (VISU)

void Visu_allocateMemory( Visu* Visu, Grid* Grid )
{
	Visu->U             = (GLfloat*)  malloc(2*Grid->nxEC*Grid->nyEC    * sizeof( GLfloat ));
	//Visu->vertices      = (GLfloat*)  malloc(Grid->nxS*Grid->nyS*2  * sizeof( GLfloat ));
	Visu->elements      = (GLuint*)   malloc(Visu->ntrivert    		* sizeof( GLuint ));

	Visu->vertices      = (GLfloat*)  malloc(4 * 4 * sizeof( GLfloat )); // 4 corners only
	Visu->particles 	= (GLfloat*) malloc (Visu->nParticles*4*sizeof(GLfloat));
	printf("%i  \n", (Visu->particleMeshRes+1) *3);
	Visu->particleMesh 	= (GLfloat*) malloc ((Visu->particleMeshRes+2) *3*sizeof(GLfloat));




	//Visu->elements      = (GLuint*)   malloc(6  * sizeof( GLuint  )); // 2 triangles

	Visu->nGlyphs 		= (int) ceil((double)Grid->nxS/(double)Visu->glyphSamplingRateX)*ceil((double)Grid->nyS/(double)Visu->glyphSamplingRateY);
	Visu->glyphs 		= (GLfloat*) malloc ( Visu->nGlyphs *4*sizeof(GLfloat));

	Visu->imageBuffer 	= (unsigned char*) malloc(Visu->retinaScale*Visu->retinaScale*4*Visu->width*Visu->height*sizeof(unsigned char)); // does not consider image resizing


	if (Visu->glyphMeshType==Triangle) {
		Visu->nGlyphMeshVert = 3;
	}
	else if (Visu->glyphMeshType==ThinArrow) {
		Visu->nGlyphMeshVert = 6;
	}
	if (Visu->glyphMeshType==ThickArrow) {
		Visu->nGlyphMeshVert = 18;
	}

	Visu->glyphMesh 	= (GLfloat*) malloc ( Visu->nGlyphMeshVert *2*sizeof(GLfloat));
}




void Visu_freeMemory( Visu* Visu )
{
	free(Visu->elements);
	free(Visu->U);
	free(Visu->particles);
	free(Visu->particleMesh);
	free(Visu->imageBuffer);

	free(Visu->glyphs);
	free(Visu->glyphMesh);

	glDeleteProgram(Visu->ShaderProgram);
	glDeleteProgram(Visu->ParticleShaderProgram);
	glDeleteProgram(Visu->ParticleBackgroundShaderProgram);
	glDeleteProgram(Visu->GlyphShaderProgram);
	glDeleteVertexArrays(1, &Visu->VAO );
	glDeleteVertexArrays(1, &Visu->VAO_part);
	glDeleteBuffers(1, &Visu->VBO);
	glDeleteBuffers(1, &Visu->VBO_part);
	glDeleteBuffers(1, &Visu->VBO_partMesh);
	glDeleteBuffers(1, &Visu->EBO);
	glDeleteTextures(1, &Visu->TEX);


	glDeleteVertexArrays(1,&Visu->VAO_glyph);
	glDeleteBuffers(1,&Visu->VBO_glyph);
	glDeleteBuffers(1,&Visu->VBO_glyphMesh);



}


void Visu_initWindow(Visu* Visu){

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
	Visu->window = glfwCreateWindow(Visu->width, Visu->height, "StokesFD", NULL, NULL);
	if (!Visu->window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(Visu->window);
	glfwSetKeyCallback(Visu->window, key_callback);

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
//#pragma omp parallel for private(iNode, thisParticle) schedule(static,32)
	FOR_PARTICLES
	Visu->particles[C] = thisParticle->x;
	Visu->particles[C+1] = thisParticle->y;

	if (Visu->typeParticles == PartPhase) {
		Visu->particles[C+2] = thisParticle->phase;
	} else if (Visu->typeParticles == PartTemp) {
#if (HEAT)
		Visu->particles[C+2] = thisParticle->T;
#else
		Visu->particles[C+2] = 0;
#endif
	} else if (Visu->typeParticles == PartSigma_xx) {
		Visu->particles[C+2] = thisParticle->sigma_xx_0;
	} else if (Visu->typeParticles == PartSigma_xy) {
		Visu->particles[C+2] = thisParticle->sigma_xy_0;
	}
	Visu->particles[C+3] = thisParticle->passive;

	C += 4;
	END_PARTICLES

}



void Visu_glyphs(Visu* Visu, Physics* Physics, Grid* Grid, Particles* Particles)
{

	int ix, iy;
	int C = 0;
	//PhaseFlag* Phase;
	if (Visu->glyphType == DarcyGradient) {
		/*
		Phase = (PhaseFlag*) malloc(Grid->nECTot * sizeof(PhaseFlag*));
		Darcy_setPhaseFlag(Phase, Darcy->hOcean, Grid, Particles);
		Darcy_setBC(Grid, Physics, Darcy->hOcean, Phase);
		 */
	}

	for (iy = 0; iy < Grid->nyS; iy+=Visu->glyphSamplingRateY) {
		for (ix = 0; ix < Grid->nxS; ix+=Visu->glyphSamplingRateX) {

			Visu->glyphs[C+0] = Grid->xmin + ix*Grid->dx;
			Visu->glyphs[C+1] = Grid->ymin + iy*Grid->dy;
			if (Visu->glyphType == StokesVelocity) {

				Visu->glyphs[C+2] = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2.0;
				Visu->glyphs[C+3] = (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2.0;
			}
			else if (Visu->glyphType == DarcyGradient) {
				/*
				if (Phase[ix+iy*Grid->nxEC] ==Air || Phase[ix+(iy+1)*Grid->nxEC]==Air) {
					Visu->glyphs[C+2] =  0;
					Visu->glyphs[C+3] =  0;
				} else {
					Visu->glyphs[C+2] =  - (Physics->psi[ix+1+iy*Grid->nxEC] - Physics->psi[ix+iy*Grid->nxEC])/Grid->dx;
					Visu->glyphs[C+3] =  - (Physics->psi[ix+(iy+1)*Grid->nxEC] - Physics->psi[ix+(iy)*Grid->nxEC]+Grid->dy)/Grid->dy;

					if (Phase[ix+(iy+1)*Grid->nxEC]==Air && Visu->glyphs[C+3]<0) {
						if (Darcy->rainFlux<-Visu->glyphs[C+3]) {
							Visu->glyphs[C+3] = -Darcy->rainFlux;
						}

					}


				}
				 */

				//printf("GradSouth = %.1e\n", (Physics->psi[ix   + (iy+1)*Grid->nxEC]-Physics->psi[ix+iy*Grid->nxEC]+dy)/dy );
			}
			else {
				printf("error: unknown glyphType\n");
				exit(0);
			}

			C+=4;
		}
	}

	if (Visu->glyphType == DarcyGradient) {
		//free(Phase);
	}






	/*
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
	 */
}


void Visu_particleMesh(Visu* Visu)
{


	// Create the particle mesh, i.e. cone
	compute radius = Visu->particleMeshSize;

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

void Visu_glyphMesh(Visu* Visu)
{
	GLfloat size 		= 0.1;
	GLfloat width 		= 0.1 	* size;
	GLfloat headLength 	= 0.35 	* size;
	GLfloat behindHead 	= 0.0 	* size;
	GLfloat stickWidth 	= 0.06 	* size;
	GLfloat headLineThickness = stickWidth*2;

	if (Visu->glyphMeshType==Triangle) {

		Visu->glyphMesh[0] = 0.0;
		Visu->glyphMesh[1] = -width;
		Visu->glyphMesh[2] = 0.0;
		Visu->glyphMesh[3] = +width;
		Visu->glyphMesh[4] = size;
		Visu->glyphMesh[5] = 0.0;
	}
	else if (Visu->glyphMeshType==ThinArrow) {
		Visu->glyphMesh[0] = 0.0;
		Visu->glyphMesh[1] = 0.0;
		Visu->glyphMesh[2] = size-headLength;
		Visu->glyphMesh[3] = 0.0;
		Visu->glyphMesh[4] = size-headLength-behindHead;
		Visu->glyphMesh[5] = +width;
		Visu->glyphMesh[6] = size;
		Visu->glyphMesh[7] = 0.0;
		Visu->glyphMesh[8] = size-headLength-behindHead;
		Visu->glyphMesh[9] = -width;
		Visu->glyphMesh[10] = size-headLength;
		Visu->glyphMesh[11] = 0.0;
	}
	else if (Visu->glyphMeshType==ThickArrow) {

		// Stick, triangle 1
		Visu->glyphMesh[ 0] = 0.0;
		Visu->glyphMesh[ 1] = stickWidth;
		Visu->glyphMesh[ 2] = 0.0;
		Visu->glyphMesh[ 3] = -stickWidth;
		Visu->glyphMesh[ 4] = size-headLength;
		Visu->glyphMesh[ 5] = -stickWidth;

		// Stick, triangle 2
		Visu->glyphMesh[ 6] = size-headLength;
		Visu->glyphMesh[ 7] = -stickWidth;
		Visu->glyphMesh[ 8] = size-headLength;
		Visu->glyphMesh[ 9] = +stickWidth;
		Visu->glyphMesh[10] = 0.0;
		Visu->glyphMesh[11] = stickWidth;

		/*
		// Lower back point
		Visu->glyphMesh[12] = size-headLength;
		Visu->glyphMesh[13] = -stickWidth;
		Visu->glyphMesh[14] = size-headLength-headLineThickness;
		Visu->glyphMesh[15] = -stickWidth;
		Visu->glyphMesh[16] = size-headLength-behindHead;
		Visu->glyphMesh[17] = -stickWidth-width;

		// Upper back point
		Visu->glyphMesh[18] = size-headLength;
		Visu->glyphMesh[19] = stickWidth;
		Visu->glyphMesh[20] = size-headLength-headLineThickness;
		Visu->glyphMesh[21] = stickWidth;
		Visu->glyphMesh[22] = size-headLength-behindHead;
		Visu->glyphMesh[23] = stickWidth+width;
		 */
		// Lower point, 1
		double alpha = atan(((double)width+(double)stickWidth)/((double)behindHead+(double)headLength));
		double gamma = PI/2.0 - 2.0*fabs(alpha);
		double L = headLineThickness/sin(alpha);
		double x = -L*sin(gamma);
		double y =  L*cos(gamma);
		Visu->glyphMesh[12] = size-headLength-behindHead;
		Visu->glyphMesh[13] = -stickWidth-width;
		Visu->glyphMesh[14] = size;
		Visu->glyphMesh[15] = 0.0;
		Visu->glyphMesh[16] = size-L;
		Visu->glyphMesh[17] = 0.0;
		// Lower point, 2


		Visu->glyphMesh[18] = size-L;
		Visu->glyphMesh[19] = 0.0;
		Visu->glyphMesh[20] = size-headLength-behindHead;
		Visu->glyphMesh[21] = -stickWidth-width;
		Visu->glyphMesh[22] = size-headLength-behindHead-x;
		Visu->glyphMesh[23] = -stickWidth-width+y;

		// Upper point
		Visu->glyphMesh[24] = size-headLength-behindHead;
		Visu->glyphMesh[25] = +stickWidth+width;
		Visu->glyphMesh[26] = size;
		Visu->glyphMesh[27] = 0.0;
		Visu->glyphMesh[28] = size-L;
		Visu->glyphMesh[29] = 0.0;

		// Upper point, 2
		Visu->glyphMesh[30] = size-L;
		Visu->glyphMesh[31] = 0.0;
		Visu->glyphMesh[32] = size-headLength-behindHead;
		Visu->glyphMesh[33] = +stickWidth+width;
		Visu->glyphMesh[34] = size-headLength-behindHead-x;
		Visu->glyphMesh[35] = +stickWidth+width-y;



	}
	else {
		printf("error: unknown Visu->glyphMeshType");
	}



}





void Visu_init(Visu* Visu, Grid* Grid, Particles* Particles, Char* Char)
{


	// Non dimensionalization
	Visu->particleMeshSize /= Char->length;

	Visu_initWindow(Visu);

	// Create the surface for plotting grid data and fill the Particles
	Visu->elements[0] = 0;
	Visu->elements[1] = 1;
	Visu->elements[2] = 2;
	Visu->elements[3] = 2;
	Visu->elements[4] = 3;
	Visu->elements[5] = 1;

	Visu_updateVertices(Visu, Grid);
	Visu_particles(Visu, Particles, Grid);


	///Init shader
	// =======================================
	Visu->VertexShaderFile 				= "src/Shaders/shader.vs";
	Visu->FragmentShaderFile 			= "src/Shaders/shader.fs";
	Visu->ParticleVertexShaderFile 		= "src/Shaders/particleShader.vs";
	Visu->ParticleGeometryShaderFile 	= "src/Shaders/particleShader.gs";
	Visu->ParticleFragmentShaderFile 	= "src/Shaders/particleShader.fs";
	Visu->ParticleBackgroundVertexShaderFile 	= "src/Shaders/particleBackgroundShader.vs";
	Visu->ParticleBackgroundFragmentShaderFile 	= "src/Shaders/particleBackgroundShader.fs";

	Visu->GlyphVertexShaderFile 			= "src/Shaders/glyphShader.vs";
	Visu->GlyphFragmentShaderFile 			= "src/Shaders/glyphShader.fs";


	Visu->ShaderProgram = 0;
	Visu->ParticleShaderProgram = 0;
	Visu->ParticleBackgroundShaderProgram = 0;
	Visu->GlyphShaderProgram = 0;

	// Generate reference to objects (indexes that act as pointers to graphic memory)
	// =======================================
	Visu->VAO = 0; // Reference to the Vertex   array object
	Visu->VBO = 0; // Reference to the Vertex   buffer object
	Visu->EBO = 0; // Reference to the Element  buffer object
	Visu->TEX = 0; // Reference to the Element  buffer object



	// And assigned them to objects (stored in the graphic memory)
	// =======================================
	glGenVertexArrays(1, &Visu->VAO);

	glGenBuffers(1, &Visu->VBO);
	glGenBuffers(1, &Visu->EBO);
	glGenTextures(1, &Visu->TEX);

	glGenVertexArrays(1, &Visu->VAO_part);
	glGenBuffers(1, &Visu->VBO_part);
	glGenBuffers(1, &Visu->VBO_partMesh);

	glGenVertexArrays(1, &Visu->VAO_glyph);
	glGenBuffers(1, &Visu->VBO_glyph);
	glGenBuffers(1, &Visu->VBO_glyphMesh);



	// Bind Vertex Array object
	// =======================================
	glBindVertexArray(Visu->VAO);
	// compile shaders
	// =======================================
	const char* dumShaderFile = NULL;
	compileShaders(&Visu->ShaderProgram, Visu->VertexShaderFile, Visu->FragmentShaderFile, dumShaderFile, false);
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
	if (Visu->filter == Linear) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	} else if (Visu->filter == Nearest) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	}


	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG, Grid->nxEC, Grid->nyEC, 0, GL_RG, GL_FLOAT, Visu->U);
	glBindTexture(GL_TEXTURE_2D, 0);










	// Declare the initial values of uniforms
	// =======================================
	int width, height;
	glfwGetWindowSize(Visu->window, &width, &height);
	GLfloat ratio = (GLfloat)width/(GLfloat)height;

	if ((Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0])>(Grid->ymax-Grid->ymin)*(1+2*Visu->shiftFac[1])/ratio){
		Visu->scale = 2.0/(1.05*(Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0]));
	}
	else {
		Visu->scale = 2.0/(1.05*(Grid->ymax-Grid->ymin)*(1+2*Visu->shiftFac[1])*ratio);
	}

	// Visu->scale = 2.0/(0.85*(Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0]));


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
	compileShaders(&Visu->ParticleBackgroundShaderProgram, Visu->ParticleBackgroundVertexShaderFile, Visu->ParticleBackgroundFragmentShaderFile, dumShaderFile, false);
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



	// =======================================
	// Glyph
	// =======================================
	glBindVertexArray(Visu->VAO_glyph);

	compileShaders(&Visu->GlyphShaderProgram, Visu->GlyphVertexShaderFile, Visu->GlyphFragmentShaderFile, dumShaderFile, false);
	glUseProgram(Visu->GlyphShaderProgram);

	GLint GlyphVertAttrib    	= glGetAttribLocation(Visu->GlyphShaderProgram,"glyphVertex");
	GLint GlyphData    	 		= glGetAttribLocation(Visu->GlyphShaderProgram,"dataVector");
	GLint GlyphMeshVertAttrib  	= glGetAttribLocation(Visu->GlyphShaderProgram,"glyphMeshVertex");

	glBindVertexArray(0);




	// Create Mesh

	glBindVertexArray(Visu->VAO_glyph);
	Visu_glyphMesh(Visu);



	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_glyphMesh);

	glBufferData(GL_ARRAY_BUFFER, 2*Visu->nGlyphMeshVert*sizeof(GLfloat), Visu->glyphMesh, GL_STATIC_DRAW);
	glVertexAttribPointer(GlyphMeshVertAttrib , 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), 0);
	glEnableVertexAttribArray(GlyphMeshVertAttrib);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_glyph);
	glBufferData(GL_ARRAY_BUFFER, 4*Visu->nGlyphs*sizeof(GLfloat), NULL, GL_STREAM_DRAW);
	glVertexAttribPointer(GlyphVertAttrib , 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), 0);
	glEnableVertexAttribArray(GlyphVertAttrib );
	glVertexAttribPointer(GlyphData , 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), (void*)(2*sizeof(GLfloat)));
	glEnableVertexAttribArray(GlyphData );
	glBindBuffer(GL_ARRAY_BUFFER, 0);


	glVertexAttribDivisor(GlyphMeshVertAttrib , 0); // never changes
	glVertexAttribDivisor(GlyphVertAttrib, 1); // counter of +1 per instance
	glVertexAttribDivisor(GlyphData, 1);

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
	compute xmin = Grid->xmin;//-0.5*Grid->dx;
	compute ymin = Grid->ymin;//-0.5*Grid->dy;

	compute Ratio = (Grid->xmax-Grid->xmin)/(Grid->ymax-Grid->ymin);

	int ix, iy;
	int C = 0;
	compute signX[2] = {1.0,-1.0};
	compute signY[2] = {1.0,-1.0};
	for (iy = 0; iy < 2; ++iy) {
		for (ix = 0; ix < 2; ++ix) {
			Visu->vertices[C  ] = xmin + ix*(Grid->xmax-xmin) ;
			Visu->vertices[C+1] = ymin + iy*(Grid->ymax-ymin);

			// Showing the sides row and columns
			//Visu->vertices[C+2] = 1.0*ix;
			//Visu->vertices[C+3] = 1.0*iy;

			// Without showing the sides row and column
			Visu->vertices[C+2] = 1.0*ix+signX[ix]*((float)Grid->nxC/(float)Grid->nxEC)*Grid->dx/(Grid->xmax-xmin);
			Visu->vertices[C+3] = 1.0*iy+signY[iy]*((float)Grid->nyC/(float)Grid->nyEC)*Grid->dy/(Grid->ymax-ymin);

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
#pragma omp parallel for private(iy, ix, I, iNW, iNE, iSW, iSE) schedule(static,32)
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			I = 2* (ix + iy*Grid->nxEC);
			Visu->U[I] = CellValue[ix + iy*Grid->nxEC];//(CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}
}

void Visu_updateCenterValuei(Visu* Visu, Grid* Grid, int* CellValue, int BCType)
{
	// UC is a scalar CellValue defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	int iNW, iNE, iSW, iSE;
	// CellValue interpolated on the center nodes
	// ======================================
#pragma omp parallel for private(iy, ix, I, iNW, iNE, iSW, iSE) schedule(static,32)
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			I = 2* (ix + iy*Grid->nxEC);
			Visu->U[I] = CellValue[ix + iy*Grid->nxEC];//(CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
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
	/*
#pragma omp parallel for private(iy, ix, I, A, B) schedule(static,32)
	for (iy=0; iy<Grid->nyEC; iy++){
		for (ix=0; ix<Grid->nxEC; ix++) {
			I = 2*(ix+iy*Grid->nxEC);
			Visu->U[I] = 0;
			//Visu->U[I]  = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2;
			//Visu->U[I] += (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2;
		}
		//printf("\n");
	}
	*/


#pragma omp parallel for private(iy, ix, I, A, B) schedule(static,32)
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = 2*(ix+iy*Grid->nxEC);


			A  = (Physics->Vx[ix-1  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy  )*Grid->nxVx])/2.0;
			B  = 0.0;//(Physics->Vy[ix    +(iy-1)*Grid->nxVy] + Physics->Vy[ix  +(iy  )*Grid->nxVy])/2.0;
			Visu->U[I] = sqrt(A*A + B*B);



		}
	}
}


void Visu_divV(Visu* Visu, Grid* Grid, Physics* Physics) {
	int iy, ix;
	int I = 0;

	compute dx, dy, divV;

//#pragma omp parallel for private(iy, ix, I, dx, dy, divV) schedule(static,32)
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = 2*(ix+iy*Grid->nxEC);


			dx = Grid->DXS[ix-1];
			dy = Grid->DYS[iy-1];
			divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/dx;
			divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/dy;

			Visu->U[I] = divV;

			//printf("divV[%i, %i] = %.2e, dy = %.2e, dx = %.2e\n",iy, ix, divV, dy, dx);

		}
	}

}


void Visu_stress(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC)
{

	int iy, ix;
	int I = 0;
	//compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	compute sigma_xy;
	int nxS = Grid->nxS;
	//Visu_updateCenterValue (Visu, Grid, Physics->sigma_xx_0, BC->SetupType);
#pragma omp parallel for private(iy, ix, I, sigma_xy) schedule(static,32)
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = (ix+iy*Grid->nxEC);

			sigma_xy = 0.25* ( Physics->sigma_xy_0[(ix-1)+(iy-1)*nxS] + Physics->sigma_xy_0[(ix-1)+(iy)*nxS] + Physics->sigma_xy_0[(ix)+(iy)*nxS] + Physics->sigma_xy_0[(ix)+(iy-1)*nxS] );
			// second invariant
			Visu->U[2*I] = sqrt( Physics->sigma_xx_0[I]*Physics->sigma_xx_0[I] + sigma_xy*sigma_xy );

			//Visu->U[2*I] = sqrt(  Physics->sigma_xy_0[I]*Physics->sigma_xy_0[I]   +   Physics->sigma_xx_0[I]*Physics->sigma_xx_0[I]  );
		}
		//printf("\n");
	}


	// Replace boundary values by their neighbours
	int INeigh;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {

		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];

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


void Visu_SIIOvYield(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC, Numerics* Numerics) {
	Visu_stress(Visu, Grid, Physics, BC);
	compute sigma_xy, sigma_xx, sigmaII;


	compute sigma_y, Pe;
	compute phi = 0.0;
	int iCell;
	compute phiCrit = Numerics->phiCrit;
	int ix, iy;
	//printf("=== Check sigmaII grid  ===\n");
	for (iy=1; iy<Grid->nyEC-1; ++iy) {
		for (ix=1; ix<Grid->nxEC-1; ++ix) {
			iCell = ix+iy*Grid->nxEC;
#if (DARCY)

			phi = Physics->phi[iCell];
			if (phi>=phiCrit) {
				Pe 		= Physics->Pc[iCell];
			} else {
				Pe 		= Physics->P [iCell];
			}

#else
			Pe = Physics->P[iCell];
#endif

			/*
			sigma_xy  = Physics->sigma_xy_0[ix-1 + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix   + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix-1 + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
			sigma_xy += Physics->sigma_xy_0[ix   + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
			sigma_xy /= 4.0;

			sigma_xx = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];



			// Get invariants EII and SigmaII
			//Physics_computeStrainInvariantForOneCell(Physics, Grid, ix,iy, &EII);
			sigmaII = (1.0-phi) * sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);
			//sigmaII = (1.0-phi) * sigma_xx;

			*/



			sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Pe * sin(Physics->frictionAngle[iCell]);


			Physics_computeStressInvariantForOneCell(Physics, Grid, ix, iy, &sigmaII);






			//
			Visu->U[2*iCell] = sigmaII/sigma_y;
			//Visu->U[2*iCell] = Visu->U[2*iCell]/sigma_y;
			//printf("%.2e  ", Pe);
		}
		//printf("\n");
	}

// Replace boundary values by their neighbours
	int INeigh, I;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {

		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];

	}


}



void Visu_PeOvYield(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC, Numerics* Numerics) {
	Visu_stress(Visu, Grid, Physics, BC);
	compute sigma_xy, sigma_xx, sigmaII;



	compute Py, Pe;
	compute phi = 0.0;
	compute sigmaT;

	compute R = 2.0;

	int iCell;

	int ix, iy;
	for (iy=1; iy<Grid->nyEC-1; ++iy) {
		for (ix=1; ix<Grid->nxEC-1; ++ix) {
			iCell = ix+iy*Grid->nxEC;
#if (DARCY)
			compute phiCrit = Numerics->phiCrit;
			phi = Physics->phi[iCell];
			if (phi>=phiCrit) {
				Pe 		= Physics->Pc[iCell];
			} else {
				Pe 		= Physics->P [iCell];
			}
			//Pe 		= Physics->Pc[iCell];
#else
			Pe = Physics->P[iCell];
#endif


			sigmaT = Physics->cohesion[iCell]/R;

			Physics_computeStressInvariantForOneCell(Physics, Grid, ix, iy, &sigmaII);

			Py = sigmaII - sigmaT;




			//sigma_y = Physics->cohesion[iCell] * cos(Physics->frictionAngle[iCell])   +   Pe * sin(Physics->frictionAngle[iCell]);





			Visu->U[2*iCell] = (Pe-Py)/Py;


			/*
			if (phi>=phiCrit) {
				Visu->U[2*iCell] = 1.0;
			} else {
				Visu->U[2*iCell] = 0.0;
			}
			*/

			//Visu->U[2*iCell] = Visu->U[2*iCell]/sigma_y;

		}
	}

// Replace boundary values by their neighbours
	int INeigh, I;
	// lower boundary
	iy = 0;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy+1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy+1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy+1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}




	// upper boundary
	iy = Grid->nyEC-1;
	for (ix = 0; ix<Grid->nxEC; ix++) {
		I = ix + iy*Grid->nxEC;
		if (ix==0) {
			INeigh =   ix+1 + (iy-1)*Grid->nxEC  ;
		} else if (ix==Grid->nxEC-1) {
			INeigh =   ix-1 + (iy-1)*Grid->nxEC  ;
		} else {
			INeigh =   ix + (iy-1)*Grid->nxEC  ;
		}
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// left boundary
	ix = 0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {

		I = ix + iy*Grid->nxEC;
		INeigh =   ix+1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];
	}
	// right boundary
	ix = Grid->nxEC-1;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		I = ix + iy*Grid->nxEC;
		INeigh =   ix-1 + (iy)*Grid->nxEC  ;
		Visu->U[2*I] = Visu->U[2*INeigh];

	}


}













void Visu_alphaValue(Visu* Visu, Grid* Grid, Physics* Physics) {
	// Based on phase
	//compute y, depth;
	//compute hOcean = Grid->ymin + (Grid->ymax-Grid->ymin)*0.35;

	float alpha;
	/*
	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, alpha) schedule(static,32)
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		thisParticle = Particles->linkHead[iNode];
		alpha = 1.0;
		while (thisParticle != NULL && alpha >0) {
			if (thisParticle->phase==0) {
				alpha = 0.0;
			}
			else if (thisParticle->phase==1 && Visu->type != FluidPressure) {
				alpha = 0.0;
			}
			thisParticle = thisParticle->next;
		}

		Visu->U[2*iNode+1] = alpha;

	}
	*/
	int i;
	for (i = 0; i < Grid->nECTot; ++i) {
		Visu->U[2*i+1] = 1.0;
		if ( Physics->phase[i] == Physics->phaseAir || Physics->phase[i] == Physics->phaseAir ) {
			Visu->U[2*i+1] = 0.0;
		}
	}

}








void Visu_updateUniforms(Visu* Visu)
{
	int width, height;
	glfwGetWindowSize(Visu->window, &width, &height);
	GLfloat ratio = (GLfloat)width/(GLfloat)height;
	//printf("ratio = %.2f, scale = %.2f\n\n\n\n",ratio, Visu->scale);

	GLfloat Transform[] = {Visu->scale,0.0f,0.0f,0.0f , 0.0f,Visu->scale*ratio,0.0f,0.0f , 0.0f,0.0f,1.0f,0.0f , Visu->shift[0],Visu->shift[1]*ratio,Visu->shift[2],1.0f};
	GLuint loc = glGetUniformLocation(Visu->ShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->ParticleBackgroundShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->GlyphShaderProgram, "transform");
	glUniformMatrix4fv(loc, 1, GL_FALSE, &Transform[0]);

	loc = glGetUniformLocation(Visu->GlyphShaderProgram, "glyphScale");
	glUniform1f(loc, Visu->glyphScale);



	loc = glGetUniformLocation(Visu->ParticleShaderProgram, "size");
	//printf("scale: %.3f\n",Visu->scale);
	glUniform1f(loc, 1.0*Visu->scale);


	//printf("scale: %.3f\n",Visu->scale);
	GLfloat type;
	if (Visu->typeParticles == PartPhase ) {
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

	loc = glGetUniformLocation(Visu->ShaderProgram, "alphaOnValue");
	glUniform1i(loc, Visu->alphaOnValue);

}




void Visu_update(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC, Char* Char, MatProps* MatProps, EqSystem* EqStokes, EqSystem* EqThermal, Numbering* NumStokes, Numbering* NumThermal, Numerics* Numerics)
{
	int i;
	char title[1024];
	switch (Visu->type) {
	case Viscosity:
		glfwSetWindowTitle(Visu->window, "Viscosity");
		Visu->valueScale = 1.0;//MatProps->eta0[0];//Char->viscosity;
		Visu->valueShift = 0;
		Visu_updateCenterValue(Visu, Grid, Physics->eta, BC->SetupType);


		Visu->colorScale[0] = -4.0;
		Visu->colorScale[1] =  4.0;
		Visu->log10_on = true;
		break;
	case Khi:
		glfwSetWindowTitle(Visu->window, "Khi");
		Visu->valueScale = 1.0;//MatProps->eta0[0];//Char->viscosity;
		Visu->valueShift = 0;
		Visu_updateCenterValue(Visu, Grid, Physics->khi, BC->SetupType);


		Visu->colorScale[0] = -4.0;
		Visu->colorScale[1] =  4.0;
		Visu->log10_on = true;
		break;

	case Khib:
#if (DARCY)
		glfwSetWindowTitle(Visu->window, "Khi_b");
		Visu->valueScale = 1.0;//MatProps->eta0[0];//Char->viscosity;
		Visu->valueShift = 0;
		Visu_updateCenterValue(Visu, Grid, Physics->khi_b, BC->SetupType);


		Visu->colorScale[0] = -4.0;
		Visu->colorScale[1] =  4.0;
		Visu->log10_on = true;
		break;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
		break;
#endif

	case StrainRate:
		glfwSetWindowTitle(Visu->window, "StrainRate");
		Visu->valueScale = Physics->epsRef;
		Visu->valueShift = 0.0;
		Visu_strainRate(Visu, Grid, Physics, BC);

		Visu->colorScale[0] = -1.0;
		Visu->colorScale[1] =  1.0;
		Visu->log10_on = true;
		break;
	case Stress:
		glfwSetWindowTitle(Visu->window, "Stress");
		Visu->valueScale = 1.0;
		Visu->valueShift = 0.0;
		Visu_stress(Visu, Grid, Physics, BC);

		Visu->colorScale[0] = -20.0;
		Visu->colorScale[1] =  20.0;
		Visu->log10_on = false;
		break;
	case Velocity:
		glfwSetWindowTitle(Visu->window, "Velocity");
		Visu_velocity(Visu, Grid, Physics);
		Visu->valueScale = 0.2*Physics->maxV;//(Physics->epsRef*Grid->xmax);

		sprintf(title,"Velocity, scale = %.2e",Visu->valueScale);
		glfwSetWindowTitle(Visu->window, title);
		Visu->valueShift = 0;
		Visu->colorScale[0] = -1;
		Visu->colorScale[1] =  1;
		Visu->log10_on = true;
		break;
	case VelocityDiv:
		glfwSetWindowTitle(Visu->window, "Velocity divergence, /!\\ values are computed using the updated dx, dy (i.e. values appear much larger)");
		Visu_divV(Visu, Grid, Physics);
		Visu->valueScale = 1e-6;//(Physics->epsRef*Grid->xmax);
		Visu->valueShift = 0;
		Visu->colorScale[0] = -4.;
		Visu->colorScale[1] =  4.;
		Visu->log10_on = true;
		break;
	case SIIOvYield:
		glfwSetWindowTitle(Visu->window, "Stress_II/Stress_y");
		Visu_SIIOvYield(Visu, Grid, Physics, BC, Numerics);
		Visu->valueScale = 1.0;//(Physics->epsRef*Grid->xmax);
		Visu->valueShift = -1.0;
		Visu->colorScale[0] = -0.5;
		Visu->colorScale[1] =  0.5;
		Visu->log10_on = false;
		break;

	case PeOvYield:
		glfwSetWindowTitle(Visu->window, "Pe/Py");
		Visu_PeOvYield(Visu, Grid, Physics, BC, Numerics);
		Visu->valueScale = 1.0;//(Physics->epsRef*Grid->xmax);
		Visu->valueShift = -1.0;
		Visu->colorScale[0] = -1.0;
		Visu->colorScale[1] =  1.0;
		Visu->log10_on = false;
		break;

	case Pressure:
		glfwSetWindowTitle(Visu->window, "Pressure");
		Visu_updateCenterValue(Visu, Grid, Physics->P, BC->SetupType);

		Visu->valueScale = 1.0;//Char->stress;
		Visu->valueShift = 0;
		Visu->colorScale[0] = -2.0;
		Visu->colorScale[1] =  2.0;
		Visu->log10_on = false;
		break;
	case Density:
		glfwSetWindowTitle(Visu->window, "Density*g");
		Visu_updateCenterValue(Visu, Grid, Physics->rho_g, BC->SetupType);
		Visu->valueScale = MatProps->rho0_g[0];
		Visu->valueShift = 0;
		Visu->colorScale[0] = -0.0002;
		Visu->colorScale[1] =  0.0002;
		Visu->log10_on = true;
		break;
	case Temperature:
#if (HEAT)
		glfwSetWindowTitle(Visu->window, "Temperature");
		Visu_updateCenterValue(Visu, Grid, Physics->T, BC->SetupType); // Not optimal but good enough for the moment
		Visu->valueScale = 1.0;

		Visu->colorScale[0] = -1.0;
		Visu->colorScale[1] =  1.0;
		Visu->valueShift = 0;//1*Visu->colorScale[0];
		Visu->log10_on = false;
#else
		glfwSetWindowTitle(Visu->window, "Temperature is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif

		break;
	case FluidPressure:

			glfwSetWindowTitle(Visu->window, "Fluid pressure");
#if (DARCY)

			//printf("Visu Psi[0] = %.1e\n", Physics->psi[0]);
			Visu_updateCenterValue(Visu, Grid, Physics->Pf, BC->SetupType); // Not optimal but good enough for the moment
			//free(dum);
			Visu->valueScale = 1.0;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nSTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
			Visu->colorScale[0] = -1.0;
			Visu->colorScale[1] =  1.0;
			Visu->valueShift = 0.0*Visu->colorScale[0];
			Visu->log10_on = false;


		break;
	case CompactionPressure:
		glfwSetWindowTitle(Visu->window, "Compaction pressure");
#if (DARCY)

			//printf("Visu Psi[0] = %.1e\n", Physics->psi[0]);
			Visu_updateCenterValue(Visu, Grid, Physics->Pc, BC->SetupType); // Not optimal but good enough for the moment
			//free(dum);
			Visu->valueScale = 0.1;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
			Visu->colorScale[0] = -1.0;
			Visu->colorScale[1] =  1.0;
			Visu->valueShift = 0.0*Visu->colorScale[0];
			Visu->log10_on = false;
			break;

	case Permeability:
		glfwSetWindowTitle(Visu->window, "Permeability/eta_f");
#if (DARCY)

			//printf("Visu Psi[0] = %.1e\n", Physics->psi[0]);
			Visu_updateCenterValue(Visu, Grid, Physics->perm_eta_f, BC->SetupType); // Not optimal but good enough for the moment
			//free(dum);
			Visu->valueScale = Physics->perm_eta_f[0];
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
			Visu->colorScale[0] = -1.0;
			Visu->colorScale[1] =  1.0;
			Visu->valueShift = 0.0*Visu->colorScale[0];
			Visu->log10_on = true;


		break;




	case Porosity:
		glfwSetWindowTitle(Visu->window, "Porosity");
#if (DARCY)

			//printf("Visu Psi[0] = %.1e\n", Physics->psi[0]);
			Visu_updateCenterValue(Visu, Grid, Physics->phi, BC->SetupType); // Not optimal but good enough for the moment
			//free(dum);
			Visu->valueScale = 1.0;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
			Visu->colorScale[0] = -0.01;
			Visu->colorScale[1] =  0.01;
			Visu->valueShift = -0.1;//0.0*Visu->colorScale[0];
			Visu->log10_on = false;


		break;
	case Phase:
		glfwSetWindowTitle(Visu->window, "Phase");
		Visu_updateCenterValuei(Visu, Grid, Physics->phase, BC->SetupType);

		Visu->valueScale = 1.0;//Char->stress;
		Visu->valueShift = -2.0;
		Visu->colorScale[0] = -2.;
		Visu->colorScale[1] =  2.;
		Visu->log10_on = false;
		break;


	case VxRes:
	case VyRes:
	case PRes:
	case PfRes:
	case PcRes:

		if 		 (Visu->type==VxRes) {
			glfwSetWindowTitle(Visu->window, "Vx residual");
		} else if(Visu->type==VyRes) {
			glfwSetWindowTitle(Visu->window, "Vy residual");
		} else if(Visu->type==PRes) {
			glfwSetWindowTitle(Visu->window, "P residual");
		} else if(Visu->type==PfRes) {
			glfwSetWindowTitle(Visu->window, "Pf residual");
		} else if(Visu->type==PcRes) {
			glfwSetWindowTitle(Visu->window, "Pc Residual");
		}

		Visu_residual(Visu, Grid, EqStokes, NumStokes);

		Visu->valueScale = 1e-7;
		Visu->valueShift = 0.0;
		Visu->colorScale[0] = -2.;
		Visu->colorScale[1] =  2.;
		Visu->log10_on = true;
		break;

	case TRes:
		glfwSetWindowTitle(Visu->window, "T residual");
		Visu_residual(Visu, Grid, EqThermal, NumThermal);

		Visu->valueScale = 1e-4;
		Visu->valueShift = 0.0;
		Visu->colorScale[0] = -2.;
		Visu->colorScale[1] =  2.;
		Visu->log10_on = true;
		break;

	case Blank:
		glfwSetWindowTitle(Visu->window, "Blank");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
		Visu->valueScale = 1.0;
		Visu->colorScale[0] = -1;
		Visu->colorScale[1] =  1;
		Visu->valueShift = 0;
		Visu->log10_on = false;
		break;
	default:
		printf("Error: unknown Visu->type: %i",Visu->type);
	}

	switch (Visu->typeParticles) {
	case PartPhase:
		Visu->partColorScale[0] = -3;
		Visu->partColorScale[1] =  3;
		break;
	case PartTemp:
		Visu->partColorScale[0] = -1.0;
		Visu->partColorScale[1] =  1.0;
		break;
	case PartSigma_xx:
		Visu->partColorScale[0] = -0.25;
		Visu->partColorScale[1] =  0.25;
		break;
	case PartSigma_xy:
		Visu->partColorScale[0] = -0.25;
		Visu->partColorScale[1] =  0.25;
		break;
	default:
		printf("Error: unknown Visu->typeParticles: %i",Visu->typeParticles);
	}




	Visu_updateUniforms(Visu);
}


void Visu_checkInput(Visu* Visu)
{
	// Check keyboard events
	if (glfwGetKey(Visu->window, GLFW_KEY_1) == GLFW_PRESS) {
		Visu->type = Viscosity;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_2) == GLFW_PRESS) {
		Visu->type = StrainRate;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_3) == GLFW_PRESS) {
		Visu->type = Stress;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_4) == GLFW_PRESS) {
		Visu->type = Pressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_5) == GLFW_PRESS) {
		Visu->type = Density;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_6) == GLFW_PRESS) {
		Visu->type = Temperature;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_7) == GLFW_PRESS) {
		Visu->type = Velocity;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_8) == GLFW_PRESS) {
		Visu->type = FluidPressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_9) == GLFW_PRESS) {
		Visu->type = Permeability;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_0) == GLFW_PRESS) {
		Visu->type = Porosity;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_C) == GLFW_PRESS) {
		Visu->type = CompactionPressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_V) == GLFW_PRESS) {
		Visu->type = Phase;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_Y) == GLFW_PRESS) {
		Visu->type = VelocityDiv;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_U) == GLFW_PRESS) {
		Visu->type = SIIOvYield;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_I) == GLFW_PRESS) {
		Visu->type = PeOvYield;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_A) == GLFW_PRESS) {
		Visu->type = Khi;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_S) == GLFW_PRESS) {
		Visu->type = Khib;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_M) == GLFW_PRESS) {
		Visu->type = Blank;
		Visu->update = true;
	}



	// Residuals
#if (HEAT)
	else if (glfwGetKey(Visu->window, GLFW_KEY_D) == GLFW_PRESS) {
		Visu->type = TRes;
		Visu->update = true;
	}
#endif
	else if (glfwGetKey(Visu->window, GLFW_KEY_F) == GLFW_PRESS) {
		Visu->type = VxRes;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_G) == GLFW_PRESS) {
		Visu->type = VyRes;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_H) == GLFW_PRESS) {
#if (DARCY)
		Visu->type = PfRes;
#else
		Visu->type = PRes;
#endif
		Visu->update = true;
	}
#if (DARCY)
	else if (glfwGetKey(Visu->window, GLFW_KEY_J) == GLFW_PRESS) {
		Visu->type = PcRes;
		Visu->update = true;
	}
#endif




	// Particles
	else if (glfwGetKey(Visu->window, GLFW_KEY_Q) == GLFW_PRESS) {
		Visu->typeParticles = PartPhase;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_W) == GLFW_PRESS) {
		Visu->typeParticles = PartTemp;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_E) == GLFW_PRESS) {
		Visu->typeParticles = PartSigma_xx;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_R) == GLFW_PRESS) {
		Visu->typeParticles = PartSigma_xy;
		Visu->update = true;
	}





	else if (glfwGetKey(Visu->window, GLFW_KEY_P) == GLFW_PRESS) {
		Visu->paused = true;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_O) == GLFW_PRESS) {
		Visu->paused = false;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_L) == GLFW_PRESS) {
		Visu->showParticles = true;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_K) == GLFW_PRESS) {
		Visu->showParticles = false;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		Visu->initPassivePart = true;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_X) == GLFW_PRESS) {
		Visu->transparency = true;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_Z) == GLFW_PRESS) {
		Visu->transparency = false;
		Visu->update = true;
	}



	// Check mouse events

	// Left click - shift visu
	if (glfwGetMouseButton(Visu->window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
		Visu->update = true;
		double xpos, ypos;
		glfwGetCursorPos(Visu->window, &xpos, &ypos);
		if (!Visu->mouse1Pressed) {
			Visu->mouse1BeginDrag[0] = xpos;
			Visu->mouse1BeginDrag[1] = ypos;
			glfwSetCursor(Visu->window,Visu->handCursor);
		}
		if (Visu->mouse1Pressed) {
			Visu->mouse1EndDrag[0] = xpos;
			Visu->mouse1EndDrag[1] = ypos;

			int width, height;
			glfwGetWindowSize(Visu->window, &width, &height);
			Visu->shift[0] += (Visu->mouse1EndDrag[0] - Visu->mouse1BeginDrag[0])/width*2.0;
			Visu->shift[1] -= (Visu->mouse1EndDrag[1] - Visu->mouse1BeginDrag[1])/height*1.0;

			Visu->mouse1BeginDrag[0] = xpos;
			Visu->mouse1BeginDrag[1] = ypos;

		}
		Visu->mouse1Pressed = true;
	}

	// Righr click - zoom
	else if (glfwGetMouseButton(Visu->window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS){
		Visu->update = true;
		double xpos, ypos;
		glfwGetCursorPos(Visu->window, &xpos, &ypos);
		int width, height;
		glfwGetWindowSize(Visu->window, &width, &height);

		if (!Visu->mouse2Pressed) {
			Visu->mouse2BeginDrag[0] = xpos/width*2.-1.;
			Visu->mouse2BeginDrag[1] = ypos/height*2.-1.;
			glfwSetCursor(Visu->window,Visu->handCursor);
		}

		if (Visu->mouse2Pressed) {
			Visu->mouse2EndDrag[0] = xpos/width*2.-1.;
			Visu->mouse2EndDrag[1] = ypos/height*2.-1.;


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
		glfwSetCursor(Visu->window,NULL);
	}




}


/*
void Visu_SaveToImageFile(Visu* Visu) {

}
 */




void Visu_main(Visu* Visu, Grid* Grid, Physics* Physics, Particles* Particles, Numerics* Numerics, BC* BCStokes, Char* Char, MatProps* MatProps, EqSystem* EqStokes, EqSystem* EqThermal, Numbering* NumStokes, Numbering* NumThermal)
{
	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                 VISUALIZATION                              //
	//                                                                            //
	//============================================================================//
	//============================================================================//
	GLfloat shiftIni[3];
	do  {

		glfwPollEvents();
		Visu_checkInput(Visu);
		if (Visu->update) {
			//printf("Updating the plot\n");
			glClearColor(1, 1, 1, 0.0); // black

			glEnable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
			glStencilMask(0xFF);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			glStencilMask(0x00);




			if (Visu->initPassivePart) {
				Particles_initPassive(Particles, Grid);
				Visu->initPassivePart = false;
			}


			shiftIni[0] = Visu->shift[0];
			shiftIni[1] = Visu->shift[1];
			shiftIni[2] = Visu->shift[2];

			Visu->shift[0] -= (Grid->xmax_ini-Grid->xmin_ini)*Visu->shiftFac[0]*Visu->scale;
			Visu->shift[1] += (Grid->ymax_ini-Grid->ymin_ini)*Visu->shiftFac[1]*Visu->scale;
			Visu->shift[2] +=                 1.0*Visu->shiftFac[2];

			// Update the grid
			if (Visu->updateGrid) {
				if (BCStokes->SetupType==PureShear || BCStokes->SetupType==Sandbox) {
					Visu_updateVertices(Visu, Grid);
					glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
					glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu->vertices, GL_STATIC_DRAW);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
				}
			}

			//============================================================================
			// 								PLOT PARTICLE
			if (Visu->showParticles) {
				glEnable(GL_STENCIL_TEST);
				// Draw the box in black and use it to set the stencil
				// Particles fragment will be drawn only where the stencil is 1 (i.e. inside the box)
				glStencilFunc(GL_ALWAYS, 1, 0xFF); // All fragments should update the stencil buffer
				glStencilMask(0xFF); // Enable writing to the stencil buffer
				glBindVertexArray(Visu->VAO);
				glUseProgram(Visu->ParticleBackgroundShaderProgram);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);
				Visu_updateUniforms(Visu);
				//Visu_update(Visu, window, Grid, Physics, BCStokes, Char);
				glDrawElements(GL_TRIANGLES, Visu->ntrivert, GL_UNSIGNED_INT, 0);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
				glUseProgram(0);
				glBindVertexArray(0);

				glStencilFunc(GL_EQUAL, 1, 0xFF);
				glStencilMask(0x00);  // disable writing to the buffer
				//glDisable(GL_DEPTH_TEST);



				glBindVertexArray(Visu->VAO_part);
				glUseProgram(Visu->ParticleShaderProgram);


				glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_part);
				//glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_partMesh);
				// update the buffer containing the particles
				Visu_particles(Visu, Particles, Grid);
				Visu_updateUniforms(Visu);
				glBufferSubData(GL_ARRAY_BUFFER, 0, 4*Particles->n*sizeof(GLfloat), Visu->particles);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				//glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_partMesh);
				//glBindBuffer(GL_ARRAY_BUFFER, 0);
				glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, Visu->particleMeshRes+2, Particles->n);
				//printf("Visu->particleMeshRes= %i\n",Visu->particleMeshRes);
				//glDrawArraysInstanced(GL_TRIANGLES, 0, 3, Particles.n);
				//
				glUseProgram(0);
				glBindVertexArray(0);
				glDisable(GL_STENCIL_TEST);
			}
			// 								PLOT PARTICLE
			//============================================================================



			Visu->shift[0] += 2*(Grid->xmax_ini-Grid->xmin_ini)*Visu->shiftFac[0]*Visu->scale;
			Visu->shift[1] -= 2*(Grid->ymax_ini-Grid->ymin_ini)*Visu->shiftFac[1]*Visu->scale;
			Visu->shift[2] -=                   2.0*Visu->shiftFac[2];

			//============================================================================
			// 								PLOT GRID DATA


			//glDisable(GL_DEPTH_TEST);


			// ****** Bind shader textures, arrays and buffers
			glBindVertexArray(Visu->VAO);
			glUseProgram(Visu->ShaderProgram);
			glBindTexture(GL_TEXTURE_2D, Visu->TEX);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);

			// 1. Update data
			Visu_update(Visu, Grid, Physics, BCStokes, Char, MatProps, EqStokes, EqThermal, NumStokes, NumThermal, Numerics);
			Visu_alphaValue(Visu, Grid, Physics);
			// update the content of Visu->U
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RG, Grid->nxEC, Grid->nyEC, 0, GL_RG, GL_FLOAT, Visu->U);	// load the updated Visu->U in the texture
			// 2. Draw
			glDrawElements(GL_TRIANGLES, Visu->ntrivert, GL_UNSIGNED_INT, 0);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			glBindTexture(GL_TEXTURE_2D, 0);
			glUseProgram(0);
			glBindVertexArray(0);


			//glEnable(GL_DEPTH_TEST);

			// ****** Unbind textures, arrays and buffers


			// 								PLOT GRID DATA
			//============================================================================





			//============================================================================
			// 								PLOT GLYPH
			if (Visu->showGlyphs) {


				glDisable(GL_DEPTH_TEST);


				glBindVertexArray(Visu->VAO_glyph);
				glUseProgram(Visu->GlyphShaderProgram);


				glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_glyph);

				// update the buffer containing the particles



				Visu_glyphs(Visu, Physics, Grid, Particles);
				Visu_updateUniforms(Visu);
				glBufferSubData(GL_ARRAY_BUFFER, 0, 4*Visu->nGlyphs*sizeof(GLfloat), Visu->glyphs);
				glBindBuffer(GL_ARRAY_BUFFER, 0);

				if (Visu->glyphMeshType==ThinArrow) {
					glDrawArraysInstanced(GL_LINE_STRIP, 0, Visu->nGlyphMeshVert, Visu->nGlyphs);
				} else {
					glDrawArraysInstanced(GL_TRIANGLES, 0, Visu->nGlyphMeshVert, Visu->nGlyphs);
				}

				glUseProgram(0);
				glBindVertexArray(0);

				glEnable(GL_DEPTH_TEST);
			}
			// 								PLOT GLYPH
			//============================================================================






			Visu->shift[0] = shiftIni[0];
			Visu->shift[1] = shiftIni[1];
			Visu->shift[2] = shiftIni[2];


			//============================================================================
			// 							  SAVE TO IMAGE FILE

			if (Visu->writeImages) {
				FILE *fptr;
				char fname[1024];
				char ftitle[1024];
				sprintf(fname,"%sFrame_%05i.png",Visu->outputFolder,Numerics->timeStep);
				sprintf(ftitle,"time_%5.5e.png",Physics->time);
				//sprintf(fname,"Frame_%04i.raw",timeStep);
				if ((fptr = fopen(fname,"w")) == NULL) {
					fprintf(stderr,"Failed to open the file for window dump\n");
					exit(0);
				}
				glPixelStorei(GL_PACK_ALIGNMENT,1);
				glReadBuffer(GL_BACK);
				glReadPixels(0,0,Visu->retinaScale*Visu->width,Visu->retinaScale*Visu->height,GL_RGBA,GL_UNSIGNED_BYTE,Visu->imageBuffer);
				//fwrite(Visu->imageBuffer,Visu->width*Visu->height*3,1,fptr);
				int result = writePNGImage(fname, Visu->retinaScale*Visu->width, Visu->retinaScale*Visu->height, Visu->imageBuffer, ftitle);
				if (result!=0) {
					printf("error: couldn't write png file\n");
					exit(0);
				}
				fclose(fptr);

			}


			// 							  SAVE TO IMAGE FILE
			//============================================================================

			glfwSwapBuffers(Visu->window);

		}
		Visu->update = false;
		Visu_checkInput(Visu);

		if (glfwWindowShouldClose(Visu->window))
			break;
		if (Numerics->timeStep==Numerics->nTimeSteps-1)
			Visu->paused = true;
	} while (Visu->paused);




	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          END OF VISUALIZATION                              //
	//                                                                            //
	//============================================================================//
	//============================================================================//
}


void Visu_residual(Visu* Visu, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering)
{
	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq, iEqStart, iEqEnd;
	int ixECStart, ixECEnd;
	int iyECStart, iyECEnd;
	int J,i;
	int ix, iy, I;
	int xLength, iGrid0;
	//EqSystem->normResidual = 0;


	if (Visu->type==TRes) {
		iEqStart = Numbering->subEqSystem0[0];
		iEqEnd   = Numbering->subEqSystem0[1];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[0];
	} else if (Visu->type==VxRes) {
		iEqStart = Numbering->subEqSystem0[0];
		iEqEnd   = Numbering->subEqSystem0[1];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxVx;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyVx;

		xLength 	= Grid->nxVx;
		iGrid0  	= Numbering->subEqSystem0Dir[0];
	} else if (Visu->type==VyRes) {
		iEqStart = Numbering->subEqSystem0[1];
		iEqEnd   = Numbering->subEqSystem0[2];
		ixECStart 	= 0;
		ixECEnd 	= Grid->nxVy;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyVy;

		xLength 	= Grid->nxVy;
		iGrid0  	= Numbering->subEqSystem0Dir[1];
	} else if (Visu->type==PRes || Visu->type==PfRes) {
		iEqStart = Numbering->subEqSystem0[2];
		iEqEnd   = Numbering->subEqSystem0[3];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[2];
	} else if (Visu->type==PcRes) {
		iEqStart = Numbering->subEqSystem0[3];
		iEqEnd   = Numbering->subEqSystem0[4];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[3];
	}


	// Could be optimized by not looping over everything (be careful to the lower trianuglar contributions though; that's why I didn't do it yet)
#pragma omp parallel for private(iEq, i, J) schedule(static,32)
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		Residual[iEq] = EqSystem->b[iEq];
		for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
			J = EqSystem->J[i];
			Residual[iEq] += - (EqSystem->V[i]*EqSystem->x[J]);
			/*
			if (UPPER_TRI) {
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}*/
		}
	}

	if (UPPER_TRI) {

//#pragma omp parallel for private(iEq, i, J) schedule(static,32)
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
				J = EqSystem->J[i];
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}
		}

	}




	int C = iEqStart;
	int iGrid;
	for (iy = iyECStart; iy < iyECEnd; ++iy) {
		for (ix = ixECStart; ix < ixECEnd; ++ix) {
			iGrid = ix+iy*xLength + iGrid0;
			I = 2*(ix+iy*Grid->nxEC);
			if (Numbering->map[iGrid]>=0) {
				Visu->U[I] = Residual[C]/EqSystem->norm_b;
				C++;
			} else {
				Visu->U[I] = 1.0;
			}

		}
	}










	free(Residual);
}




#endif
