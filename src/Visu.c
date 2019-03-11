/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stddef.h>
#if (VISU)


bool shiftMod;
void Visu_Memory_allocate( Visu* Visu, Grid* Grid )
{
	Visu->U             = (GLfloat*)  malloc(2*Grid->nxEC*Grid->nyEC    * sizeof( GLfloat ));
	//Visu->vertices      = (GLfloat*)  malloc(Grid->nxS*Grid->nyS*2  * sizeof( GLfloat ));
	Visu->elements      = (GLuint*)   malloc(Visu->ntrivert    		* sizeof( GLuint ));

	Visu->vertices      = (GLfloat*)  malloc(4 * 4 * sizeof( GLfloat )); // 4 corners only
	Visu->particles 	= (GLfloat*) malloc (Visu->nParticles*4*sizeof(GLfloat));
	Visu->particleMesh 	= (GLfloat*) malloc ((Visu->particleMeshRes+2) *3*sizeof(GLfloat));
	

	Visu->nGlyphs 		= (int) ceil((double)Grid->nxS/(double)Visu->glyphSamplingRateX)*ceil((double)Grid->nyS/(double)Visu->glyphSamplingRateY);
	if (Visu->glyphSamplingRateX<1 || Visu->glyphSamplingRateY<1 ) {
		printf("warning!! Visu->Visu->glyphSamplingRateX<1 or Visu->Visu->glyphSamplingRateY<1\n");
		Visu->glyphSamplingRateX = INT_MAX;
		Visu->glyphSamplingRateY = INT_MAX;
		Visu->nGlyphs 		= 1;
	}


	Visu->glyphs 		= (GLfloat*) malloc ( Visu->nGlyphs *4*sizeof(GLfloat));
	Visu->imageBuffer 	= (unsigned char*) malloc(Visu->retinaScale*Visu->retinaScale*4*Visu->width*Visu->height*sizeof(unsigned char)); // does not consider image resizing


	if (Visu->glyphMeshType==VisuGlyphMeshType_Triangle) {
		Visu->nGlyphMeshVert = 3;
	}
	else if (Visu->glyphMeshType==VisuGlyphMeshType_ThinArrow) {
		Visu->nGlyphMeshVert = 6;
	}
	else if (Visu->glyphMeshType==VisuGlyphMeshType_ThickArrow) {
		Visu->nGlyphMeshVert = 18;
	} else if (Visu->glyphMeshType==VisuGlyphMeshType_TensorCross) {
			Visu->nGlyphMeshVert = 18;
	} else {
		printf("error in Visu_Memory_allocate: unknwon glyphMeshType");
	}

	Visu->glyphMesh 	= (GLfloat*) malloc ( Visu->nGlyphMeshVert *2*sizeof(GLfloat));

}




void Visu_Memory_free( Visu* Visu )
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



	glfwSetInputMode(Visu->window, GLFW_STICKY_KEYS, 1);


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
//#pragma omp parallel for private(iNode, thisParticle) OMP_SCHEDULE
	FOR_PARTICLES
	Visu->particles[C] = thisParticle->x;
	Visu->particles[C+1] = thisParticle->y;

	if (Visu->typeParticles == VisuType_PartPhase) {
		Visu->particles[C+2] = thisParticle->phase;
	} else if (Visu->typeParticles == VisuType_PartTemp) {
#if (HEAT)
		Visu->particles[C+2] = thisParticle->T;
#else
		Visu->particles[C+2] = 0;
#endif
	} else if (Visu->typeParticles == VisuType_PartSigma_xx) {
		Visu->particles[C+2] = thisParticle->sigma_xx_0;
	} else if (Visu->typeParticles == VisuType_PartSigma_xy) {
		Visu->particles[C+2] = thisParticle->sigma_xy_0;
	} else if (Visu->typeParticles == VisuType_PartDeltaP) {
#if (DARCY)
		Visu->particles[C+2] = thisParticle->DeltaP0;
#else
		Visu->particles[C+2] = 0.0;
#endif
	} else if (Visu->typeParticles == VisuType_PartPorosity) {
#if (DARCY)
		Visu->particles[C+2] = thisParticle->phi;
#else
		Visu->particles[C+2] = 0.0;
#endif
	} else if (Visu->typeParticles == VisuType_PartStrain) {
#if (STORE_PLASTIC_STRAIN)
		Visu->particles[C+2] = thisParticle->strain;
#else
		Visu->particles[C+2] = 0.0;
#endif
	} else if (Visu->typeParticles == VisuType_PartExtraField) {
#if (EXTRA_PART_FIELD)
		Visu->particles[C+2] = thisParticle->extraField;
#else
		Visu->particles[C+2] = 0.0;
#endif
	}
	Visu->particles[C+3] = thisParticle->passive;

	C += 4;
	END_PARTICLES

}



void Visu_glyphs(Model* Model)
{
	Visu* Visu 				= &(Model->Visu);
	Physics* Physics 		= &(Model->Physics);
	Grid* Grid 				= &(Model->Grid);


	int ix, iy, iCell;
	int C = 0;

#if (DARCY)
	compute perm_eta_f, phi, dPfdx, dPfdy;
#endif



	int n = 0;
	if (Visu->glyphType == VisuGlyphType_StokesVelocity) {
		for (iy = 0; iy < Grid->nyS; iy+=Visu->glyphSamplingRateY) {
			for (ix = 0; ix < Grid->nxS; ix+=Visu->glyphSamplingRateX) {
				iCell = ix + iy*Grid->nxEC; // Cell at the left of the lowest Vx node

				if (Physics->phase[iCell]!=Physics->phaseAir && Physics->phase[iCell]!=Physics->phaseWater) {
					Visu->glyphs[C+0] = Grid->xmin + ix*Grid->dx;
					Visu->glyphs[C+1] = Grid->ymin + iy*Grid->dy;



					Visu->glyphs[C+2] = (Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2.0;
					Visu->glyphs[C+3] = (Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2.0;
					C+=4;
					n++;
				}


			}


		}
	}


	else if (Visu->glyphType == VisuGlyphType_DarcyGradient) {
#if (DARCY)
		for (iy = 0; iy < Grid->nyS; iy+=Visu->glyphSamplingRateY) {
			for (ix = 0; ix < Grid->nxS; ix+=Visu->glyphSamplingRateX) {
				iCell = ix + iy*Grid->nxEC; // Cell at the left of the lowest Vx node

				if (Physics->phase[iCell]!=Physics->phaseAir && Physics->phase[iCell]!=Physics->phaseWater) {
					Visu->glyphs[C+0] = Grid->xmin + ix*Grid->dx;
					Visu->glyphs[C+1] = Grid->ymin + iy*Grid->dy;



					perm_eta_f = Physics->perm_eta_f[iCell];
					phi = Physics->phi[iCell];
					perm_eta_f = Physics->perm_eta_f[iCell];
					dPfdx = (Physics->Pf[ix+1 + iy*Grid->nxEC] - Physics->Pf[ix-1 + iy*Grid->nxEC])/2.0/Grid->dx;
					dPfdy = (Physics->Pf[ix + (iy+1)*Grid->nxEC] - Physics->Pf[ix + (iy-1)*Grid->nxEC])/2.0/Grid->dy;

					Visu->glyphs[C+2] = perm_eta_f * (-dPfdx + Physics->rho_f*Physics->g[0]); // DarcyVelX
					Visu->glyphs[C+3] = perm_eta_f * (-dPfdy + Physics->rho_f*Physics->g[1]); // DarcyVelY

									C+=4;
					n++;
				}
			}
		}
#endif
		//printf("GradSouth = %.1e\n", (Physics->psi[ix   + (iy+1)*Grid->nxEC]-Physics->psi[ix+iy*Grid->nxEC]+dy)/dy );
	} else if (Visu->glyphType == VisuGlyphType_DeviatoricStressTensor) {
		compute Tau, psi, Sxy, SII; // Tau is some non dimensional stress and spi is the angle between sigma1 and x
		for (iy = 1; iy < Grid->nyEC-1; iy+=Visu->glyphSamplingRateY) {
			for (ix = 1; ix < Grid->nxEC-1; ix+=Visu->glyphSamplingRateX) {
				iCell = ix + iy*Grid->nxEC; // Cell at the left of the lowest Vx node

				if (Physics->phase[iCell]!=Physics->phaseAir && Physics->phase[iCell]!=Physics->phaseWater) {
					Visu->glyphs[C+0] = Grid->xmin-Grid->dx/2.0 + ix*Grid->dx;
					Visu->glyphs[C+1] = Grid->ymin-Grid->dy/2.0 + iy*Grid->dy;

					Sxy = Interp_ECVal_Cell2Node_Local(Physics->sigma_xy_0, ix, iy, Grid->nxEC);
					Tau = Physics->sigma_xx_0[iCell] / Sxy;
					SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);

					//if (Physics->sigma_xx_0[iCell]<0.0) { // need to check if it's the proper condition for the switch
					if (Sxy>0.0){
						psi = atan(-Tau+sqrt(Tau*Tau+1));
					} else {
						psi = atan(-Tau-sqrt(Tau*Tau+1));
					}

					//psi = atan(-Tau+sqrt(Tau*Tau+1));

					Visu->glyphs[C+2] = cos(psi) * SII;//(Physics->Vx[ix  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy+1)*Grid->nxVx])/2.0;
					Visu->glyphs[C+3] = sin(psi) * SII;//(Physics->Vy[ix  +(iy  )*Grid->nxVy] + Physics->Vy[ix+1+(iy  )*Grid->nxVy])/2.0;
					C+=4;
					n++;
				}


			}


		}
	}





	Visu->nGlyphs = n;


	/*
	int C = 0;
	INIT_PARTICLE
	FOR_PARTICLES
		Visu->particles[C] = thisParticle->x;
		Visu->particles[C+1] = thisParticle->y;

		if (Visu->typeParticles == VisuType_Phase) {
			Visu->particles[C+2] = thisParticle->phase;
		} else if (Visu->typeParticles == VisuType_PartTemp) {
			Visu->particles[C+2] = thisParticle->T;
		} else if (Visu->typeParticles == VisuType_PartSigma_xx) {
			Visu->particles[C+2] = thisParticle->sigma_xx_0;
		} else if (Visu->typeParticles == VisuType_PartSigma_xy) {
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

	if (Visu->glyphMeshType==VisuGlyphMeshType_Triangle) {

		Visu->glyphMesh[0] = 0.0;
		Visu->glyphMesh[1] = -width;
		Visu->glyphMesh[2] = 0.0;
		Visu->glyphMesh[3] = +width;
		Visu->glyphMesh[4] = size;
		Visu->glyphMesh[5] = 0.0;
	}
	else if (Visu->glyphMeshType==VisuGlyphMeshType_ThinArrow) {
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
	else if (Visu->glyphMeshType==VisuGlyphMeshType_ThickArrow) {

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



	} else if (Visu->glyphMeshType==VisuGlyphMeshType_TensorCross) {

		// Stick, triangle 1
		Visu->glyphMesh[ 0] = 0.0 - size/2.0;
		Visu->glyphMesh[ 1] = stickWidth;
		Visu->glyphMesh[ 2] = 0.0 - size/2.0;
		Visu->glyphMesh[ 3] = -stickWidth;
		Visu->glyphMesh[ 4] = size - size/2.0;
		Visu->glyphMesh[ 5] = -stickWidth;

		// Stick, triangle 2
		Visu->glyphMesh[ 6] = size - size/2.0;
		Visu->glyphMesh[ 7] = -stickWidth;
		Visu->glyphMesh[ 8] = size - size/2.0;
		Visu->glyphMesh[ 9] = +stickWidth;
		Visu->glyphMesh[10] = 0.0 - size/2.0;
		Visu->glyphMesh[11] = stickWidth;



	}
	else {
		printf("error: unknown Visu->glyphMeshType");
	}



}





void Visu_init(Visu* Visu, Grid* Grid, Particles* Particles, Char* Char, Input* Input)
{
	Visu->timeStep_residual = -1;
	Visu->stepsSinceLastRender = 0;
    Visu->timeSinceLastRender = 0.0;

	Visu->renderCounter = 0;

	if (Visu->renderTimeFrequency>0.0) {
		Visu->useTimeFrequency = true;
	} else {
		Visu->useTimeFrequency = false;
	}

	Visu->updateGrid = false;
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
	/*
	sprintf(Visu->VertexShaderFile 						,"%s../Shaders/Default/shader.vs"					, (Input->currentFolder));
	sprintf(Visu->FragmentShaderFile 					,"%s../Shaders/Default/shader.fs"					, (Input->currentFolder));
	sprintf(Visu->ParticleVertexShaderFile 				,"%s../Shaders/Default/particleShader.vs"			, (Input->currentFolder));
	sprintf(Visu->ParticleGeometryShaderFile 			,"%s../Shaders/Default/particleShader.gs"			, (Input->currentFolder));
	sprintf(Visu->ParticleFragmentShaderFile 			,"%s../Shaders/Default/particleShader.fs"			, (Input->currentFolder));
	sprintf(Visu->ParticleBackgroundVertexShaderFile 	,"%s../Shaders/Default/particleBackgroundShader.vs"	, (Input->currentFolder));
	sprintf(Visu->ParticleBackgroundFragmentShaderFile 	,"%s../Shaders/Default/particleBackgroundShader.fs"	, (Input->currentFolder));
	sprintf(Visu->GlyphVertexShaderFile 				,"%s../Shaders/Default/glyphShader.vs"				, (Input->currentFolder));
	sprintf(Visu->GlyphFragmentShaderFile 				,"%s../Shaders/Default/glyphShader.fs"				, (Input->currentFolder));
	*/

	sprintf(Visu->VertexShaderFile 						,"%s%s/shader.vs"					, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->FragmentShaderFile 					,"%s%s/shader.fs"					, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->ParticleVertexShaderFile 				,"%s%s/particleShader.vs"			, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->ParticleGeometryShaderFile 			,"%s%s/particleShader.gs"			, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->ParticleFragmentShaderFile 			,"%s%s/particleShader.fs"			, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->ParticleBackgroundVertexShaderFile 	,"%s%s/particleBackgroundShader.vs"	, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->ParticleBackgroundFragmentShaderFile 	,"%s%s/particleBackgroundShader.fs"	, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->GlyphVertexShaderFile 				,"%s%s/glyphShader.vs"				, Input->currentFolder,Visu->shaderFolder);
	sprintf(Visu->GlyphFragmentShaderFile 				,"%s%s/glyphShader.fs"				, Input->currentFolder,Visu->shaderFolder);

	Visu->ShaderProgram 					= 0;
	Visu->ParticleShaderProgram 			= 0;
	Visu->ParticleBackgroundShaderProgram 	= 0;
	Visu->GlyphShaderProgram 				= 0;

	// Generate reference to objects (indexes that act as pointers to graphic memory)
	// =======================================
	Visu->VAO = 0; // Reference to the Vertex   array object
	Visu->VBO = 0; // Reference to the Vertex   buffer object
	Visu->EBO = 0; // Reference to the Element  buffer object
	Visu->TEX = 0; // Reference to the Element  buffer object



	// And assigned them to objects (stored in the graphic memory)
	// =======================================
	glGenVertexArrays 	(1, &Visu->VAO);

	glGenBuffers 		(1, &Visu->VBO);
	glGenBuffers 		(1, &Visu->EBO);
	glGenTextures		(1, &Visu->TEX);

	glGenVertexArrays	(1, &Visu->VAO_part);
	glGenBuffers		(1, &Visu->VBO_part);
	glGenBuffers		(1, &Visu->VBO_partMesh);

	glGenVertexArrays	(1, &Visu->VAO_glyph);
	glGenBuffers		(1, &Visu->VBO_glyph);
	glGenBuffers		(1, &Visu->VBO_glyphMesh);



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

	if (Visu->filter == VisuFilterType_Linear) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	} else if (Visu->filter == VisuFilterType_Nearest) {
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
		Visu->scale = 2.0/(1.05*(Grid->xmax-Grid->xmin)*(1.0+2.0*Visu->shiftFac[0]));
	}
	else {
		Visu->scale = 2.0/(1.05*(Grid->ymax-Grid->ymin)*(1.0+2.0*Visu->shiftFac[1])*ratio);
	}

	// Visu->scale = 2.0/(0.85*(Grid->xmax-Grid->xmin)*(1+2*Visu->shiftFac[0]));


	GLint loc = glGetUniformLocation(Visu->ShaderProgram, "one_ov_log_of_10");
	glUniform1f(loc, 1.0/log(10));

	Visu->colorScale[0] = -0.5;
	Visu->colorScale[1] =  0.5;
	Visu->partColorScale[0] = -0.5;
	Visu->partColorScale[1] =  0.5;
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
	Visu->initPassivePart = false;
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

	if (mods == GLFW_MOD_SHIFT) {
		shiftMod = true;
	} else {
		shiftMod = false;
	}
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

	//compute Ratio = (Grid->xmax-Grid->xmin)/(Grid->ymax-Grid->ymin);

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
void Visu_ECVal_updateGlobal(Visu* Visu, Grid* Grid, compute* CellValue)
{
	// UC is a scalar CellValue defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	//int iNW, iNE, iSW, iSE;
	// CellValue interpolated on the center nodes
	// ======================================
#pragma omp parallel for private(iy, ix, I) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			I = 2* (ix + iy*Grid->nxEC);
			Visu->U[I] = CellValue[ix + iy*Grid->nxEC];//(CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}
}

void Visu_ECVal_updateGlobal_i(Visu* Visu, Grid* Grid, int* CellValue)
{
	// UC is a scalar CellValue defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	//int iNW, iNE, iSW, iSE;
	// CellValue interpolated on the center nodes
	// ======================================
#pragma omp parallel for private(iy, ix, I) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			I = 2* (ix + iy*Grid->nxEC);
			Visu->U[I] = CellValue[ix + iy*Grid->nxEC];//(CellValue[iNW] + CellValue[iNE] + CellValue[iSW] + CellValue[iSE])/4;
		}
	}
}








void Visu_strainRate(Model* Model)
{

	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);


	int iy, ix;
	int I = 0;
	//compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	compute EII;
	//Visu_ECVal_updateGlobal (Visu, Grid, Physics->sigma_xx_0, BC->SetupType);
#pragma omp parallel for private(iy, ix, I, EII) OMP_SCHEDULE
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = (ix+iy*Grid->nxEC);
			Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);
			Visu->U[2*I] = EII;
		}
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
}


void Visu_rotationRate(Visu* Visu, Grid* Grid, Physics* Physics)
{

	int iy, ix;
	int I = 0;
	//compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	//Visu_ECVal_updateGlobal (Visu, Grid, Physics->sigma_xx_0, BC->SetupType);
#pragma omp parallel for private(iy, ix, I) OMP_SCHEDULE
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = (ix+iy*Grid->nxEC);
			compute dVxdy, dVydx;

			int iNode, Ix, Iy;
			int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
			int IyMod[4] = {0,0,1,1};

			// Method A: using the averageing of derivatives on the four nodes
			// Compute Eps_xy at the four nodes of the cell
			// 1. Sum contributions
			dVxdy = 0;
			dVydx = 0;
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];

				dVxdy += .25*( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
									  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


				dVydx += .25*( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
									  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;


			}
			Visu->U[2*I] = 0.5*(dVydx-dVxdy);
		}
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
}


void Visu_velocity(Visu* Visu, Grid* Grid, Physics* Physics)
{
	int iy, ix;
	int I = 0;
	compute A, B;



#pragma omp parallel for private(iy, ix, I, A, B) OMP_SCHEDULE
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = 2*(ix+iy*Grid->nxEC);


			A  = (Physics->Vx[ix-1  +(iy  )*Grid->nxVx] + Physics->Vx[ix  +(iy  )*Grid->nxVx])/2.0;
			B  = (Physics->Vy[ix    +(iy-1)*Grid->nxVy] + Physics->Vy[ix  +(iy  )*Grid->nxVy])/2.0;
			Visu->U[I] = sqrt(A*A + B*B);

		}
	}
}


void Visu_divV(Visu* Visu, Grid* Grid, Physics* Physics) {
	int iy, ix;
	int I = 0;

	compute dx, dy, divV;

//#pragma omp parallel for private(iy, ix, I, dx, dy, divV) OMP_SCHEDULE
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = 2*(ix+iy*Grid->nxEC);

			dx = Grid->DXS[ix-1];
			dy = Grid->DYS[iy-1];
			divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/dx;
			divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/dy;

			Visu->U[I] = divV;

		}
	}

}


void Visu_stress(Model* Model)
{


	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);

	int iy, ix;
	int I = 0;
	//compute A, B;
	// Loop through Vx nodes
	//printf("=== Visu Vel ===\n");
	compute SII;
	//Visu_ECVal_updateGlobal (Visu, Grid, Physics->sigma_xx_0, BC->SetupType);
#pragma omp parallel for private(iy, ix, I, SII) OMP_SCHEDULE
	for (iy=1; iy<Grid->nyEC-1; iy++){
		for (ix=1; ix<Grid->nxEC-1; ix++) {
			I = (ix+iy*Grid->nxEC);
			SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);
			Visu->U[2*I] = SII;
			//Visu->U[2*I] = Physics->sigma_xx_0[I];
			/*
			//  Compute sigmaII0
			compute sq_sigma_xy0, sigma_xx0, SII0;
			sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
			sigma_xx0  = Physics->sigma_xx_0[I];// + Physics->Dsigma_xx_0[iCell];
			SII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

			
			Visu->U[2*I] = (SII-SII0)/SII*30.0;
			*/
		}
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
}


void Visu_SIIOvYield(Model* Model) 
{

	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Physics* Physics 		= &(Model->Physics);

	Visu_stress(Model);
	compute sigmaII;


	compute sigma_y, Pe;

	int iCell;
#if (DARCY)
	compute phi = 0.0;
	compute phiCrit = Numerics->phiCrit;
#endif
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


			sigma_y = Physics->Tau_y[iCell];

			sigmaII = Physics_StressInvariant_getLocalCell(Model, ix, iy);
			
			Visu->U[2*iCell] = sigmaII/sigma_y;
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

void Visu_POvPlitho(Model* Model) 
{

	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);

	int iy, ix, iCell, iCellS, iCellN, iCellW, iCellE;
	compute rho_g_h;
#if (DARCY)
	compute rho_f_g_h;
#endif
	//int ixStart, ixEnd, ixInc;
	//int iyStart, iyEnd, iyInc;

	compute* Plitho = (compute*) malloc(Grid->nECTot * sizeof(compute));
#if (DARCY)
	compute* Phydro = (compute*) malloc(Grid->nECTot * sizeof(compute));
#endif
	//printf("enter Plitho\n");

	// Contribution of gy
	if (Physics->g[1]>0){
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				iCell = ix + iy*Grid->nxEC;
				iCellS = ix + (iy-1)*Grid->nxEC;
				if (iy==0) {
					rho_g_h = Physics->rho[iCell] * Physics->g[1] * (-0.5*Grid->DYEC[iy] );
#if (DARCY)
					rho_f_g_h = Physics->rho_f * Physics->g[1] * (-0.5*Grid->DYEC[iy] );
#endif
				} else {
					rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellS]) * Physics->g[1] * Grid->DYEC[iy-1] ;
#if (DARCY)
					rho_f_g_h += Physics->rho_f * Physics->g[1] * (Grid->DYEC[iy-1] );
#endif
				}
				Plitho[iCell] = rho_g_h;
#if (DARCY)
				Phydro[iCell] = rho_f_g_h;
#endif
			}
		}

	} else {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			for (iy = Grid->nyEC-1; iy >= 0; --iy) {

				iCell = ix + iy*Grid->nxEC;
				iCellN = ix + (iy+1)*Grid->nxEC;
				iCellS = ix + (iy-1)*Grid->nxEC;
				if (iy==Grid->nyEC-1) {
					rho_g_h = Physics->rho[iCell] * -Physics->g[1] * (-0.5*Grid->DYEC[iy-1] );
#if (DARCY)
					rho_f_g_h = Physics->rho_f * -Physics->g[1] * (-0.5*Grid->DYEC[iy-1] );
#endif
				} else {
					rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellN]) * -Physics->g[1] * Grid->DYEC[iy] ;
#if (DARCY)
					rho_f_g_h += Physics->rho_f * -Physics->g[1] * (Grid->DYEC[iy] );
#endif
				}
				//printf("ix = %i, iy = %i, rhogh = %.2e, Physics->rho[iCell] = %.2e\n", ix, iy, rho_g_h,Physics->rho[iCell]);

				Plitho[iCell] = rho_g_h;
#if (DARCY)
				Phydro[iCell] = rho_f_g_h;
#endif
			}
		}
	}


	if (abs(Physics->g[0])>1E-8) {
		// Contribution of gx
		if (Physics->g[0]>0){
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					iCellW = ix-1 + (iy)*Grid->nxEC;
					if (ix==0) {
						rho_g_h = Physics->rho[iCell] * Physics->g[0] * (-0.5*Grid->DXEC[ix] );
#if (DARCY)
						rho_f_g_h = Physics->rho_f * Physics->g[0] * (-0.5*Grid->DXEC[ix] );
#endif
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellW]) * Physics->g[0] * Grid->DXEC[ix-1] ;
#if (DARCY)
						rho_f_g_h += Physics->rho_f * Physics->g[0] * (Grid->DXEC[ix-1] );
#endif
					}
					Plitho[iCell] += rho_g_h;
#if (DARCY)
					Phydro[iCell] += rho_f_g_h;
#endif
				}
			}
		} else {

			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = Grid->nxEC-1; ix >= 0; --ix) {
					iCell = ix + iy*Grid->nxEC;
					iCellE = ix+1 + (iy)*Grid->nxEC;
					iCellW = ix-1 + (iy)*Grid->nxEC;
					if (ix==Grid->nxEC-1) {
						rho_g_h = Physics->rho[iCell] * -Physics->g[0] * (-0.5*Grid->DXEC[ix-1] );
#if (DARCY)
						rho_f_g_h = Physics->rho_f * -Physics->g[0] * (-0.5*Grid->DXEC[ix-1] );
#endif
					} else {
						rho_g_h += 0.5*(Physics->rho[iCell]+Physics->rho[iCellE]) * -Physics->g[0] * Grid->DXEC[ix] ;
#if (DARCY)
						rho_f_g_h += Physics->rho_f * -Physics->g[0] * (Grid->DXEC[ix] );
#endif
					}
					Plitho[iCell] += rho_g_h;
#if (DARCY)
					Phydro[iCell] += rho_f_g_h;
#endif
				}
			}
		}
	}

	compute SII;
	for (iy=0; iy<Grid->nyEC; ++iy) {
		for (ix=0; ix<Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			// P Ov Plitho
			//Visu->U[2*iCell] = Physics->P[iCell]/Plitho[iCell];


			// For frictionAngle = 30 deg, Sigma_n = (Sigma3+P)/2.0
			SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);
			//Sigma3 = (-SII+Physics->P[iCell]);
#if (DARCY)
			if (Physics->phi[iCell]>Numerics->phiCrit) {
				Sigma_n = (-SII/2.0+Physics->Pc[iCell]); // actually -SII*sin(phi) + P
			} else {
				Sigma_n = (-SII/2.0+Physics->P[iCell]);  // actually -SII*sin(phi) + P
			}
#else
			//Sigma_n = (-SII/2.0+Physics->P[iCell]); // actually -SII*sin(phi) + P
#endif

			//Sigma_v = (-Physics->sigma_xx_0[iCell]+Physics->P[iCell]);
			//Visu->U[2*iCell] = Sigma_n/Plitho[iCell];
			Visu->U[2*iCell] = Physics->P[iCell]/Plitho[iCell];
			//Visu->U[2*iCell] = Sigma3/Plitho[iCell];
			//Visu->U[2*iCell] = Sigma_v/Plitho[iCell];


			//Lambda = (Physics->Pf[iCell]-Phydro[iCell])/(Physics->P[iCell]-Phydro[iCell]);
			//Lambda = (Physics->Pf[iCell]-Phydro[iCell])/(Plitho[iCell]-Phydro[iCell]);
			//Lambda = ((Plitho[iCell]-Physics->Pc[iCell])-Phydro[iCell])/(Plitho[iCell]-Phydro[iCell]);
			//Visu->U[2*iCell] = Lambda;
			//if (ix == 50) {
			//	printf("ix = %i, iy = %i, Lambda = %.2e, Phydro = %.2e, Plitho = %.2e\n", ix, iy, Lambda, Phydro[iCell],Plitho[iCell]);
			//}
		}
	}



	free(Plitho);
#if (DARCY)
	free(Phydro);
#endif
}

void Visu_PeOvYield(Model* Model) 
{
	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);


	Visu_stress(Model);
	compute sigmaII;



	compute Py, Pe;
#if (DARCY)
	compute phi = 0.0;
	compute phiCrit = Numerics->phiCrit;
#endif
	compute sigmaT;

	compute R = 2.0;

	int iCell;

	int ix, iy;
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
			//Pe 		= Physics->Pc[iCell];
#else
			Pe = Physics->P[iCell];
#endif


			sigmaT = Physics->cohesion[iCell]/R;

			sigmaII = Physics_StressInvariant_getLocalCell(Model, ix, iy);

			Py = sigmaII - sigmaT;





			Visu->U[2*iCell] = (Pe-Py)/Py;



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
/*
	float alpha;

	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, alpha) OMP_SCHEDULE
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
		if ( Physics->phase[i] == Physics->phaseAir || Physics->phase[i] == Physics->phaseWater ) {
			Visu->U[2*i+1] = 0.0;
		}
	}
	

	/*
	int type = 2;
	compute lowerThreshold = .1*Visu->colorScale[1];
	//compute upperThreshold = 1.0*Visu->colorScale[1];

		//if (Visu->alphaAbsThreshold>0.0) {
			for (i = 0; i < Grid->nECTot; ++i) {
				if (type == 0) {
					Visu->U[2*i+1] = 1.0;
				} else if (type == 1) {
					Visu->U[2*i+1] = (    (Visu->U[2*i]) + Visu->valueShift)/Visu->colorScale[1];
				} else if (type == 2) {
					Visu->U[2*i+1] = (fabs(Visu->U[2*i]) + Visu->valueShift - lowerThreshold)/Visu->colorScale[1];
				} else {
					printf("unknwon visu alpha type\n");
					exit(0);
				}


				Visu->U[2*i+1] = fmax(Visu->U[2*i+1],0.0);
				Visu->U[2*i+1] = fmin(Visu->U[2*i+1],1.0);



				if ( Physics->phase[i] == Physics->phaseAir || Physics->phase[i] == Physics->phaseWater ) {
					Visu->U[2*i+1] = 0.0;
				}


			}

	//}
	 */
	/*
	for (i = 0; i < Grid->nECTot; ++i) {
		Visu->U[2*i+1] = 0.0;
	}
	*/




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
	if (Visu->typeParticles == VisuType_PartPhase ) {
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







void Visu_update(Model* Model)
{
	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);


		Visu->valueScale 	=  Visu->colorMap[Visu->type].scale;
		Visu->valueShift 	= -Visu->colorMap[Visu->type].center;
		Visu->colorScale[0] =  0.0; // dummy
		Visu->colorScale[1] =  Visu->colorMap[Visu->type].max-Visu->colorMap[Visu->type].center;
		Visu->log10_on 		=  Visu->colorMap[Visu->type].log10on;
		Visu->alphaAbsThreshold = Visu->colorMap[Visu->type].alphaAbsThreshold;
		//printf("Visu->type = %d, Visu->valueScale = %.2e, backSR = %.2e, Visu->log10_on = %i, Visu->valueShift = %.2e\n",Visu->type, Visu->valueScale,Physics->epsRef,Visu->log10_on, Visu->valueShift);

		int i;
	char title[1024];
	switch (Visu->type) {
	case VisuType_Viscosity:
		glfwSetWindowTitle(Visu->window, "Viscosity");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->eta);

		break;
	case VisuType_Khi:
		glfwSetWindowTitle(Visu->window, "Khi");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->khi);
		break;
	case VisuType_Khib:
#if (DARCY)
		glfwSetWindowTitle(Visu->window, "Khi_b");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->khi_b);
		break;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
		break;
#endif

	case VisuType_StrainRate:
		glfwSetWindowTitle(Visu->window, "StrainRate");
		Visu_strainRate(Model);
		break;
	case VisuType_Stress:
		glfwSetWindowTitle(Visu->window, "Stress");
		Visu_stress(Model);
		break;
	case VisuType_Velocity:
		//glfwSetWindowTitle(Visu->window, "Velocity");
		Visu_velocity(Visu, Grid, Physics);

		sprintf(title,"Velocity, scale = %.2e",Visu->valueScale);
		glfwSetWindowTitle(Visu->window, title);
		/*
		Visu->valueScale 	= 0.2*sqrt(Physics->maxVx*Physics->maxVx+Physics->maxVy*Physics->maxVy);//(Physics->epsRef*Grid->xmax);
		Visu->valueShift 	= 0;
		Visu->colorScale[0] = -1;
		Visu->colorScale[1] =  1;
		Visu->log10_on 		= false	;
		*/
		break;
	case VisuType_VelocityDiv:
		glfwSetWindowTitle(Visu->window, "Velocity divergence, /!\\ values are computed using the updated dx, dy (i.e. values appear much larger)");
		Visu_divV(Visu, Grid, Physics);
		break;
	case VisuType_SIIOvYield:
		glfwSetWindowTitle(Visu->window, "Stress_II/Stress_y");
		Visu_SIIOvYield(Model);
		break;

	case VisuType_PeOvYield:
#if (DARCY)
		glfwSetWindowTitle(Visu->window, "Pe/Py");
		Visu_PeOvYield(Model);
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;

	case VisuType_Pressure:
		glfwSetWindowTitle(Visu->window, "Pressure");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->P);
		break;
	case VisuType_Density:
		//glfwSetWindowTitle(Visu->window, "Density*g, MatProps->rho0_g[0] = %.2e", MatProps->rho0_g[0]);
		sprintf(title,"Density");
		glfwSetWindowTitle(Visu->window, title);
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->rho);
		break;
	case VisuType_Temperature:
#if (HEAT)
		glfwSetWindowTitle(Visu->window, "Temperature");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->T); // Not optimal but good enough for the moment
#else
		glfwSetWindowTitle(Visu->window, "Temperature is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif

		break;
	case VisuType_FluidPressure:
			glfwSetWindowTitle(Visu->window, "Fluid pressure");
#if (DARCY)
			Visu_ECVal_updateGlobal(Visu, Grid, Physics->Pf); // Not optimal but good enough for the moment
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;
	case VisuType_CompactionPressure:
		glfwSetWindowTitle(Visu->window, "Compaction pressure");
#if (DARCY)
			Visu_ECVal_updateGlobal(Visu, Grid, Physics->Pc); // Not optimal but good enough for the moment
			//Visu->valueScale = 0.2;
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
			break;
	case VisuType_Permeability:
		glfwSetWindowTitle(Visu->window, "Permeability/eta_f");
#if (DARCY)
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->perm_eta_f); // Not optimal but good enough for the moment
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;




	case VisuType_Porosity:
		glfwSetWindowTitle(Visu->window, "Porosity");
#if (DARCY)
			Visu_ECVal_updateGlobal(Visu, Grid, Physics->phi); // Not optimal but good enough for the moment
#else
		glfwSetWindowTitle(Visu->window, "Darcy is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;
	case VisuType_Phase:
		glfwSetWindowTitle(Visu->window, "Phase");
		Visu_ECVal_updateGlobal_i(Visu, Grid, Physics->phase);
		break;


	case VisuType_VxRes:
	case VisuType_VyRes:
	case VisuType_PRes:
	case VisuType_PfRes:
	case VisuType_PcRes:

		if 		 (Visu->type==VisuType_VxRes) {
			glfwSetWindowTitle(Visu->window, "Vx residual");
		} else if(Visu->type==VisuType_VyRes) {
			glfwSetWindowTitle(Visu->window, "Vy residual");
		} else if(Visu->type==VisuType_PRes) {
			glfwSetWindowTitle(Visu->window, "P residual");
		} else if(Visu->type==VisuType_PfRes) {
			glfwSetWindowTitle(Visu->window, "Pf residual");
		} else if(Visu->type==VisuType_PcRes) {
			glfwSetWindowTitle(Visu->window, "Pc Residual");
		}

		Visu_residual(Model);

		break;

	case VisuType_TRes:
		glfwSetWindowTitle(Visu->window, "T residual");
		Visu_residual(Model);

		break;
	case VisuType_Strain:
#if (STORE_PLASTIC_STRAIN)
		glfwSetWindowTitle(Visu->window, "Strain");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->strain);
#else
		glfwSetWindowTitle(Visu->window, "Strain softening is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;
	case VisuType_Vorticity:
		glfwSetWindowTitle(Visu->window, "Vorticity");
		Visu_rotationRate(Visu, Grid, Physics);
		break;
	case VisuType_POvPlitho:
		glfwSetWindowTitle(Visu->window, "POvPlitho");
		Visu_POvPlitho(Model);
		break;
	case VisuType_Blank:
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
	case VisuType_EffectiveViscosity:
		glfwSetWindowTitle(Visu->window, "Effective Viscosity");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->Z);
		break;
	case VisuType_ShearModulus:
		glfwSetWindowTitle(Visu->window, "Shear Modulus");
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->G);
		break;
	case VisuType_ExtraField:
		glfwSetWindowTitle(Visu->window, "extraField");
#if (EXTRA_PART_FIELD)
		Visu_ECVal_updateGlobal(Visu, Grid, Physics->extraField); // Not optimal but good enough for the moment
		
#else
		glfwSetWindowTitle(Visu->window, "EXTRA_PART_FIELD is switched off");
		for (i=0;i<Grid->nECTot;i++) {
			Visu->U[2*i] = 0;
		}
#endif
		break;

	default:
		printf("Error: unknown Visu->type: %i",Visu->type);
	}

	switch (Visu->typeParticles) {
	case VisuType_PartPhase:
		Visu->partColorScale[0] = -3;
		Visu->partColorScale[1] =  3;
		break;
	case VisuType_PartTemp:
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_Temperature].max-Visu->colorMap[VisuType_Temperature].center)*Visu->colorMap[VisuType_Temperature].scale;
		break;
	case VisuType_PartSigma_xx:
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_Stress].max-Visu->colorMap[VisuType_Stress].center)*Visu->colorMap[VisuType_Stress].scale;
		break;
	case VisuType_PartSigma_xy:
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_Stress].max-Visu->colorMap[VisuType_Stress].center)*Visu->colorMap[VisuType_Stress].scale;
		break;
	case VisuType_PartDeltaP:
#if (DARCY)
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_CompactionPressure].max-Visu->colorMap[VisuType_CompactionPressure].center)*Visu->colorMap[VisuType_CompactionPressure].scale;
#endif
		break;
	case VisuType_PartPorosity:
#if (DARCY)
		Visu->partColorScale[0] = 0;
		Visu->partColorScale[1] = 1.0;
#endif
		break;
	case VisuType_PartStrain:
#if (STORE_PLASTIC_STRAIN)
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_Strain].max-Visu->colorMap[VisuType_Strain].center)*Visu->colorMap[VisuType_Strain].scale;
#endif
		break;
	case VisuType_PartExtraField:
#if (STORE_PLASTIC_STRAIN)
		Visu->partColorScale[0] =  0.0; // dummy
		Visu->partColorScale[1] =  (Visu->colorMap[VisuType_ExtraField].max-Visu->colorMap[VisuType_ExtraField].center)*Visu->colorMap[VisuType_ExtraField].scale;
#endif
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
		if (!shiftMod) {
			Visu->type =  VisuType_Viscosity;
			Visu->update = true;
		} else {
			Visu->type =  VisuType_EffectiveViscosity;
			Visu->update = true;
		}
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_2) == GLFW_PRESS) {
		Visu->type =  VisuType_StrainRate;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_3) == GLFW_PRESS) {
		Visu->type =  VisuType_Stress;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_4) == GLFW_PRESS) {
		Visu->type =  VisuType_Pressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_5) == GLFW_PRESS) {
		Visu->type =  VisuType_Density;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_6) == GLFW_PRESS) {
		Visu->type =  VisuType_Temperature;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_7) == GLFW_PRESS) {
		Visu->type =  VisuType_Velocity;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_8) == GLFW_PRESS) {
		Visu->type =  VisuType_FluidPressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_9) == GLFW_PRESS) {
		Visu->type =  VisuType_Permeability;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_0) == GLFW_PRESS) {
		Visu->type =  VisuType_Porosity;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_C) == GLFW_PRESS) {
		Visu->type =  VisuType_CompactionPressure;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_V) == GLFW_PRESS) {
		Visu->type =  VisuType_Phase;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_B) == GLFW_PRESS) {
		if (!shiftMod) {
			Visu->type =  VisuType_Strain;
			Visu->update = true;
		} else {
			Visu->typeParticles =  VisuType_PartStrain;
			Visu->update = true;
		}
		
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_N) == GLFW_PRESS) {
		Visu->type =  VisuType_Vorticity;
		Visu->update = true;
	}
	
	
	
	else if (glfwGetKey(Visu->window, GLFW_KEY_U) == GLFW_PRESS) {
		if (!shiftMod) {
			Visu->type =  VisuType_SIIOvYield;
			Visu->update = true;
		} else {
			Visu->type =  VisuType_POvPlitho;
			Visu->update = true;
		}
	}
	
	else if (glfwGetKey(Visu->window, GLFW_KEY_I) == GLFW_PRESS) {
		Visu->type =  VisuType_VelocityDiv;
		Visu->update = true;
	}

	//else if (glfwGetKey(Visu->window, GLFW_KEY_I) == GLFW_PRESS) {
		//Visu->type =  VisuType_PeOvYield;
		//Visu->update = true;
	//}
	/*
	else if (glfwGetKey(Visu->window, GLFW_KEY_I) == GLFW_PRESS) {
		if (!shiftMod) {
			Visu->type =  VisuType_ExtraField;
			Visu->update = true;
		} else {
			Visu->typeParticles =  VisuType_PartExtraField;
			
			Visu->update = true;
		}
	}
	*/
	else if (glfwGetKey(Visu->window, GLFW_KEY_A) == GLFW_PRESS) {
		if (!shiftMod) {
			Visu->type =  VisuType_Khi;
			Visu->update = true;
		} else {
			Visu->type =  VisuType_ShearModulus;
			Visu->update = true;
		}
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_S) == GLFW_PRESS) {
		Visu->type =  VisuType_Khib;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_M) == GLFW_PRESS) {
		Visu->type =  VisuType_Blank;
		Visu->update = true;
	}



	// Residuals
#if (HEAT)
	else if (glfwGetKey(Visu->window, GLFW_KEY_D) == GLFW_PRESS) {
		Visu->type =  VisuType_TRes;
		Visu->update = true;
	}
#endif
	else if (glfwGetKey(Visu->window, GLFW_KEY_F) == GLFW_PRESS) {
		Visu->type =  VisuType_VxRes;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_G) == GLFW_PRESS) {
		Visu->type =  VisuType_VyRes;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_H) == GLFW_PRESS) {
#if (DARCY)
		Visu->type =  VisuType_PfRes;
#else
		Visu->type =  VisuType_PRes;
#endif
		Visu->update = true;
	}
#if (DARCY)
	else if (glfwGetKey(Visu->window, GLFW_KEY_J) == GLFW_PRESS) {
		Visu->type =  VisuType_PcRes;
		Visu->update = true;
	}
#endif


	// Particles
	else if (glfwGetKey(Visu->window, GLFW_KEY_Q) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartPhase;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_W) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartTemp;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_E) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartSigma_xx;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_R) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartSigma_xy;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_T) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartDeltaP;
		Visu->update = true;
	}
	else if (glfwGetKey(Visu->window, GLFW_KEY_Y) == GLFW_PRESS) {
		Visu->typeParticles = VisuType_PartPorosity;
		Visu->update = true;
	}
	




	else if (glfwGetKey(Visu->window, GLFW_KEY_P) == GLFW_PRESS) {
		if (Visu->paused == false) {
			Visu->update = true;
			Visu->paused = true;
		} else {
			Visu->update = false;
			Visu->paused = true;
		}
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




void Visu_main(Model* Model)
{
	Visu* Visu 				= &(Model->Visu);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Particles* Particles 	= &(Model->Particles);
	Numerics* Numerics 		= &(Model->Numerics);

	

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
				Particles_initPassive(Particles, Grid, Physics);
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
				//if (BCStokes->SetupType==Stokes_PureShear || BCStokes->SetupType==Stokes_Sandbox) {
					Visu_updateVertices(Visu, Grid);
					glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO);
					glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), Visu->vertices, GL_STATIC_DRAW);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
				//}
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

#if (MULTI_VISU)
			int nSubOutput;
#if (DARCY)
			nSubOutput = 12;
#else
			nSubOutput = 7;
#endif
#if (HEAT)
			nSubOutput = 13;
#endif
			int iSubOutput;
			char typeName[1024];
			for (iSubOutput = 0; iSubOutput < nSubOutput; ++iSubOutput) {

				if (iSubOutput == 0) {
					Visu->type =  VisuType_StrainRate;
					//typeName = "StrainRate";
					strcpy(typeName, "StrainRate");
				} else if (iSubOutput == 1) {
					Visu->type =  VisuType_Pressure;
					strcpy(typeName, "Pressure");
				} else if (iSubOutput == 2) {
					Visu->type =  VisuType_Velocity;
					strcpy(typeName, "Velocity");
				} else if (iSubOutput == 3) {
					Visu->type =  VisuType_Stress;
					strcpy(typeName, "Stress");
				} else if (iSubOutput == 4) {
					Visu->type =  VisuType_Vorticity;
					strcpy(typeName,  "Vorticity");
				} else if (iSubOutput == 5) {
					Visu->type =  VisuType_Khi;
					strcpy(typeName, "Khi");
				} else if (iSubOutput == 6) {
					Visu->type =  VisuType_POvPlitho;
					strcpy(typeName, "POvPlitho");
				} else if (iSubOutput == 7) {
					Visu->type =  VisuType_Viscosity;
					strcpy(typeName, "Viscosity");
				} else if (iSubOutput == 8) {
					Visu->type =  VisuType_Porosity;
					strcpy(typeName, "Porosity");
				} else if (iSubOutput == 9) {
					Visu->type =  VisuType_CompactionPressure;
					strcpy(typeName, "CompactionPressure");
				} else if (iSubOutput == 10) {
					Visu->type =  VisuType_FluidPressure;
					strcpy(typeName, "FluidPressure");
				} else if (iSubOutput == 11) {
					Visu->type =  VisuType_Khib;
					strcpy(typeName,  "Khib");
				} else if (iSubOutput == 12) {
					Visu->type =  VisuType_Temperature;
					strcpy(typeName, "Temperature");
					//Visu->type =  VisuType_Permeability;
					//strcpy(typeName, "Permeability");
				}
				glDisable(GL_DEPTH_TEST);

				char fname[2048];
				sprintf(fname,"%s%s",Visu->outputFolder,typeName);
				struct stat st = {0};

				if (stat(fname, &st) == -1) {
					mkdir(fname, 0700);
				}

#endif
			//============================================================================
			// 								PLOT GRID DATA


			//glDisable(GL_DEPTH_TEST);


			// ****** Bind shader textures, arrays and buffers
			glBindVertexArray(Visu->VAO);
			glUseProgram(Visu->ShaderProgram);
			glBindTexture(GL_TEXTURE_2D, Visu->TEX);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Visu->EBO);

			// 1. Update data
			Visu_update(Model);
			Visu_alphaValue(Visu, Grid, Physics);
			// update the content of Visu->U
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RG, Grid->nxEC, Grid->nyEC, 0, GL_RG, GL_FLOAT, Visu->U);	// load the updated Visu->U in the texture
			// 2. Draw
			glDrawElements(GL_TRIANGLES, Visu->ntrivert, GL_UNSIGNED_INT, 0);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			glBindTexture(GL_TEXTURE_2D, 0);
			glUseProgram(0);
			glBindVertexArray(0);


			// ****** Unbind textures, arrays and buffers


			// 								PLOT GRID DATA
			//============================================================================



			//============================================================================
			// 								PLOT GLYPH
			if (Visu->showGlyphs) {
				//Visu->shift[0] -= 2*(Grid->xmax_ini-Grid->xmin_ini)*Visu->shiftFac[0]*Visu->scale; // to put the glyphs on the particles
				//Visu->shift[1] += 2*(Grid->ymax_ini-Grid->ymin_ini)*Visu->shiftFac[1]*Visu->scale;
				//Visu->shift[2] +=                   2.0*Visu->shiftFac[2];


				glDisable(GL_DEPTH_TEST);


				glBindVertexArray(Visu->VAO_glyph);
				glUseProgram(Visu->GlyphShaderProgram);


				glBindBuffer(GL_ARRAY_BUFFER, Visu->VBO_glyph);

				// update the buffer containing the particles



				Visu_glyphs(Model);
				Visu_updateUniforms(Visu);
				glBufferSubData(GL_ARRAY_BUFFER, 0, 4*Visu->nGlyphs*sizeof(GLfloat), Visu->glyphs);
				glBindBuffer(GL_ARRAY_BUFFER, 0);

				if (Visu->glyphMeshType==VisuGlyphMeshType_ThinArrow) {
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



			//============================================================================
			// 							  SAVE TO IMAGE FILE

			if (Visu->writeImages) {
				FILE *fptr;
				char fname[BUFFER_STRING_LENGTH];
				char ftitle[1024];
#if (MULTI_VISU)
				sprintf(fname,"%s%s/Frame_%05i.png",Visu->outputFolder,typeName,Visu->renderCounter);
#else
				sprintf(fname,"%s/Frame_%05i.png",Visu->outputFolder,Visu->renderCounter);
#endif
				sprintf(ftitle,"time_%5.5e.png",Physics->time);
				//sprintf(fname,"Frame_%04i.raw",timeStep);
				if ((fptr = fopen(fname,"w")) == NULL) {
					fprintf(stderr,"Failed to open the file for window dump\n");
					printf("%s/Frame_%05i.png",Visu->outputFolder,Numerics->timeStep);
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
#if (MULTI_VISU)
			}
#endif

			Visu->shift[0] = shiftIni[0];
			Visu->shift[1] = shiftIni[1];
			Visu->shift[2] = shiftIni[2];

			glfwSwapBuffers(Visu->window);

		}
		Visu->update = false;


		if (Visu->closeAtTheEndOfSimulation==false) {
			if (Numerics->timeStep==Numerics->nTimeSteps-1 && Visu->nonLinItisOver)
				Visu->paused = true;

			if (Physics->time+Physics->dtAdv >= Numerics->maxTime) {
				Visu->paused = true;
			}
		}

		Visu_checkInput(Visu);
		if (glfwWindowShouldClose(Visu->window))
			break;
	} while (Visu->paused);




	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                          END OF VISUALIZATION                              //
	//                                                                            //
	//============================================================================//
	//============================================================================//
}


void Visu_residual(Model* Model)
{

	Visu* Visu = &(Model->Visu);
	Grid* Grid = &(Model->Grid);
	Numerics* Numerics 		= &(Model->Numerics);
	EqSystem* EqSystem = NULL;
	Numbering* Numbering = NULL;
	if (Visu->type==VisuType_TRes) { 
		EqSystem = &(Model->EqThermal);
		Numbering  = &(Model->NumThermal);
	} else if (Visu->type==VisuType_VxRes || Visu->type==VisuType_VyRes || Visu->type==VisuType_PRes || Visu->type==VisuType_PfRes || Visu->type==VisuType_PcRes) {
		EqSystem = &(Model->EqStokes);
		Numbering  = &(Model->NumStokes);
	}

	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq, iEqStart;//, iEqEnd;
	int ixECStart, ixECEnd;
	int iyECStart, iyECEnd;
	int J,i;
	int ix, iy, I;
	int xLength, iGrid0;
	//EqSystem->normResidual = 0;
	if (Visu->timeStep_residual != Numerics->timeStep) {
		for (i=0; i<EqSystem->nEq; ++i) {
			EqSystem->x[i] /= EqSystem->S[i];
			EqSystem->b[i] *= EqSystem->S[i];
		}
		Visu->timeStep_residual = Numerics->timeStep;
	}

	if (Visu->type==VisuType_TRes) {
		iEqStart = Numbering->subEqSystem0[0];
		//iEqEnd   = Numbering->subEqSystem0[1];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[0];
	} else if (Visu->type==VisuType_VxRes) {
		iEqStart = Numbering->subEqSystem0[0];
		//iEqEnd   = Numbering->subEqSystem0[1];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxVx;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyVx;

		xLength 	= Grid->nxVx;
		iGrid0  	= Numbering->subEqSystem0Dir[0];
	} else if (Visu->type==VisuType_VyRes) {
		iEqStart = Numbering->subEqSystem0[1];
		//iEqEnd   = Numbering->subEqSystem0[2];
		ixECStart 	= 0;
		ixECEnd 	= Grid->nxVy;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyVy;

		xLength 	= Grid->nxVy;
		iGrid0  	= Numbering->subEqSystem0Dir[1];
	} else if (Visu->type==VisuType_PRes || Visu->type==VisuType_PfRes) {
		iEqStart = Numbering->subEqSystem0[2];
		//iEqEnd   = Numbering->subEqSystem0[3];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[2];
	} else if (Visu->type==VisuType_PcRes) {
		iEqStart = Numbering->subEqSystem0[3];
		//iEqEnd   = Numbering->subEqSystem0[4];

		ixECStart 	= 0;
		ixECEnd 	= Grid->nxEC;

		iyECStart 	= 0;
		iyECEnd 	= Grid->nyEC;

		xLength 	= Grid->nxEC;
		iGrid0  	= Numbering->subEqSystem0Dir[3];
	} else {
		printf("Unauthorized Visu->type in VIsu_Residual.\n");
		exit(0);
	}


	// Could be optimized by not looping over everything (be careful to the lower trianuglar contributions though; that's why I didn't do it yet)
#pragma omp parallel for private(iEq, i, J) OMP_SCHEDULE
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

//#pragma omp parallel for private(iEq, i, J) OMP_SCHEDULE
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
