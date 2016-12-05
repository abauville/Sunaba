/*
 * Input.c

 *
 *  Created on: Jun 6, 2016
 *      Author: abauville
 */


#include "MiscLibraries/jsmn.h"
#include "stokes.h"
#define TOKEN(string) jsoneq(JSON_STRING, &t[i], string) == 0
#define VALUE(string) jsoneq(JSON_STRING, &t[i+1], string) == 0

#define NUM_TOKEN 2048

static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
	if ((int) strlen(s) == tok->end - tok->start &&
			strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
		return 0;
	}
	return -1;
}

// Define geometry types
typedef struct Circle {
	int phase;
	compute cx, cy, radius;
} Circle;
typedef struct Rect {
	int phase;
	compute llx, lly, width, height;
} Rect;
typedef struct Line {
	int phase;
	compute a, b, min, max;
	int condition, definedFor;
} Line;
typedef struct Sine {
	int phase;
	compute amplitude, base, wavelength, wavephase, min, max;
	int condition, definedFor;
} Sine;

void assignCircle(Particles* Particles, Grid* Grid, Circle* Circle);
void assignRect(Particles* Particles, Grid* Grid, Rect* Rect);
void assignLine(Particles* Particles, Grid* Grid, Line* Line);
void assignSine(Particles* Particles, Grid* Grid, Sine* Sine);

void get_ixmin_ixmax_iymin_iymax (Grid* Grid, compute coordLimits[4], int indexLimits[4]);



void Input_read(Input* Input, Grid* Grid, Numerics* Numerics, Physics* Physics, MatProps* MatProps, Particles* Particles, Char* Char, BC* BCStokes, BC* BCThermal)
{
	// ===================================================
	// 				INIT OPTIONAL VALUES
	// Init phaseAir, phaseWater, in case it is not assigned
	Physics->phaseAir = -1;
	Physics->phaseWater = -1;

	// ===================================================
	// 				INIT OPTIONAL VALUES





	// ===================================================
	// 				LOAD AND PARSE THE FILE

	char* JSON_STRING = readFile(Input->inputFile);

	int i;
	int r;
	jsmn_parser p;
	jsmntok_t t[NUM_TOKEN]; /* We expect no more than 128 tokens */

	jsmn_init(&p);
	r = jsmn_parse(&p, JSON_STRING, strlen(JSON_STRING), t, sizeof(t)/sizeof(t[0]));


	// Error checking
	if (r < 0) {
		printf("Failed to parse JSON: %d. Try increasing NUM_TOKEN\n", r);
		exit(0);
	}

	/* Assume the top-level element is an object */
	if (r < 1 || t[0].type != JSMN_OBJECT) {
		printf("Object expected\n");
		exit(0);
	}


	// 				LOAD AND PARSE THE FILE
	// ===================================================



	// Loop over all keys of the root object
	i = 1;
	int iSub,iPhaseAttr, iPhase;
	int j;
	int size, subSize;

	char* strValue = NULL; // adress where to fetch a value;
	char* strToken = NULL; // adress where to fetch a value;

	bool Stop = false;

	while (i<r) {



		//printf("Token #%i: type = %i, start = %i, end = %i, size = %i\n", i, t[i].type, t[i].start, t[i].end, t[i].size);
		//printf("Token #%i: %.*s\n",i, t[i].end-t[i].start, JSON_STRING + t[i].start);

		if (TOKEN("Numerics")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("nTimeSteps") ) {
					Numerics->nTimeSteps = atoi(strValue);
				} else if (  TOKEN("nLineSearch") ) {
					Numerics->nLineSearch = atoi(strValue);
				} else if  (  TOKEN("maxNonLinearIter") ) {
					Numerics->maxNonLinearIter = atoi(strValue);
				} else if  (  TOKEN("minNonLinearIter") ) {
					Numerics->minNonLinearIter = atoi(strValue);
				} else if  (  TOKEN("relativeTolerance") ) {
					Numerics->relativeTolerance = atof(strValue);
				} else if  (  TOKEN("absoluteTolerance") ) {
					Numerics->absoluteTolerance = atof(strValue);
				} else if  (  TOKEN("maxCorrection") ) {
					Numerics->maxCorrection = atof(strValue);
				} else if  (  TOKEN("CFL_fac_Stokes") ) {
					Numerics->CFL_fac_Stokes = atof(strValue);
				} else if  (  TOKEN("CFL_fac_Thermal") ) {
					Numerics->CFL_fac_Thermal = atof(strValue);
				} else if  (  TOKEN("CFL_fac_Darcy") ) {
					Numerics->CFL_fac_Darcy = atof(strValue);
				} else if  (  TOKEN("etaMin") ) {
					Numerics->etaMin = atof(strValue);
				} else if  (  TOKEN("etaMax") ) {
					Numerics->etaMax = atof(strValue);
				} else if  (  TOKEN("dtMin") ) {
					Numerics->dtMin = atof(strValue);
				} else if  (  TOKEN("dtMax") ) {
					Numerics->dtMax = atof(strValue);
				} else {
					printf("Unexpected key in Numerics: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}




		else if (TOKEN("Grid")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key
			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("xmin") ) {
					Grid->xmin = atof(strValue);
				} else if (  TOKEN("xmax") ) {
					Grid->xmax = atof(strValue);
				} else if  (  TOKEN("ymin") ) {
					Grid->ymin = atof(strValue);
				} else if  (  TOKEN("ymax") ) {
					Grid->ymax = atof(strValue);
				} else if  (  TOKEN("nxC") ) {
					Grid->nxC = atoi(strValue);
					if (t[i+1].size==0) {
						Grid->nSegX = 1;
						Grid->nxC = atoi(strValue);
						Grid->nxPerSeg = (int*) malloc(Grid->nSegX * sizeof(int) );
						Grid->nxPerSeg[0] = Grid->nxC;
					} else {
						printf("***************** here, size = %i\n",t[i+1].size);
						int I = i;
						Grid->nxC = 0;
						Grid->nSegX = t[i+1].size;

						Grid->nxPerSeg = (int*) malloc(Grid->nSegX * sizeof(int));

						for (j = 0; j < t[I+1].size; j++) {
							strValue = JSON_STRING+t[I+1+j+1].start;
							Grid->nxC += atoi(strValue);
							Grid->nxPerSeg[j] = atoi(strValue);
							printf("*********Value in nxC array:: %.*s\n", t[I+1+j+1].end-t[I+1+j+1].start, JSON_STRING + t[I+1+j+1].start);
							i+=1;
						}
					}
					printf("Grid->nxC = %i\n",Grid->nxC);


				} else if  (  TOKEN("nyC") ) {
					if (t[i+1].size==0) {
						Grid->nSegY = 1;
						Grid->nyPerSeg = (int*) malloc(Grid->nSegY * sizeof(int));
						Grid->nyC = atoi(strValue);
					} else  {
						printf("***************** here, size = %i\n",t[i+1].size);
						Grid->nSegY = t[i+1].size;
						Grid->nyPerSeg = (int*) malloc(Grid->nSegY * sizeof(int));
						int I = i;
						Grid->nyC = 0;
						for (j = 0; j < t[I+1].size; j++) {
							strValue = JSON_STRING+t[I+1+j+1].start;
							Grid->nyC += atoi(strValue);
							Grid->nyPerSeg[j] = atoi(strValue);
							printf("*********Value in nyC array:: %.*s\n", t[I+1+j+1].end-t[I+1+j+1].start, JSON_STRING + t[I+1+j+1].start);
							i+=1;
						}
					}
					printf("Grid->nyC = %i\n",Grid->nyC);

				} else if  (  TOKEN("xSeg") ) {
					Grid->nxC = atoi(strValue);
					if (t[i+1].size==0) {
						Grid->nSegX = 1;
						Grid->nxC = atoi(strValue);
					} else {
						printf("***************** here, size = %i\n",t[i+1].size);
						int I = i;
						Grid->nxC = 0;
						for (j = 0; j < t[I+1].size; j++) {
							strValue = JSON_STRING+t[I+1+j+1].start;
							Grid->nxC += atoi(strValue);
							printf("*********Value in nxC array:: %.*s\n", t[I+1+j+1].end-t[I+1+j+1].start, JSON_STRING + t[I+1+j+1].start);
							i+=1;
						}
					}
					printf("Grid->nxC = %i\n",Grid->nxC);
				} else if  (  TOKEN("fixedBox") ) {
					Grid->isFixed = VALUE("true"); // returns true if true, false otherwise


				} else {
					printf("Unexpected key in Grid: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}



		else if (TOKEN("Physics")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("gx") ) {
					Physics->g[0] = atof(strValue);
				} else if (  TOKEN("gy") ) {
					Physics->g[1] = atof(strValue);
				} else if  (  TOKEN("Cp") ) {
					Physics->Cp = atof(strValue);
				} else if  (  TOKEN("eta_f") ) {
#if (DARCY)
					Physics->eta_f = atof(strValue);
#endif
				} else if  (  TOKEN("rho_f") ) {
#if (DARCY)
					Physics->rho_f = atof(strValue);
#endif
				} else if  (  TOKEN("y_oceanSurface") ) {
#if (DARCY)
					Physics->y_oceanSurface = atof(strValue);
#endif
				} else {
					printf("Unexpected key in Physics: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}



		else if (TOKEN("Particles")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("noiseFactor") ) {
					Particles->noiseFactor = atof(strValue);
				} else if (  TOKEN("nPCX") ) {
					Particles->nPCX = atoi(strValue);
				} else if  (  TOKEN("nPCY") ) {
					Particles->nPCY = atoi(strValue);
				} else if  (  TOKEN("minPartPerCellFactor") ) {
					Particles->minPartPerCellFactor = atof(strValue);
				} else if  (  TOKEN("maxPartPerCellFactor") ) {
					Particles->maxPartPerCellFactor = atof(strValue);

				} else if  (  TOKEN("passiveRes") ) {
					Particles->passiveRes = atof(strValue);


				} else if  (  TOKEN("passiveGeom") ) {
					if 		  ( VALUE("Grid")) {
						Particles->passiveGeom = PartPassive_Grid;
					} else {
						printf("Unexpected Particles.passiveGeom: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}

				} else {
					printf("Unexpected key in Particles: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}


		else if (TOKEN("Char")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("length") ) {
					Char->length = atof(strValue);
				} else if (  TOKEN("mass") ) {
					Char->mass = atof(strValue);
				} else if  (  TOKEN("time") ) {
					Char->time = atof(strValue);
				} else if  (  TOKEN("temperature") ) {
					Char->temperature = atof(strValue);
				} else {
					printf("Unexpected key in Char: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}



		else if (TOKEN("BCStokes")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("backStrainRate") ) {
					BCStokes->backStrainRate = atof(strValue);
				} else if (  TOKEN("specialPhase") ) {
					BCStokes->specialPhase = atoi(strValue);
				} else if (  TOKEN("refValue") ) {
					BCStokes->refValue = atof(strValue);
				} else if (  TOKEN("DeltaL") ) {
					BCStokes->DeltaL = atof(strValue);
				} else if (  TOKEN("SetupType") ) {
					if 		  ( VALUE("PureShear")) {
						BCStokes->SetupType = Stokes_PureShear;
					} else if ( VALUE("SimpleShear")) {
						BCStokes->SetupType = Stokes_SimpleShear;
					} else if ( VALUE("FixedLeftWall")) {
						BCStokes->SetupType = Stokes_FixedLeftWall;
					} else if ( VALUE("Sandbox")) {
						BCStokes->SetupType = Stokes_Sandbox;
					} else if ( VALUE("SandboxWeakBackstop")) {
						BCStokes->SetupType = Stokes_SandboxWeakBackstop;
					} else if ( VALUE("CornerFlow")) {
						BCStokes->SetupType = Stokes_CornerFlow;
					} else {
						printf("Unexpected BCStokes.type: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						Stop = true;
					}
				} else {
					printf("Unexpected key in BCStokes: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}


		else if (TOKEN("BCThermal")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if 		  (  TOKEN("TT") ) {
					BCThermal->TT = atof(strValue);
				} else if (  TOKEN("TB") ) {
					BCThermal->TB = atof(strValue);
				} else if (  TOKEN("refValue") ) {
					BCThermal->refValue = atof(strValue);
				} else if (  TOKEN("DeltaL") ) {
					BCStokes->DeltaL = atof(strValue);
				} else if (  TOKEN("SetupType") ) {
					if 		  ( VALUE("TT_TB_LRNoFlux")) {
						BCThermal->SetupType = Thermal_TT_TB_LRNoFlux;
					} else if 		  ( VALUE("TT_TBExternal_LRNoFlux")) {
						BCThermal->SetupType = Thermal_TT_TBExternal_LRNoFlux;
					} else {
						printf("Unexpected BCThermal.SetupType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						Stop = true;
					}

				} else {
					printf("Unexpected key in BCThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
		}


















		else if (TOKEN("MatProps")) {
			i++; // Move to the first token, which is the object
			MatProps->nPhase = t[i].size; // number of elements in the token
			//printf("nPhase = %i\n",t[i].size);
			i++; // Move to the first key


			for (iSub=0; iSub<MatProps->nPhase; iSub++) {
				strToken = JSON_STRING+t[i].start;
				iPhase = atoi(strToken);
				//printf("PhaseNumber = %i\n",PhaseNumber);
				//printf("Phase token: %.*s\n", t[i].end-t[i].start, strToken);


				i++;
				subSize = t[i].size; // number of elements in the token
				i++;
				for (iPhaseAttr=0; iPhaseAttr<subSize; iPhaseAttr++) {
					strToken = JSON_STRING+t[i].start;
					strValue = JSON_STRING+t[i+1].start;

					if 		  (  TOKEN("n") ) {
						MatProps->n[iPhase] = atof(strValue);
					} else if (  TOKEN("eta0") ) {
						MatProps->eta0[iPhase] = atof(strValue);
					} else if (  TOKEN("rho0") ) {
						MatProps->rho0[iPhase] = atof(strValue);
					} else if (  TOKEN("frictionAngle") ) {
						MatProps->frictionAngle[iPhase] = atof(strValue);
					} else if (  TOKEN("cohesion") ) {
						MatProps->cohesion[iPhase] = atof(strValue);

					} else if (  TOKEN("G") ) {
						MatProps->G[iPhase] = atof(strValue);

					} else if (  TOKEN("k") ) {
						MatProps->k[iPhase] = atof(strValue);
					} else if (  TOKEN("alpha") ) {
						MatProps->alpha[iPhase] = atof(strValue);
					} else if (  TOKEN("beta") ) {
						MatProps->beta[iPhase] = atof(strValue);
					} else if (  TOKEN("perm0") ) {
						MatProps->perm0[iPhase] = atof(strValue);
					} else if (  TOKEN("eta_b") ) {
						MatProps->eta_b[iPhase] = atof(strValue);
					} else if (  TOKEN("B") ) {
						MatProps->B[iPhase] = atof(strValue);
					} else if (  TOKEN("material") ) {
						// Place holder
					} else if (  TOKEN("name") ) {
						// Place holder
					} else if  (  TOKEN("isAir") ) {
						MatProps->isAir[iPhase] = VALUE("true"); // returns true if true, false otherwise
						if (VALUE("true"))
							Physics->phaseAir = iPhase;
					} else if  (  TOKEN("isWater") ) {
						MatProps->isWater[iPhase] = VALUE("true"); // returns true if true, false otherwise
						if (VALUE("true"))
							Physics->phaseWater = iPhase;
					} else if  (  TOKEN("isRef") ) {
						if (VALUE("true"))
							Physics->phaseRef = iPhase;
					} else {
						printf("Unexpected key in MatProps: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
						Stop = true;
					}


					//printf("Inside phase token: %.*s\n", t[i].end-t[i].start, strToken);
					i+=2;
				}


				//i+=2;
			}
		}


		else if (TOKEN("Visu")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key
			//printf("In Visu: size = %i\n",size);
			for (iSub=0; iSub<size; iSub++) {
				i+=size;
			}
		}

		else if (TOKEN("Geometry")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key
			//printf("In Visu: size = %i\n",size);
			for (iSub=0; iSub<size; iSub++) {

				i++;
				subSize = t[i].size; // number of elements in the token
				i++;
				int iSubSub;
				for (iSubSub=0; iSubSub<subSize; iSubSub++) {
					i+=2;
				}
				//i+=2;
			}
		}

		else if (TOKEN("Description")) {
			i++; // Move to the first token, which is the object
			i++; // Move to the first key
		}


		else {
			printf("Unexpected key: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
			Stop = true;
			i++;

		}

		//i++;


	}



	if (Stop) {
		printf("error: unexpected keys were found. Some variables could not be assigned values");
		exit(0);
	}



	free(JSON_STRING);


}


#if (VISU)

void Input_readVisu(Input* Input, Visu* Visu)
{
	// ===================================================
	// 				LOAD AND PARSE THE FILE

	char* JSON_STRING = readFile(Input->inputFile);

	int i;
	int r;
	jsmn_parser p;
	jsmntok_t t[NUM_TOKEN]; /* We expect no more than 128 tokens */

	jsmn_init(&p);
	r = jsmn_parse(&p, JSON_STRING, strlen(JSON_STRING), t, sizeof(t)/sizeof(t[0]));


	// Error checking
	if (r < 0) {
		printf("Failed to parse JSON: %d\n", r);
		exit(0);
	}

	/* Assume the top-level element is an object */
	if (r < 1 || t[0].type != JSMN_OBJECT) {
		printf("Object expected\n");
		exit(0);
	}


	// 				LOAD AND PARSE THE FILE
	// ===================================================



	// Loop over all keys of the root object
	i = 1;
	int iSub;
	int size;

	char* strValue = NULL; // adress where to fetch a value;

	bool FoundVisu = false;
	bool Stop = false;


	while (i<r) {
		if (TOKEN("Visu")) {
			FoundVisu = true;
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				printf("This key is: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
				strValue = JSON_STRING+t[i+1].start;
				if 		  (  TOKEN("shiftFacX") ) {
					Visu->shiftFac[0] = atof(strValue);
				} else if (  TOKEN("shiftFacY") ) {
					Visu->shiftFac[1] = atof(strValue);
				} else if (  TOKEN("shiftFacZ") ) {
					Visu->shiftFac[2] = atof(strValue);
				} else if  (  TOKEN("glyphSamplingRateX") ) {
					Visu->glyphSamplingRateX = atoi(strValue);
				} else if  (  TOKEN("glyphSamplingRateY") ) {
					Visu->glyphSamplingRateY = atoi(strValue);

				} else if  (  TOKEN("alphaOnValue") ) {
					Visu->alphaOnValue = VALUE("true"); // returns true if true, false otherwise
				} else if  (  TOKEN("width") ) {
					Visu->width = atoi(strValue);
				} else if  (  TOKEN("height") ) {
					Visu->height = atoi(strValue);
				} else if  (  TOKEN("particleMeshRes") ) {
					Visu->particleMeshRes = atoi(strValue);
				} else if  (  TOKEN("particleMeshSize") ) {
					Visu->particleMeshSize = atof(strValue);
				} else if  (  TOKEN("glyphScale") ) {
					Visu->glyphScale = atof(strValue);
				} else if  (  TOKEN("showParticles") ) {
					Visu->showParticles = VALUE("true");
				} else if  (  TOKEN("showGlyphs") ) {
					Visu->showGlyphs = VALUE("true");
				} else if  (  TOKEN("writeImages") ) {
					Visu->writeImages = VALUE("true");
				} else if  (  TOKEN("transparency") ) {
					Visu->transparency = VALUE("true");
				} else if  (  TOKEN("outputFolder") ) {
					if (t[i+1].end-t[i+1].start>MAX_STRING_LENGTH) {
						printf("the Visu.outputFolder string is too long, maximum authorized: %i. Please change your folder or increase the value of the macro MAX_STRING_LENGTH", MAX_STRING_LENGTH);
					}
					strncpy(Visu->outputFolder, strValue, t[i+1].end-t[i+1].start);
				} else if  (  TOKEN("retinaScale") ) {
					Visu->retinaScale = atoi(strValue);








				} else if  (  TOKEN("glyphType") ) {

					if 		  ( VALUE("StokesVelocity")) {
						Visu->glyphType = StokesVelocity;
					} else if ( VALUE("DarcyGradient")) {
						Visu->glyphType = StokesVelocity;
					} else {
						printf("Unexpected Visu.glyphType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}

				} else if  (  TOKEN("filter") ) {

					if 		  ( VALUE("Linear")) {
						Visu->filter = Linear;
					} else if ( VALUE("Nearest")) {
						Visu->filter = Nearest;
					} else {
						printf("Unexpected Visu.glyphType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}





				} else if  (  TOKEN("glyphMeshType") ) {
					if 		  ( VALUE("Triangle")) {
						Visu->glyphMeshType = Triangle;
					} else if ( VALUE("ThinArrow")) {
						Visu->glyphMeshType = ThinArrow;
					} else if ( VALUE("ThickArrow")) {
						Visu->glyphMeshType = ThickArrow;
					} else {
						printf("Unexpected Visu.particleType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}


				} else if  (  TOKEN("typeParticles") ) {
					if 		  ( VALUE("PartPhase")) {
						Visu->typeParticles = PartPhase;
					} else if ( VALUE("PartTemp")) {
						Visu->typeParticles = PartTemp;
					} else if ( VALUE("PartSigma_xx")) {
						Visu->typeParticles = PartSigma_xx;
					} else if ( VALUE("PartSigma_xy")) {
						Visu->typeParticles = PartSigma_xy;
					} else {
						printf("Unexpected Visu.particleType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}




				} else if  (  TOKEN("type") ) {
					// not use, but kept to make the input file human readable
				} else if (  TOKEN("typeNumber") ) {
					Visu->type = atoi(strValue);




				} else if (  TOKEN("colorMap") ) {
					//printf("Unexpected Visu.particleType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
					// This is the name
					i++; // Move to the first token, which is the object
					int size2 = t[i].size; // number of elements in the token
					int iSub2;
					i++; // Move to the first key

					//printf("size2 = %i\n",size2);
					for (iSub2=0; iSub2<size2; iSub2++) {
						//printf("00 This key: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
						i++; // Move to the first token, which is the object
						int size3 = t[i].size; // number of elements in the token
						//printf("size3 = %i\n",size3);
						int iSub3;
						i++; // Move to the first key
						int thisType = 0;
						for (iSub3=0; iSub3<size3; iSub3++) {
							//printf("iSub3 =%i\n",iSub3);
							//printf("This sub key: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							strValue = JSON_STRING+t[i+1].start;

							if 			(  TOKEN("a0number") ) {
								thisType = atoi(strValue);
							} else if	(  TOKEN("center") ) {
								Visu->colorMap[thisType].center = atof(strValue);
								printf("thisType = %i\n",thisType);
							} else if	(  TOKEN("colorMapRes") ) {
								Visu->colorMap[thisType].colorMapRes = atoi(strValue);
							} else if	(  TOKEN("colorMap") ) {
								Visu->colorMap[thisType].colorMap = 0.0;// Dummy for the moment
							} else if	(  TOKEN("log10on") ) {
								Visu->colorMap[thisType].log10on = VALUE("true");
							} else if	(  TOKEN("max") ) {
								Visu->colorMap[thisType].max = atof(strValue);
							} else if	(  TOKEN("scale") ) {
								Visu->colorMap[thisType].scale = atof(strValue);
							} else if	(  TOKEN("type") ) {
								// Dummy - not implemented yet
							} else {
								printf("Unexpected Visu.colorMap attribute: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
								exit(0);
							}
							i+=2;
						}
						//i+=2;
					}

					i-=2;


				} else {
					printf("Unexpected Key in Visu: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}
				i+=2;
			}
		}

		else {
			i++;
		}

	}


	if (!FoundVisu) {
		printf("error: Visu was not found in the input file\n");
		exit(0);
	}

	if (Stop) {
		printf("error: unexpected keys were found. Some variables could not be assigned values");
		exit(0);
	}


	free(JSON_STRING);

}

#endif





void Input_assignPhaseToParticles(Input* Input, Particles* Particles, Grid* Grid, Char* Char)
{
	// ===================================================
	// 				LOAD AND PARSE THE FILE

	char* JSON_STRING = readFile(Input->inputFile);

	int i;
	int r;
	jsmn_parser p;
	jsmntok_t t[NUM_TOKEN]; /* We expect no more than 128 tokens */

	jsmn_init(&p);
	r = jsmn_parse(&p, JSON_STRING, strlen(JSON_STRING), t, sizeof(t)/sizeof(t[0]));


	// Error checking
	if (r < 0) {
		printf("Failed to parse JSON: %d\n", r);
		exit(0);
	}

	/* Assume the top-level element is an object */
	if (r < 1 || t[0].type != JSMN_OBJECT) {
		printf("Object expected\n");
		exit(0);
	}


	// 				LOAD AND PARSE THE FILE
	// ===================================================



	// Loop over all keys of the root object
	i = 1;
	int iSub, iSubSub;
	int size, subSize;

	char* strValue = NULL; // adress where to fetch a value;
	//char* strToken = NULL; // adress where to fetch a value;

	bool FoundGeometry = false;


	bool Stop = false;

	// Variables for the different types of geometry
	Circle circle;
	// Rect
	Rect rect;
	// Line
	Line line;
	// Sine
	Sine sine;



	while (i<r) {




		if (TOKEN("Geometry")) {
			FoundGeometry = true;
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			//printf("nPhase = %i\n",t[i].size);
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				t[i].start += 6;
				//strToken = JSON_STRING+t[i].start;

				//iPhase = atoi(strToken);
				//printf("PhaseNumber = %i\n",PhaseNumber);
				//printf("Phase token: %.*s\n", t[i].end-t[i].start, strToken);
				//printf("Token: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);


				if 		  (  TOKEN("circle") ) {
					//printf("Found circle\n");
					// Get attributes
					// ========================
					i++;
					subSize = t[i].size; // number of elements in the token
					i++;
					for (iSubSub=0; iSubSub<subSize; iSubSub++) {

						//strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							circle.phase = atoi(strValue);
						} else if (  TOKEN("cx") ) {
							circle.cx = atof(strValue)/Char->length;
						} else if (  TOKEN("cy") ) {
							circle.cy = atof(strValue)/Char->length;
						} else if (  TOKEN("radius") ) {
							circle.radius = atof(strValue)/Char->length;
						} else {
							printf("Unexpected Key in Geometry circle: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							Stop = true;
						}
						i+=2;
					}

					assignCircle(Particles, Grid, &circle);

				} else if (  TOKEN("rect") ) {
					//printf("Found rect\n");
					// Get attributes
					// ========================
					i++;
					subSize = t[i].size; // number of elements in the token
					i++;
					for (iSubSub=0; iSubSub<subSize; iSubSub++) {

						//strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							rect.phase = atoi(strValue);
						} else if (  TOKEN("llx") ) {
							rect.llx = atof(strValue)/Char->length;
						} else if (  TOKEN("lly") ) {
							rect.lly = atof(strValue)/Char->length;
						} else if (  TOKEN("width") ) {
							rect.width = atof(strValue)/Char->length;
						} else if (  TOKEN("height") ) {
							rect.height = atof(strValue)/Char->length;
						} else {
							printf("Unexpected Key in Geometry rect: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							Stop = true;
						}
						i+=2;
					}

					assignRect(Particles, Grid, &rect);

				} else if (  TOKEN("line") ) {
					//printf("Found line\n");
					// Get attributes
					// ========================
					i++;
					subSize = t[i].size; // number of elements in the token
					i++;
					for (iSubSub=0; iSubSub<subSize; iSubSub++) {

						//strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							line.phase = atoi(strValue);
						} else if (  TOKEN("a") ) {
							line.a = atof(strValue);
						} else if (  TOKEN("b") ) {
							line.b = atof(strValue)/Char->length;
						} else if (  TOKEN("min") ) {
							line.min = atof(strValue)/Char->length;
						} else if (  TOKEN("max") ) {
							line.max = atof(strValue)/Char->length;
						} else if (  TOKEN("definedFor") ) {
							if (VALUE("x")) {
								line.definedFor = 0;
							} else if (VALUE("y")) {
								line.definedFor = 1;
							} else {
								printf("Unexpected Value for line.definedFor: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
								Stop = true;
							}
						} else if (  TOKEN("condition") ) {
							if (VALUE("<")) {
								line.condition = 0;
							} else if (VALUE(">")) {
								line.condition = 1;
							} else {
								printf("Unexpected Value for line.condition: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
								Stop = true;
							}
						} else {
							printf("Unexpected Key in Geometry line: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							Stop = true;
						}

						i+=2;
					}

					assignLine(Particles, Grid, &line);

				} else if (  TOKEN("sine") ) {
					//printf("Found sine\n");
					i++;
					subSize = t[i].size; // number of elements in the token
					i++;
					for (iSubSub=0; iSubSub<subSize; iSubSub++) {

						//strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							sine.phase = atoi(strValue);
						} else if (  TOKEN("amplitude") ) {
							sine.amplitude = atof(strValue)/Char->length;
						} else if (  TOKEN("base") ) {
							sine.base = atof(strValue)/Char->length;
						} else if (  TOKEN("min") ) {
							sine.min = atof(strValue)/Char->length;
						} else if (  TOKEN("max") ) {
							sine.max = atof(strValue)/Char->length;
						} else if (  TOKEN("wavelength") ) {
							sine.wavelength = atof(strValue)/Char->length;
						} else if (  TOKEN("wavephase") ) {
							sine.wavephase = atof(strValue);
						} else if (  TOKEN("definedFor") ) {
							if (VALUE("x")) {
								sine.definedFor = 0;
							} else if (VALUE("y")) {
								sine.definedFor = 1;
							} else {
								printf("Unexpected Value for line.definedFor: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
								Stop = true;
							}
						} else if (  TOKEN("condition") ) {
							if (VALUE("<")) {
								sine.condition = 0;
							} else if (VALUE(">")) {
								sine.condition = 1;
							} else {
								printf("Unexpected Value for line.condition: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
								Stop = true;
							}
						} else {
							printf("Unexpected Key in Geometry sine: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							Stop = true;
						}

						i+=2;
					}
					assignSine(Particles, Grid, &sine);


				} else if (  TOKEN("polygon") ) {
					printf("error: polygons are not supported yet\n");
					exit(0);
				} else {
					printf("Unexpected Key in Geometry: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					i+=2;
					Stop = true;
				}


				//i+=2;
			}
		}
		else {
			i++;
		}
	}



	if (!FoundGeometry) {
		printf("error: no Geometry section was found in the input file\n");
		exit(0);
	}

	if (Stop) {
		printf("error: unexpected keys were found. Some variables could not be assigned values");
		exit(0);
	}

	free(JSON_STRING);

	//printf("End of assign\n");
	//exit(0);
}

void assignCircle(Particles* Particles, Grid* Grid, Circle* Circle)
{
	printf("In assign Circle\n");
	int ix, iy;
	compute x, y;
	compute r2 = Circle->radius*Circle->radius;

	// Get the bounding box
	compute xmin, xmax, ymin,ymax;
	xmin = Circle->cx - Circle->radius;
	xmax = Circle->cx + Circle->radius;
	ymin = Circle->cy - Circle->radius;
	ymax = Circle->cy + Circle->radius;

	compute coordLimits[4] = {xmin,xmax,ymin,ymax};
	int indexLimits[4];
	get_ixmin_ixmax_iymin_iymax (Grid, coordLimits, indexLimits);



	//printf("init ok assign Circle, CirclePhase = %i, iymin = %i, iymax = %i, ixmin = %i, ixmax = %i, xmin = %.3f, Grid->xmin = %.3f\n",Circle->phase, iymin, iymax, ixmin, ixmax, xmin, Grid->xmin);
	INIT_PARTICLE

	for (iy = indexLimits[2]; iy < indexLimits[3]; ++iy) {
		for (ix = indexLimits[0]; ix < indexLimits[1]; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle != NULL) {
				x = thisParticle->x-Circle->cx;
				y = thisParticle->y-Circle->cy;
				//if (sqrDistance < sqrRadius) {
				if ( ( (x*x)/(r2) + (y*y)/(r2) )<= 1 ) { // note: infinity shape condition: sqrDistance<rx*cos(alpha)*ry*sin(alpha)
					//printf("In!!\n");
					thisParticle->phase = Circle->phase;
				}


				thisParticle = thisParticle->next;
			}
		}
	}

	printf("out of assign Circle\n");
}



void assignRect(Particles* Particles, Grid* Grid, Rect* Rect) {
	printf("In assign Rect\n");
	int ix, iy;
	compute x, y;


	compute coordLimits[4] = {Rect->llx,Rect->llx+Rect->width,Rect->lly,Rect->lly+Rect->height};
	int indexLimits[4];
	get_ixmin_ixmax_iymin_iymax (Grid, coordLimits, indexLimits);





	//printf("ixmin = %i, ixmax = %i, iymin = %i, iymax = %i, nxS = %i\n", ixmin,ixmax, iymin, iymax, Grid->nxS);

	INIT_PARTICLE

	for (iy = indexLimits[2]; iy < indexLimits[3]; ++iy) {
		for (ix = indexLimits[0]; ix < indexLimits[1]; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle != NULL) {
				x = thisParticle->x-Rect->llx;
				y = thisParticle->y-Rect->lly;
				//if (sqrDistance < sqrRadius) {
				if ( x>=0 && x<=Rect->width && y>=0 && y<=Rect->height) { // note: infinity shape condition: sqrDistance<rx*cos(alpha)*ry*sin(alpha)
					//printf("In!!\n");
					thisParticle->phase = Rect->phase;
				}


				thisParticle = thisParticle->next;
			}
		}
	}

	printf("out of assign Rect\n");
}



void assignLine(Particles* Particles, Grid* Grid, Line* Line) {
	printf("In assign Line\n");
	int ix, iy;
	compute x, y;
	compute a = Line->a;
	compute b = Line->b;

	//printf("A\n");
	compute xmin,xmax, ymin, ymax;
	if (Line->definedFor == 1) {
		xmin = Line->min;
		xmax = Line->max;
		if (a>0) {
			if (Line->condition == 1) {
				ymin = a*xmin + b;
				ymax = Grid->ymax;
			}
			else if (Line->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xmax + b;
			}

		} else {
			if (Line->condition == 1) {
				ymin = a*xmax + b;
				ymax = Grid->ymax;
			}
			else if (Line->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xmin + b;
			}

		}
	} else if (Line->definedFor == 0) {
		ymin = Line->min;
		ymax = Line->max;
		if (a>0) {
			if (Line->condition == 1) {
				xmin = a*ymin + b;
				xmax = Grid->xmax;
			}
			else if (Line->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*ymax + b;
			}

		} else {
			if (Line->condition == 1) {
				xmin = a*ymax + b;
				xmax = Grid->xmax;
			}
			else if (Line->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*ymin + b;
			}
		}
	}

	//printf("B\n");

	compute coordLimits[4] = {xmin,xmax,ymin,ymax};
	int indexLimits[4];
	get_ixmin_ixmax_iymin_iymax (Grid, coordLimits, indexLimits);

	//printf("C\n");



	INIT_PARTICLE

	for (iy = indexLimits[2]; iy < indexLimits[3]; ++iy) {
		for (ix = indexLimits[0]; ix < indexLimits[1]; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;

			thisParticle = Particles->linkHead[iNode];
			//printf("D\n");
			while (thisParticle != NULL) {
				//printf("E\n");
				//printf("Line->phase = %i\n", Line->phase);
				//printf("thisPart->phase = %i\n", thisParticle->phase);
				x = thisParticle->x;
				y = thisParticle->y;
				//if (sqrDistance < sqrRadius) {
				if (Line->definedFor == 1) {
					if ( Line->condition == 1 ) { // >
						//printf(">\n");
						if ( y > a*x + b) {
							thisParticle->phase = Line->phase;
						}
					} else if ( Line->condition == 0 ) { // <
						//printf("<\n");
						if ( y < a*x + b) {
							thisParticle->phase = Line->phase;
						}
					}
				} else if (Line->definedFor == 0) {
					if ( Line->condition == 1 ) { // >
						if ( x > a*y + b) {
							thisParticle->phase = Line->phase;
						}
					} else if ( Line->condition == 0 ) { // <
						if ( x < a*y + b) {
							thisParticle->phase = Line->phase;
						}
					}
				}



				thisParticle = thisParticle->next;
			}
		}
	}

	printf("out of assign Line\n");
}
void assignSine(Particles* Particles, Grid* Grid, Sine* Sine) {
	printf("In assign Sine\n");
	int ix, iy;
	compute x, y;


	compute xmin,xmax, ymin, ymax;

	if (Sine->definedFor == 1) {
		xmin = Sine->min;
		xmax = Sine->max;
		if (Sine->condition == 1) {
			ymin = Sine->base-Sine->amplitude;
			ymax = Grid->ymax;
		}
		else if (Sine->condition == 0) {
			ymax = Sine->base+Sine->amplitude;
			ymin = Grid->ymin;
		}

	} else if (Sine->definedFor == 0) {
		ymin = Sine->min;
		ymax = Sine->max;
		if (Sine->condition == 1) {
			xmin = Sine->base-Sine->amplitude;
			xmax = Grid->xmax;
		}
		else if (Sine->condition == 0) {
			xmax = Sine->base+Sine->amplitude;
			xmin = Grid->xmin;
		}
	}


	compute coordLimits[4] = {xmin,xmax,ymin,ymax};
	int indexLimits[4];
	get_ixmin_ixmax_iymin_iymax (Grid, coordLimits, indexLimits);



	INIT_PARTICLE

	for (iy = indexLimits[2]; iy < indexLimits[3]; ++iy) {
		for (ix = indexLimits[0]; ix < indexLimits[1]; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle != NULL) {

				//if (sqrDistance < sqrRadius) {
				if (Sine->definedFor == 1) {
					x = (thisParticle->x-Grid->xmin);///(Grid->xmax-Grid->xmin);
					y = thisParticle->y;
					if ( Sine->condition == 1 ) { // >
						if ( y > Sine->base + Sine->amplitude*sin(1.0/Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( y < Sine->base + Sine->amplitude*sin(1.0/Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					}
				} else if (Sine->definedFor == 0) {
					x = thisParticle->x;///(Grid->xmax-Grid->xmin);
					y = thisParticle->y-Grid->ymin;
					if ( Sine->condition == 1 ) { // >
						if ( x > Sine->base + Sine->amplitude*sin(1.0/Sine->wavelength*y*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( x < Sine->base + Sine->amplitude*sin(1.0/Sine->wavelength*y*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					}
				}



				thisParticle = thisParticle->next;
			}
		}
	}

	printf("out of assign Sine\n");

}


void get_ixmin_ixmax_iymin_iymax (Grid* Grid, compute coordLimits[4], int indexLimits[4])
{
	// coordLimits[4] = {xmin, xmax, ymin, ymax};
	// indexlimits[4] = {ixmin, ixmax. iymin, iymax}

	int ix, iy;
	ix = 0;
	while (ix < Grid->nxS-1 && Grid->X[ix]<coordLimits[0]  ) {
		ix++;
	}
	if (ix!=0)
		ix = ix-1;
	indexLimits[0] = ix;
	while (ix < Grid->nxS-1 && Grid->X[ix]<coordLimits[1]) {
		ix++;
	}
	//if (ix==Grid->nxS)
	//	ix = ix-1;

	indexLimits[1] = ix+1;

	iy = 0;
	while (iy < Grid->nyS-1 && Grid->Y[iy]<coordLimits[2]) {
		iy++;
	}
	if (iy!=0)
	iy = iy-1;
	indexLimits[2] = iy;
	while (iy < Grid->nyS-1 && Grid->Y[iy]<coordLimits[3]) {
		iy++;
	}
	//if (iy==Grid->nyS)
	//	iy = iy-1;
	indexLimits[3] = iy+1;



}

