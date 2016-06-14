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





void Input_read(Input* Input, Grid* Grid, Numerics* Numerics, Physics* Physics, MatProps* MatProps, Particles* Particles, Char* Char, BC* BCStokes, BC* BCThermal)
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
				} else if  (  TOKEN("CFL_fac") ) {
					Numerics->CFL_fac = atof(strValue);
				} else if  (  TOKEN("etaMin") ) {
					Numerics->etaMin = atof(strValue);
				} else if  (  TOKEN("etaMax") ) {
					Numerics->etaMax = atof(strValue);
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
				} else if  (  TOKEN("nyC") ) {
					Grid->nyC = atoi(strValue);
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
				} else if (  TOKEN("SetupType") ) {
					if 		  ( VALUE("PureShear")) {
						BCStokes->SetupType = PureShear;
					} else if ( VALUE("SimpleShearPeriodic")) {
						BCStokes->SetupType = SimpleShearPeriodic;
					} else if ( VALUE("FixedLeftWall")) {
						BCStokes->SetupType = FixedLeftWall;
					} else if ( VALUE("SandBox")) {
						BCStokes->SetupType = Sandbox;
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
				} else if (  TOKEN("SetupType") ) {
					if 		  ( VALUE("PureShear")) {
						BCThermal->SetupType = PureShear;
					} else if ( VALUE("SimpleShearPeriodic")) {
						BCThermal->SetupType = SimpleShearPeriodic;
					} else if ( VALUE("FixedLeftWall")) {
						BCThermal->SetupType = FixedLeftWall;
					} else if ( VALUE("SandBox")) {
						BCThermal->SetupType = Sandbox;
					} else {
						printf("Unexpected BCThermal.type: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
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
					} else if (  TOKEN("material") ) {
						// Place holder
					} else if (  TOKEN("name") ) {
						// Place holder
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
				i+=2;
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
				} else if  (  TOKEN("retinaScale") ) {
					Visu->glyphSamplingRateY = atoi(strValue);
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
					if 		  ( VALUE("Phase")) {
						Visu->typeParticles = Phase;
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
					if 		  ( VALUE("Blank")) {
						Visu->type = Blank;
					} else if ( VALUE("Viscosity")) {
						Visu->type = Viscosity;
					} else if ( VALUE("StrainRate")) {
						Visu->type = StrainRate;
					} else if ( VALUE("Velocity")) {
						Visu->type = Velocity;
					} else if ( VALUE("Pressure")) {
						Visu->type = Pressure;
					} else if ( VALUE("Density")) {
						Visu->type = Density;
					} else if ( VALUE("Temperature")) {
						Visu->type = Temperature;
					} else if ( VALUE("Stress")) {
						Visu->type = Stress;
					} else if ( VALUE("WaterPressureHead")) {
						Visu->type = WaterPressureHead;
					} else if ( VALUE("Permeability")) {
						Visu->type = Permeability;
					} else {
						printf("Unexpected Visu.type: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						Stop = true;
					}

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





void Input_assignPhaseToParticles(Input* Input, Particles* Particles, Grid* Grid)
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
	char* strToken = NULL; // adress where to fetch a value;

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
				strToken = JSON_STRING+t[i].start;

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

						strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							circle.phase = atoi(strValue);
						} else if (  TOKEN("cx") ) {
							circle.cx = atof(strValue);
						} else if (  TOKEN("cy") ) {
							circle.cy = atof(strValue);
						} else if (  TOKEN("radius") ) {
							circle.radius = atof(strValue);
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

						strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							rect.phase = atoi(strValue);
						} else if (  TOKEN("llx") ) {
							rect.llx = atof(strValue);
						} else if (  TOKEN("lly") ) {
							rect.lly = atof(strValue);
						} else if (  TOKEN("width") ) {
							rect.width = atof(strValue);
						} else if (  TOKEN("height") ) {
							rect.height = atof(strValue);
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

						strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							line.phase = atoi(strValue);
						} else if (  TOKEN("a") ) {
							line.a = atof(strValue);
						} else if (  TOKEN("b") ) {
							line.b = atof(strValue);
						} else if (  TOKEN("min") ) {
							line.min = atof(strValue);
						} else if (  TOKEN("max") ) {
							line.max = atof(strValue);
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

						strToken = JSON_STRING+t[i].start;
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("phase") ) {
							sine.phase = atoi(strValue);
						} else if (  TOKEN("amplitude") ) {
							sine.amplitude = atof(strValue);
						} else if (  TOKEN("base") ) {
							sine.base = atof(strValue);
						} else if (  TOKEN("min") ) {
							sine.min = atof(strValue);
						} else if (  TOKEN("max") ) {
							sine.max = atof(strValue);
						} else if (  TOKEN("wavelength") ) {
							sine.wavelength = atof(strValue);
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

	int ixmin, ixmax, iymin, iymax;
	ixmin = floor((xmin-Grid->xmin)/Grid->dx);
	ixmax = ceil((xmax-Grid->xmin)/Grid->dx);
	iymin = floor((ymin-Grid->ymin)/Grid->dy);
	iymax = ceil((ymax-Grid->ymin)/Grid->dy);

	printf("init ok assign Circle, CirclePhase = %i, iymin = %i, iymax = %i, ixmin = %i, ixmax = %i, xmin = %.3f, Grid->xmin = %.3f\n",Circle->phase, iymin, iymax, ixmin, ixmax, xmin, Grid->xmin);
	INIT_PARTICLE

	for (iy = iymin; iy < iymax; ++iy) {
		for (ix = ixmin; ix < ixmax; ++ix) {
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


	int ixmin, ixmax, iymin, iymax;
	ixmin = floor((Rect->llx-Grid->xmin)/Grid->dx);
	ixmax = ceil(((Rect->llx+Rect->width)-Grid->xmin)/Grid->dx);
	iymin = floor((Rect->lly-Grid->ymin)/Grid->dy);
	iymax = ceil(((Rect->llx+Rect->height)-Grid->ymin)/Grid->dy);

	INIT_PARTICLE

	for (iy = iymin; iy < iymax; ++iy) {
		for (ix = ixmin; ix < ixmax; ++ix) {
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

	printf("out of assign Circle\n");
}



void assignLine(Particles* Particles, Grid* Grid, Line* Line) {
	printf("In assign Line\n");
	int ix, iy;
	compute x, y;
	compute a = Line->a;
	compute b = Line->b;

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



	int ixmin, ixmax, iymin, iymax;
	ixmin = floor((xmin-Grid->xmin)/Grid->dx);
	ixmax = floor((xmax-Grid->xmin)/Grid->dx)+1;
	iymin = floor((ymin-Grid->ymin)/Grid->dy);
	iymax = floor((ymax-Grid->ymin)/Grid->dy)+1;



	INIT_PARTICLE

	for (iy = iymin; iy < iymax; ++iy) {
		for (ix = ixmin; ix < ixmax; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle != NULL) {
				x = thisParticle->x;
				y = thisParticle->y;
				//if (sqrDistance < sqrRadius) {
				if (Line->definedFor == 1) {
					if ( Line->condition == 1 ) { // >
						if ( y > a*x + b) {
							thisParticle->phase = Line->phase;
						}
					} else if ( Line->condition == 0 ) { // <
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
	printf("In assign Line\n");
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



	int ixmin, ixmax, iymin, iymax;
	ixmin = floor((xmin-Grid->xmin)/Grid->dx);
	ixmax = floor((xmax-Grid->xmin)/Grid->dx)+1;
	iymin = floor((ymin-Grid->ymin)/Grid->dy);
	iymax = floor((ymax-Grid->ymin)/Grid->dy)+1;



	INIT_PARTICLE

	for (iy = iymin; iy < iymax; ++iy) {
		for (ix = ixmin; ix < ixmax; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle != NULL) {
				x = thisParticle->x;
				y = thisParticle->y;
				//if (sqrDistance < sqrRadius) {
				if (Sine->definedFor == 1) {
					if ( Sine->condition == 1 ) { // >
						if ( y > Sine->base + Sine->amplitude*sin(Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( y < Sine->base + Sine->amplitude*sin(Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					}
				} else if (Sine->definedFor == 0) {
					if ( Sine->condition == 1 ) { // >
						if ( x > Sine->base + Sine->amplitude*sin(Sine->wavelength*y*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( x > Sine->base + Sine->amplitude*sin(Sine->wavelength*y*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					}
				}



				thisParticle = thisParticle->next;
			}
		}
	}

	printf("out of assign Line\n");

}

