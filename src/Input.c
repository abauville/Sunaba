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
						exit(0);
					}

                } else {
                    printf("Unexpected key in BCStokes: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
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
						exit(0);
					}

                } else {
                    printf("Unexpected key in BCSThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
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
            printf("In Visu: size = %i\n",size);
        	for (iSub=0; iSub<size; iSub++) {
                i+=2;
            }
        }



        else {
        	printf("Unexpected key: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
        	i++;

        }

        //i++;


    }










}




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
						Visu->type = Triangle;
					} else if ( VALUE("ThinArrow")) {
						Visu->type = ThinArrow;
					} else if ( VALUE("ThickArrow")) {
						Visu->type = ThickArrow;
					} else {
						printf("Unexpected Visu.particleType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}


				} else if  (  TOKEN("typeParticles") ) {
					if 		  ( VALUE("Phase")) {
						Visu->type = Phase;
					} else if ( VALUE("PartTemp")) {
						Visu->type = PartTemp;
					} else if ( VALUE("PartSigma_xx")) {
						Visu->type = PartSigma_xx;
					} else if ( VALUE("PartSigma_xy")) {
						Visu->type = PartSigma_xy;
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
						exit(0);
					}

				} else {
					printf("Unexpected Key in Visu: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
				}
				i+=2;
			}
		}

		else {
			i++;
		}

	}


	if (FoundVisu == false) {
		printf("error: Visu was not found in the input file\n");
		exit(0);
	}

}


