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

#define NUM_TOKEN 4096

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
	compute amplitude, base, wavelength, wavephase, min, max, slope;
	int condition, definedFor;
} Sine;

void assignCircle(Particles* Particles, Grid* Grid, Circle* Circle);
void assignRect(Particles* Particles, Grid* Grid, Rect* Rect);
void assignLine(Particles* Particles, Grid* Grid, Line* Line);
void assignSine(Particles* Particles, Grid* Grid, Sine* Sine);

void get_ixmin_ixmax_iymin_iymax (Grid* Grid, compute coordLimits[4], int indexLimits[4]);



void Input_read(Model* Model)
{
	// Declare structures
	// =================================
	Input* Input 			= &(Model->Input);
	Grid* Grid 				= &(Model->Grid);
	Numerics* Numerics 		= &(Model->Numerics);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Char* Char 				= &(Model->Char);
	BC* BCStokes 			= &(Model->BCStokes);
	BC* BCThermal 			= &(Model->BCThermal);
	IC* ICDarcy 			= &(Model->ICDarcy);
	IC* ICThermal 			= &(Model->ICThermal);
	Output* Output 			= &(Model->Output);
	Breakpoint* Breakpoint 	= &(Model->Breakpoint);
	



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
	int iSub, iSub2, iPhaseAttr, iPhase;
	int j;
	int size, size2, subSize;

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
				} else if (  TOKEN("maxTime") ) {
					Numerics->maxTime = atof(strValue);
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
				} else if  (  TOKEN("dtAlphaCorr") ) {
					Numerics->dtAlphaCorrIni = atof(strValue);
					Numerics->dtAlphaCorr = Numerics->dtAlphaCorrIni;
				} else if  (  TOKEN("CFL_fac_Stokes") ) {
					Numerics->CFL_fac_Stokes = atof(strValue);
				} else if  (  TOKEN("CFL_fac_Thermal") ) {
					Numerics->CFL_fac_Thermal = atof(strValue);
				} else if  (  TOKEN("CFL_fac_Darcy") ) {
					Numerics->CFL_fac_Darcy = atof(strValue);
				} else if  (  TOKEN("use_dtMaxwellLimit") ) {
					Numerics->use_dtMaxwellLimit = VALUE("true");
				} else if  (  TOKEN("etaMin") ) {
					Numerics->etaMin = atof(strValue);
				} else if  (  TOKEN("etaMax") ) {
					Numerics->etaMax = atof(strValue);
				} else if  (  TOKEN("phiMin") ) {
					Numerics->phiMin = atof(strValue);
				} else if  (  TOKEN("phiMax") ) {
					Numerics->phiMax = atof(strValue);
				} else if  (  TOKEN("phiCrit") ) {
					Numerics->phiCrit = atof(strValue);
				} else if  (  TOKEN("dtMin") ) {
					Numerics->dtMin = atof(strValue);
				} else if  (  TOKEN("dtMax") ) {
					Numerics->dtMax = atof(strValue);
				} else if  (  TOKEN("dtVep") ) {
					Numerics->dtVep = atof(strValue);
				} else if  (  TOKEN("dtMaxwellFac_EP_ov_E") ) {
					Numerics->dtMaxwellFac_EP_ov_E = atof(strValue);
				} else if  (  TOKEN("dtMaxwellFac_VP_ov_E") ) {
					Numerics->dtMaxwellFac_VP_ov_E = atof(strValue);
				} else if  (  TOKEN("dtMaxwellFac_VP_ov_EP") ) {
					Numerics->dtMaxwellFac_VP_ov_EP = atof(strValue);
				} else if  (  TOKEN("dt_stressFac") ) {
					Numerics->dt_stressFac = atof(strValue);
				} else if  (  TOKEN("dt_plasticFac") ) {
					Numerics->dt_plasticFac = atof(strValue);	
				} else if  (  TOKEN("deltaSigmaMin") ) {
					Numerics->deltaSigmaMin = atof(strValue);
				} else if  (  TOKEN("yieldComputationType") ) {
					Numerics->yieldComputationType = atoi(strValue);	
				} else if  (  TOKEN("invariantComputationType") ) {
					Numerics->invariantComputationType = atoi(strValue);	
				} else if (   TOKEN("stickyAirSwitchingDepth") ) {
					Numerics->stickyAirSwitchingDepth = atof(strValue);
				} else if (   TOKEN("stickyAirSwitchPhaseTo") ) {
					Numerics->stickyAirSwitchPhaseTo = atoi(strValue);
				} else if (   TOKEN("stickyAirSwitchPassiveTo") ) {
					Numerics->stickyAirSwitchPassiveTo = atoi(strValue);
				} else if (   TOKEN("stickyAirTimeSwitchPassive") ) {
					Numerics->stickyAirTimeSwitchPassive = atof(strValue);
				} else if (   TOKEN("stressSubGridDiffFac") ) {
					Numerics->stressSubGridDiffFac = atof(strValue);
				} else if (   TOKEN("dtIni") ) {
					Numerics->dtIni = atof(strValue);
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
				} else if  (  TOKEN("Pback") ) {
					Physics->Pback = atof(strValue);
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

				} else if  (  TOKEN("passiveDx") ) {
					Particles->passiveDx = atof(strValue);
				} else if  (  TOKEN("passiveDy") ) {
					Particles->passiveDy = atof(strValue);

				} else if  (  TOKEN("passiveGeom") ) {
					if 		  ( VALUE("Grid")) {
						Particles->passiveGeom = PartPassive_Grid;
					} else if 		  ( VALUE("Grid_w_Layers")) {
						Particles->passiveGeom = PartPassive_Grid_w_Layers;
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



		else if (TOKEN("BC")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;
				if (TOKEN("Stokes")) {
					i++; // Move to the first token, which is the object
					size2 = t[i].size; // number of elements in the token
					i++; // Move to the first key

					for (iSub2=0; iSub2<size2; iSub2++) {
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("backStrainRate") ) {
							BCStokes->backStrainRate = atof(strValue);
						} else if (  TOKEN("specialPhase") ) {
							BCStokes->specialPhase = atoi(strValue);
						} else if (  TOKEN("refValue") ) {
							BCStokes->refValue = atof(strValue);
						} else if (  TOKEN("DeltaL") ) {
							BCStokes->DeltaL = atof(strValue);

						} else if (  TOKEN("Sandbox_TopSeg00") ){
							BCStokes->Sandbox_TopSeg00 = atof(strValue);
						} else if (  TOKEN("Sandbox_TopSeg01") ){
							BCStokes->Sandbox_TopSeg01 = atof(strValue);
						} else if (  TOKEN("Sandbox_NoSlipWall") ){
							BCStokes->Sandbox_NoSlipWall = VALUE("true");
						} else if (  TOKEN("Corner_SubductionAngle") ){
							BCStokes->Corner_SubductionAngle = atof(strValue);
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
							} else if ( VALUE("WindTunnel")) {
								BCStokes->SetupType = Stokes_WindTunnel;

							} else {
								printf("\n   ### ERROR ###: Unexpected BCStokes.type: %.*s\n\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
								Stop = true;
							}
						} else {
							printf("Unexpected key in BCStokes: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
							Stop = true;
						}

						i+=2;
					}
					//i-=2;
				}
				else if (TOKEN("Thermal")) {
					i++; // Move to the first token, which is the object
					size2 = t[i].size; // number of elements in the token
					i++; // Move to the first key

					for (iSub2=0; iSub2<size2; iSub2++) {
						strValue = JSON_STRING+t[i+1].start;

						if 		  (  TOKEN("TT") ) {
							BCThermal->TT = atof(strValue);
						} else if (  TOKEN("TB") ) {
							BCThermal->TB = atof(strValue);
						} else if (  TOKEN("refValue") ) {
							BCThermal->refValue = atof(strValue);
						} else if (  TOKEN("DeltaL") ) {
							BCThermal->DeltaL = atof(strValue);
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
					//i-=2;
				} else {
					printf("Unexpected key in BC: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}
			}
		}



		else if (TOKEN("IC")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key

			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;

				if (TOKEN("Thermal")) {
					//#if (HEAT)
					i++; // Move to the first token, which is the object
					size2 = t[i].size; // number of elements in the token
					i++; // Move to the first key
					strValue = JSON_STRING+t[i+1].start;
					if 		  ( VALUE("HSC")) {
						ICThermal->SetupType = IC_HSC;
						printf("ICThermal Type = OK\n");
					} else if ( VALUE("Gaussian")){
						ICThermal->SetupType = IC_Gaussian;
					} else {
						printf("Unexpected type in ICThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					}
					i+=2;
					size2 -= 1;
					if (ICThermal->SetupType == IC_HSC) {
						for (iSub2=0; iSub2<size2; iSub2++) {
							strValue = JSON_STRING+t[i+1].start;

							if (  TOKEN("noise") ) {
								ICThermal->data[0] = atof(strValue);
							} else if (  TOKEN("Tm") ) {
								ICThermal->data[1] = atof(strValue);
							} else if (  TOKEN("age") ) {
								ICThermal->data[2] = atof(strValue);
							} else {
								printf("Unexpected key in ICThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
								Stop = true;
							}
							i+=2;
						}

					} else if (ICThermal->SetupType == IC_Gaussian) {
						for (iSub2=0; iSub2<size2; iSub2++) {
							strValue = JSON_STRING+t[i+1].start;

							if (  TOKEN("noise") ) {
								ICThermal->data[0] = atof(strValue);
							} else if (  TOKEN("background") ) {
								ICThermal->data[1] = atof(strValue);
							} else if (  TOKEN("Amp") ) {
								ICThermal->data[2] = atof(strValue);
							} else if (  TOKEN("xc") ) {
								ICThermal->data[3] = atof(strValue);
							} else if (  TOKEN("yc") ) {
								ICThermal->data[4] = atof(strValue);
							} else if (  TOKEN("wx") ) {
								ICThermal->data[5] = atof(strValue);
							} else if (  TOKEN("wy") ) {
								ICThermal->data[6] = atof(strValue);
							} else {
								printf("Unexpected key in ICThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
								Stop = true;
							}
							i+=2;
						}

					}
					//i-=2;
					//#endif
				}
				else if (TOKEN("Darcy")) {
					//#if (DARCY)
					i++; // Move to the first token, which is the object
					size2 = t[i].size; // number of elements in the token
					i++; // Move to the first key
					printf("koko\n");
					strValue = JSON_STRING+t[i+1].start;
					if 		  ( VALUE("HSC")) {
						printf("error: Initial conditional for Darcy is set to HSC (half space cooling), which cannot be applied\n");
						exit(0);
					} else if ( VALUE("Gaussian")){
						ICDarcy->SetupType = IC_Gaussian;
					} else {
						printf("Unexpected type in ICThermal: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					}
					i+=2;
					size2 -= 1;
					if (ICDarcy->SetupType == IC_Gaussian) {
						for (iSub2=0; iSub2<size2; iSub2++) {
							strValue = JSON_STRING+t[i+1].start;

							if (  TOKEN("noise") ) {
								ICDarcy->data[0] = atof(strValue);
							} else if (  TOKEN("background") ) {
								ICDarcy->data[1] = atof(strValue);
							} else if (  TOKEN("Amp") ) {
								ICDarcy->data[2] = atof(strValue);
							} else if (  TOKEN("xc") ) {
								ICDarcy->data[3] = atof(strValue);
							} else if (  TOKEN("yc") ) {
								ICDarcy->data[4] = atof(strValue);
							} else if (  TOKEN("wx") ) {
								ICDarcy->data[5] = atof(strValue);
							} else if (  TOKEN("wy") ) {
								ICDarcy->data[6] = atof(strValue);
							} else {
								printf("Unexpected key in ICDarcy: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
								Stop = true;
							}
							i+=2;
						}

					}
					//#endif
					//i-=2;
				}
			} // for
		}















		else if (TOKEN("MatProps")) {
			i++; // Move to the first token, which is the object
			MatProps->nPhase = t[i].size; // number of elements in the token
			//printf("nPhase = %i\n",t[i].size);
			i++; // Move to the first key


			for (iSub=0; iSub<MatProps->nPhase; iSub++) {
				strToken = JSON_STRING+t[i].start;
				iPhase = atoi(strToken);

				if (iPhase>NB_PHASE_MAX) {
					printf("The number of phases is greater than NB_PHASE_MAX. Reduce the number of phases or increase NB_PHASE_MAX (%i).\n", NB_PHASE_MAX);
					exit(0);
				}

				//printf("PhaseNumber = %i\n",PhaseNumber);
				//printf("Phase token: %.*s\n", t[i].end-t[i].start, strToken);


				i++;
				subSize = t[i].size; // number of elements in the token
				i++;
				for (iPhaseAttr=0; iPhaseAttr<subSize; iPhaseAttr++) {
					strToken = JSON_STRING+t[i].start;
					strValue = JSON_STRING+t[i+1].start;

					if 		  (  TOKEN("rho0") ) {
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
					} else if (  TOKEN("phiIni") ) {
						MatProps->phiIni[iPhase] = atof(strValue);

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

					} else if ( 	TOKEN("use_dtMaxwellLimit") ) {
						MatProps->use_dtMaxwellLimit[iPhase] = VALUE("true");

					} else if ( 	TOKEN("staticPfFac") ) {
						MatProps->staticPfFac[iPhase] = atof(strValue);
					} else if ( 	TOKEN("staticPfFacWeakFac") ) {
						MatProps->staticPfFacWeakFac[iPhase] = atof(strValue);
					} else if ( 	TOKEN("frictionAngleWeakFac") ) {
						MatProps->frictionAngleWeakFac[iPhase] = atof(strValue);
					} else if ( 	TOKEN("cohesionWeakFac") ) {
						MatProps->cohesionWeakFac[iPhase] = atof(strValue);
					} else if ( 	TOKEN("strainWeakStart") ) {
						MatProps->strainWeakStart[iPhase] = atof(strValue);
					} else if ( 	TOKEN("strainWeakEnd") ) {
						MatProps->strainWeakEnd[iPhase] = atof(strValue);
					} else if  (  TOKEN("vDiff") ) {

						i++;
						size2 = t[i].size; // number of elements in the token
						i++;
						for (iSub2=0; iSub2<size2; iSub2++) {
							strToken = JSON_STRING+t[i].start;
							strValue = JSON_STRING+t[i+1].start;
							//["flowLaw","A","E","V","tensorCorrection","MPa","d0","p","C_OH_0","r","isActive"]
							if 		(  TOKEN("flowLaw") ) {
								// Place holder
							} else if 	(  TOKEN("A") ) {
								// Place holder
							} else if 	(  TOKEN("B") ) {
								MatProps->vDiff[iPhase].B = atof(strValue);
							} else if 	(  TOKEN("E") ) {
								MatProps->vDiff[iPhase].E = atof(strValue);
							} else if 	(  TOKEN("V") ) {
								MatProps->vDiff[iPhase].V = atof(strValue);
							} else if 	(  TOKEN("d0") ) {
								// Place holder
							} else if 	(  TOKEN("p") ) {
								// Place holder
							} else if 	(  TOKEN("C_OH_0") ) {
								// Place holder
							} else if 	(  TOKEN("r") ) {
								// Place holder
							} else if 	(  TOKEN("isActive") ) {
								MatProps->vDiff[iPhase].isActive = VALUE("true");
							} else if 	(  TOKEN("tensorCorrection") ) {
								// Place holder
							}
							i += 2;
						}
						i-=2;
					} else if  (  TOKEN("vDisl") ) {
						i++;
						size2 = t[i].size; // number of elements in the token
						i++;
						for (iSub2=0; iSub2<size2; iSub2++) {
							strToken = JSON_STRING+t[i].start;
							strValue = JSON_STRING+t[i+1].start;
							//["flowLaw","A","E","V","tensorCorrection","MPa","d0","p","C_OH_0","r","isActive"]
							if 	(  TOKEN("flowLaw") ) {
								// Place holder
							} else if 	(  TOKEN("A") ) {
								// Place holder
							} else if 	(  TOKEN("B") ) {
								MatProps->vDisl[iPhase].B = atof(strValue);
								//printf("Assigned vDisl of iPhase = %i, = %.2e\n",iPhase, MatProps->vDisl[iPhase].B);
							} else if 	(  TOKEN("E") ) {
								MatProps->vDisl[iPhase].E = atof(strValue);
							} else if 	(  TOKEN("V") ) {
								MatProps->vDisl[iPhase].V = atof(strValue);
							} else if 	(  TOKEN("n") ) {
								MatProps->vDisl[iPhase].n = atof(strValue);
							} else if 	(  TOKEN("C_OH_0") ) {
								// Place holder
							} else if 	(  TOKEN("r") ) {
								// Place holder
							} else if 	(  TOKEN("isActive") ) {
								MatProps->vDisl[iPhase].isActive = VALUE("true");
							} else if 	(  TOKEN("tensorCorrection") ) {
								// Place holder

							}
							i += 2;
						}
						i-=2;

					} else if  (  TOKEN("vPei") ) {

						i++;
						size2 = t[i].size; // number of elements in the token
						i++;
						for (iSub2=0; iSub2<size2; iSub2++) {
							strToken = JSON_STRING+t[i].start;
							strValue = JSON_STRING+t[i+1].start;
							//["flowLaw","A","E","V","tensorCorrection","MPa","d0","p","C_OH_0","r","isActive"]
							if 	(  TOKEN("flowLaw") ) {
								// Place holder
							}else if  	(  TOKEN("B") ) {
								MatProps->vPei[iPhase].B = atof(strValue);
							} else if 	(  TOKEN("E") ) {
								MatProps->vPei[iPhase].E = atof(strValue);
							} else if 	(  TOKEN("V") ) {
								MatProps->vPei[iPhase].V = atof(strValue);
							} else if 	(  TOKEN("tau") ) {
								MatProps->vPei[iPhase].tau = atof(strValue);
							} else if 	(  TOKEN("gamma") ) {
								MatProps->vPei[iPhase].gamma = atof(strValue);
							} else if 	(  TOKEN("q") ) {
								MatProps->vPei[iPhase].q = atof(strValue);
							} else if 	(  TOKEN("isActive") ) {
								MatProps->vPei[iPhase].isActive = VALUE("true");
							}
							i += 2;
						}
						i-=2;
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


			/*
			int charLength = t[i+1].end-t[i+1].start+1;
			Output->ModelDescription = (char*) malloc(charLength * sizeof(char));

			strValue = JSON_STRING+t[i+1].start;
			strncpy(Output->ModelDescription, strValue, charLength-1);
			Output->ModelDescription[charLength-1] = '\0';
			*/
			i++; // Move to the first token, which is the object
			i++; // Move to the first key
		}


		else if (TOKEN("Output")) {
			i++; // Move to the first token, which is the object
			size = t[i].size; // number of elements in the token
			i++; // Move to the first key
			Output->nTypes = 0;
			Output->nPartTypes = 0;
			for (iSub=0; iSub<size; iSub++) {
				strValue = JSON_STRING+t[i+1].start;
				if  (  TOKEN("folder") ) {
					if (t[i+1].end-t[i+1].start>MAX_STRING_LENGTH) {
						printf("the Visu.outputFolder string is too long, maximum authorized: %i. Please change your folder or increase the value of the macro MAX_STRING_LENGTH", MAX_STRING_LENGTH);
					}

					strncpy(Output->outputFolder, strValue, t[i+1].end-t[i+1].start);
					// Check for "/" at the end of the Folder name
					Output->outputFolder[t[i+1].end-t[i+1].start] = '\0';


					printf("Data Output folder: %s\n",Output->outputFolder);
				} else if  (  TOKEN("breakpointFolder") ) {
					if (t[i+1].end-t[i+1].start>MAX_STRING_LENGTH) {
						printf("the Visu.breakpointFolder string is too long, maximum authorized: %i. Please change your folder or increase the value of the macro MAX_STRING_LENGTH", MAX_STRING_LENGTH);
					}

					strncpy(Breakpoint->breakpointFolder, strValue, t[i+1].end-t[i+1].start);
					// Check for "/" at the end of the Folder name
					Breakpoint->breakpointFolder[t[i+1].end-t[i+1].start] = '\0';


					printf("Breakpoint folder: %s\n",Breakpoint->breakpointFolder);
				} else if 	(  TOKEN("Vx") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Vx;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("Vy") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Vy;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("P") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_P;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("Pf") ) {
#if (DARCY)
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Pf;
						Output->nTypes++;
					}
#endif
				} else if  	(  TOKEN("Pc") ) {
#if (DARCY)
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Pc;
						Output->nTypes++;
					}
#endif
				} else if  	(  TOKEN("eta") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Viscosity;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("porosity") ) {
#if (DARCY)
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Porosity;
						Output->nTypes++;
					}
#endif
				} else if  	(  TOKEN("Z") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Z;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("G") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_G;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("khi") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Khi;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_xx0") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Sxx0;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_xy0") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Sxy0;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_xx") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Sxx;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_xy") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Sxy;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_xy_node") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Sxy_Node;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("sigma_II") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_SII;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("strainRate") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_StrainRate;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("strain") ) {
#if (STORE_PLASTIC_STRAIN)
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Strain;
						Output->nTypes++;
					}
#endif
				} else if  	(  TOKEN("temperature") ) {
#if (HEAT)
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Temperature;
						Output->nTypes++;
					}
#endif
				} else if  	(  TOKEN("phase") ) {
					if (VALUE("true")) {
						Output->type[Output->nTypes] = Out_Phase;
						Output->nTypes++;
					}
				} else if  	(  TOKEN("frequency") ) {
					Output->frequency = atoi(strValue);
				} else if  	(  TOKEN("timeFrequency") ) {
					Output->timeFrequency = atof(strValue);
					if (Output->timeFrequency>0.0) {
						Output->useTimeFrequency = true;
					} else {
						Output->useTimeFrequency = false;
					}
				} else if  	(  TOKEN("breakpointFrequency") ) {
					Breakpoint->frequency = atoi(strValue);
				} else if  	(  TOKEN("saveFirstStep") ) {
					Output->saveFirstStep = VALUE("true");





				} else if  	(  TOKEN("particles_pos") ) {
					if (VALUE("true")) {
						Output->partType[Output->nPartTypes] = OutPart_x;
						Output->nPartTypes++;
						Output->partType[Output->nPartTypes] = OutPart_y;
						Output->nPartTypes++;

					}
				} else if  	(  TOKEN("particles_posIni") ) {
#if (STORE_PARTICLE_POS_INI)
					if (VALUE("true")) {
						Output->partType[Output->nPartTypes] = OutPart_xIni;
						Output->nPartTypes++;
						Output->partType[Output->nPartTypes] = OutPart_yIni;
						Output->nPartTypes++;
					}
#endif
				} else if  	(  TOKEN("particles_phase") ) {
					if (VALUE("true")) {
						Output->partType[Output->nPartTypes] = OutPart_Phase;
						Output->nPartTypes++;
					}
				} else if  	(  TOKEN("particles_passive") ) {
					if (VALUE("true")) {
						Output->partType[Output->nPartTypes] = OutPart_Passive;
						Output->nPartTypes++;
					}
				} else if  	(  TOKEN("particles_T") ) {
					if (VALUE("true")) {
#if (HEAT)
						Output->partType[Output->nPartTypes] = OutPart_T;
						Output->nPartTypes++;
#endif
					}
				} else if  	(  TOKEN("particles_stress") ) {
					if (VALUE("true")) {
						Output->partType[Output->nPartTypes] = OutPart_Sxx0;
						Output->nPartTypes++;
						Output->partType[Output->nPartTypes] = OutPart_Sxy0;
						Output->nPartTypes++;
#if (DARCY)
						Output->partType[Output->nPartTypes] = OutPart_DeltaP0;
						Output->nPartTypes++;
#endif
					}
				} else if  	(  TOKEN("particles_phi") ) {
					if (VALUE("true")) {
#if (DARCY)
						Output->partType[Output->nPartTypes] = OutPart_Phi;
						Output->nPartTypes++;
#endif
					}
				} else if  	(  TOKEN("particles_strain") ) {
					if (VALUE("true")) {
#if (STORE_PLASTIC_STRAIN)
						Output->partType[Output->nPartTypes] = OutPart_Strain;
						Output->nPartTypes++;
#endif
					}
				} else if  	(  TOKEN("particles_timeLastPlastic") ) {
					if (VALUE("true")) {
#if (STORE_TIME_LAST_PLASTIC)
						Output->partType[Output->nPartTypes] = OutPart_TimeLastPlastic;
						Output->nPartTypes++;
#endif
					}
				} else {
					printf("Unexpected key in Output: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
					Stop = true;
				}

				i+=2;
			}
			printf("nTypes = %i\n",Output->nTypes);
		}




		else {
			printf("Unexpected key: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
			Stop = true;
			i++;

		}


	}



	if (Stop) {
		printf("error: unexpected keys were found. Some variables could not be assigned values");
		exit(0);
	}



	free(JSON_STRING);


}


#if (VISU)

void Input_readVisu(Model* Model)
{
	Visu* Visu 				= &(Model->Visu);
	Input* Input 			= &(Model->Input);


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
				//printf("This key is: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
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
				} else if  (  TOKEN("renderFrequency") ) {
					Visu->renderFrequency = atoi(strValue);
				} else if  (  TOKEN("renderTimeFrequency") ) {
					Visu->renderTimeFrequency = atof(strValue);
				} else if  (  TOKEN("closeAtTheEndOfSimulation") ) {
					Visu->closeAtTheEndOfSimulation = VALUE("true");
				} else if  (  TOKEN("outputFolder") ) {
					if (t[i+1].end-t[i+1].start>MAX_STRING_LENGTH) {
						printf("the Visu.outputFolder string is too long, maximum authorized: %i. Please change your folder or increase the value of the macro MAX_STRING_LENGTH", MAX_STRING_LENGTH);
					}

					strncpy(Visu->outputFolder, strValue, t[i+1].end-t[i+1].start);
					Visu->outputFolder[t[i+1].end-t[i+1].start] = '\0';

					printf("%s\n",Visu->outputFolder);


				} else if  (  TOKEN("shaderFolder") ) {
					if (t[i+1].end-t[i+1].start>MAX_STRING_LENGTH) {
						printf("the Visu.outputFolder string is too long, maximum authorized: %i. Please change your folder or increase the value of the macro MAX_STRING_LENGTH", MAX_STRING_LENGTH);
					}

					strncpy(Visu->shaderFolder, strValue, t[i+1].end-t[i+1].start);
					Visu->shaderFolder[t[i+1].end-t[i+1].start] = '\0';

					printf("%s\n",Visu->shaderFolder);
					//memset(Visu->outputFolder, '\0', t[i+1].end-t[i+1].start);

				} else if  (  TOKEN("retinaScale") ) {
					Visu->retinaScale = atoi(strValue);








				} else if  (  TOKEN("glyphType") ) {

					if 		  ( VALUE("StokesVelocity")) {
						Visu->glyphType = VisuGlyphType_StokesVelocity;
					} else if ( VALUE("DarcyGradient")) {
						Visu->glyphType = VisuGlyphType_DarcyGradient;
					} else if ( VALUE("DeviatoricStressTensor")) {
						Visu->glyphType = VisuGlyphType_DeviatoricStressTensor;
					} else {
						printf("Unexpected Visu.glyphType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}

				} else if  (  TOKEN("filter") ) {

					if 		  ( VALUE("Linear")) {
						Visu->filter = VisuFilterType_Linear;
					} else if ( VALUE("Nearest")) {
						Visu->filter = VisuFilterType_Nearest;
					} else {
						printf("Unexpected Visu.glyphType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}





				} else if  (  TOKEN("glyphMeshType") ) {
					if 		  ( VALUE("Triangle")) {
						Visu->glyphMeshType = VisuGlyphMeshType_Triangle;
					} else if ( VALUE("ThinArrow")) {
						Visu->glyphMeshType = VisuGlyphMeshType_ThinArrow;
					} else if ( VALUE("ThickArrow")) {
						Visu->glyphMeshType = VisuGlyphMeshType_ThickArrow;
					} else if ( VALUE("TensorCross")) {
						Visu->glyphMeshType = VisuGlyphMeshType_TensorCross;
					} else {
						printf("Unexpected Visu.particleType: %.*s\n", t[i+1].end-t[i+1].start, JSON_STRING + t[i+1].start);
						exit(0);
					}


				} else if  (  TOKEN("typeParticles") ) {
					if 		  ( VALUE("PartPhase")) {
						Visu->typeParticles = VisuType_PartPhase;
					} else if ( VALUE("PartTemp")) {
						Visu->typeParticles = VisuType_PartTemp;
					} else if ( VALUE("PartSigma_xx")) {
						Visu->typeParticles = VisuType_PartSigma_xx;
					} else if ( VALUE("PartSigma_xy")) {
						Visu->typeParticles = VisuType_PartSigma_xy;
					} else if ( VALUE("PartStrain")) {
						Visu->typeParticles = VisuType_PartStrain;
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

							if 			(  TOKEN("A0number") ) {
								thisType = atoi(strValue);
							} else if	(  TOKEN("center") ) {
								Visu->colorMap[thisType].center = atof(strValue);
								//printf("thisType = %i\n",thisType);
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
							} else if	(  TOKEN("alphaAbsThreshold") ) {
								Visu->colorMap[thisType].alphaAbsThreshold = atof(strValue);
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





void Input_assignPhaseToParticles(Model* Model)
{
Input* Input 			= &(Model->Input);
Particles* Particles 	= &(Model->Particles);
Grid* Grid 				= &(Model->Grid);
Char* Char 				= &(Model->Char);
	


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
						} else if (  TOKEN("slope") ) {
							sine.slope = atof(strValue);
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

	compute xmin,xmax, ymin, ymax;
	compute xminL, xmaxL, yminL, ymaxL; // in the coordinate of the line (where xmin or ymin = 0)
	if (Line->definedFor == 1) {
		xmin = Line->min;
		xmax = Line->max;
		ymin = 0.0;
		ymax = 0.0;
		xminL = 0.0;
		xmaxL = Line->max-Line->min;
		if (a>0) {
			if (Line->condition == 1) {
				ymin = a*xminL + b;
				ymax = Grid->ymax;
			}
			else if (Line->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xmaxL + b;
			}
		} else {
			if (Line->condition == 1) {
				ymin = a*xmaxL + b;
				ymax = Grid->ymax;
			}
			else if (Line->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xminL + b;
			}
		}
	} else if (Line->definedFor == 0) {
		ymin = Line->min;
		ymax = Line->max;
		xmin = 0.0;
		xmax = 0.0;
		yminL = 0.0;
		ymaxL = Line->max-Line->min;
		if (a>0) {
			if (Line->condition == 1) {
				xmin = a*yminL + b;
				xmax = Grid->xmax;
			}
			else if (Line->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*ymaxL + b;
			}

		} else {
			if (Line->condition == 1) {
				xmin = a*ymaxL + b;
				xmax = Grid->xmax;
			}
			else if (Line->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*yminL + b;
			}
		}
	} else {
		printf("error: in Input assignLine Line->definedFor couldn't be understand.\n");
		exit(0);
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
					x -= Line->min;
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
					y -= Line->min;
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

	compute a = Sine->slope;
	compute b = Sine->base;
	compute Amp = Sine->amplitude;
	compute xmin,xmax, ymin, ymax;
	compute xminL, xmaxL, yminL, ymaxL;


	if (Sine->definedFor == 1) {
		xmin = Sine->min;
		xmax = Sine->max;
		ymin = 0.0;
		ymax = 0.0;
		xminL = 0.0;
		xmaxL = Sine->max-Sine->min;
		if (a>0) {
			if (Sine->condition == 1) {
				ymin = a*xminL + b - Amp;
				ymax = Grid->ymax;
			}
			else if (Sine->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xmaxL + b + Amp;
			}

		} else {
			if (Sine->condition == 1) {
				ymin = a*xmaxL + b - Amp;
				ymax = Grid->ymax;
			}
			else if (Sine->condition == 0) {
				ymin = Grid->ymin;
				ymax = a*xminL + b + Amp;
			}

		}
	} else if (Sine->definedFor == 0) {
		ymin = Sine->min;
		ymax = Sine->max;
		xmin = 0.0;
		xmax = 0.0;
		yminL = 0.0;
		ymaxL = Sine->max-Sine->min;
		if (a>0) {
			if (Sine->condition == 1) {
				xmin = a*yminL + b - Amp;
				xmax = Grid->xmax;
			}
			else if (Sine->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*ymaxL + b + Amp;
			}

		} else {
			if (Sine->condition == 1) {
				xmin = a*ymaxL + b - Amp;
				xmax = Grid->xmax;
			}
			else if (Sine->condition == 0) {
				xmin = Grid->xmin;
				xmax = a*yminL + b + Amp;
			}
		}
	} else {
		printf("error: in Input assignSine Sine->definedFor couldn't be understand.\n");
		exit(0);
	}


	/*
	if (Sine->definedFor == 1) {
		xmin = Sine->min;
		xmax = Sine->max;
		xminL = 0.0;
		xmaxL = Sine->max-Sine->min;
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
		yminL = 0.0;
		ymaxL = Sine->max-Sine->min;
		if (Sine->condition == 1) {
			xmin = Sine->base-Sine->amplitude;
			xmax = Grid->xmax;
		}
		else if (Sine->condition == 0) {
			xmax = Sine->base+Sine->amplitude;
			xmin = Grid->xmin;
		}
	}
	 */


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
					x = (thisParticle->x - Sine->min);///(Grid->xmax-Grid->xmin);
					y = thisParticle->y;
					if ( Sine->condition == 1 ) { // >
						if ( y > a*x + b + Amp*sin(1.0/Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( y < a*x + b + Amp*sin(1.0/Sine->wavelength*x*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					}
				} else if (Sine->definedFor == 0) {
					x = thisParticle->x;///(Grid->xmax-Grid->xmin);
					y = thisParticle->y-Sine->min;
					if ( Sine->condition == 1 ) { // >
						if ( x >a*y + b + Amp*sin(1.0/Sine->wavelength*y*2*PI+ Sine->wavephase)) {
							thisParticle->phase = Sine->phase;
						}
					} else if ( Sine->condition == 0 ) { // <
						if ( x < a*y + b + Amp*sin(1.0/Sine->wavelength*y*2*PI+ Sine->wavephase)) {
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
	//if (ix!=0)
	//	ix = ix-1;
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
	//if (iy!=0)
	//	iy = iy-1;
	indexLimits[2] = iy;
	while (iy < Grid->nyS-1 && Grid->Y[iy]<coordLimits[3]) {
		iy++;
	}
	//if (iy==Grid->nyS)
	//	iy = iy-1;
	indexLimits[3] = iy+1;



}


#if (EXTRA_PART_FIELD)

void Input_setFieldsOnParticles(Model* Model) {
	
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Particles* Particles 	= &(Model->Particles);


	compute* buffer = (compute*) malloc (Grid->nxEC *Grid->nyEC* sizeof(compute));
	//compute* buffer = (compute*) malloc (Grid->nxEC *Grid->nyEC* sizeof(compute));
	//unsigned char buffer[10];
	//compute* buffer;
	FILE *ptr;

	ptr = fopen("/Users/abauville/Work/ProjectWithMarcel/Input/test.bin","rb");  // r for read, b for binary
	fread(buffer,sizeof(compute),Grid->nECTot,ptr); // read 10 bytes to our buffer
	fclose(ptr);

	int ix, iy, iNode;
	/*
	
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			Physics->extraField[iCell] = buffer[iCell];
			//printf("%.2e  ", buffer[iCell]);
		}
		printf("\n");
	}
	*/

	SingleParticle* thisParticle;

	
	compute locX, locY;
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			
			while (thisParticle!=NULL) {
				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);

				thisParticle->extraField = Interp_ECVal_Cell2Particle_Local(buffer, ix, iy, Grid->nxEC, locX, locY);
				thisParticle = thisParticle->next;

			}
			
		}
	}
	
	
	/*
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix+iy*Grid->nxEC;
			printf("%.2f  ", Physics->extraField[iCell]);
		}
		printf("\n");
	}
	*/

	free(buffer);
	
}
#endif