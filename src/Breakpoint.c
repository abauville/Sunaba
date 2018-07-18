/*
 * Breapoint.c
 *
 *  Created on: 18 Jul 2018
 *      Author: abauville
 */



#include "stokes.h"
#include "MiscLibraries/jsmn.h"
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



// Works only for UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stddef.h>



void Breakpoint_writeData (Model* Model){
	Output* Output 			= &(Model->Output);
	OutType outputTypeCopy[20];
	OutPartType outputPartTypeCopy[13];
	int nTypesCopy = Output->nTypes;
	int nPartTypesCopy = Output->nPartTypes;
	char outputFolderCopy[MAX_STRING_LENGTH];


	// Make a copy of Output->type
	memcpy(outputTypeCopy, Output->type, 20 * sizeof(OutType));
	memcpy(outputPartTypeCopy, Output->partType, 13 * sizeof(OutPartType));
	memcpy(outputFolderCopy, Output->outputFolder, MAX_STRING_LENGTH * sizeof(char));
	memcpy(Output->outputFolder, Output->breakpointFolder, MAX_STRING_LENGTH * sizeof(char));

	int i;

	// Modify Output->type to include stuff useful for the restart
	Output->type[ 0] = Out_Vx;
	Output->type[ 1] = Out_Vy;
	Output->type[ 2] = Out_P;
	Output->type[ 3] = Out_Z;
	Output->nTypes = 4;

	i = 0;
	Output->partType[ i] = OutPart_x; i++;
	Output->partType[ i] = OutPart_y; i++;
	Output->partType[ i] = OutPart_Phase; i++;
	Output->partType[ i] = OutPart_Passive; i++;
	Output->partType[ i] = OutPart_Sxx0; i++;
	Output->partType[ i] = OutPart_Sxy0; i++;
	
#if (STORE_PARTICLE_POS_INI)
	Output->partType[ i] = OutPart_xIni; i++;
	Output->partType[ i] = OutPart_yIni; i++;
#endif
#if (STORE_PLASTIC_STRAIN)
	Output->partType[ i] = OutPart_Strain; i++;
#endif
#if (STORE_TIME_LAST_PLASTIC)
	Output->partType[ i] = OutPart_TimeLastPlastic; i++;
#endif
#if (HEAT)
	Output->partType[ i] = OutPart_T; i++;
#endif
#if (DARCY)
	Output->partType[ i] = OutPart_DeltaP0; i++;
	Output->partType[i] = OutPart_Phi; i++;
#endif
	Output->nPartTypes = i;

	// Call  Output_data and Output_particles
	Output_modelState(Model);
	Output_data(Model);
	Output_particles(Model);

	// Restore Output->type from the copy
	memcpy(Output->type, outputTypeCopy, 20 * sizeof(OutType));
	memcpy(Output->partType, outputPartTypeCopy, 13 * sizeof(OutPartType));
	memcpy(Output->outputFolder, outputFolderCopy, MAX_STRING_LENGTH * sizeof(char));

	Output->nTypes = nTypesCopy;
	Output->nPartTypes = nPartTypesCopy;

}


void Breakpoint_readData(Model* Model) {
	Physics* Physics 		= &(Model->Physics);
	Output* Output 			= &(Model->Output);
	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];
	sprintf(Folder_thistStep, "%sOut_%05i/", Output->breakpointFolder,Output->counter);

	
	// Matrix Data
	sprintf(fname,"%s%s.bin",Folder_thistStep, "P");
	Breakpoint_fillMatrixFromDataFile(Physics->P , fname, Model);
	sprintf(fname,"%s%s.bin",Folder_thistStep, "Vx");
	Breakpoint_fillMatrixFromDataFile(Physics->Vx, fname, Model);
	sprintf(fname,"%s%s.bin",Folder_thistStep, "Vy");
	Breakpoint_fillMatrixFromDataFile(Physics->Vy, fname, Model);
	sprintf(fname,"%s%s.bin",Folder_thistStep, "Z");
	Breakpoint_fillMatrixFromDataFile(Physics->Z , fname, Model);

}



void Breakpoint_fillMatrixFromDataFile(compute* MatrixToFill, char* fname, Model* Model) {
	FILE * fptr;
	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	// Open bin file
	fptr = fopen(fname,"rb");
	
	// Read the header
	fread(nxy , sizeof(int), 2, fptr);
	fread(&xmin, sizeof(double), 1, fptr);
	fread(&xmax, sizeof(double), 1, fptr);
	fread(&ymin, sizeof(double), 1, fptr);
	fread(&ymax, sizeof(double), 1, fptr);
	fread(&Char_quantity, sizeof(double), 1, fptr);
	
	// Collect the data and fill the matrix
	fwrite(MatrixToFill, sizeof(double), nxy[0]*nxy[1], fptr);

	// Close the file
	fclose(fptr);
}


void Breakpoint_readModelState(Model* Model) {

	// Declare structures
	// =================================
	Grid* Grid 				= &(Model->Grid);
	Numerics* Numerics 		= &(Model->Numerics);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);
	Input* Input 			= &(Model->Input);
	Output* Output 			= &(Model->Output);
	
	// ===================================================
	// 				LOAD AND PARSE THE FILE

	char* JSON_STRING = readFile(Input->breakPointFile);

	int i = 1;
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

	char* strValue = NULL; // adress where to fetch a value;
	//char* strToken = NULL; // adress where to fetch a value;

	while (i<r) {
		strValue = JSON_STRING+t[i+1].start;
		if 			(TOKEN("timeStep")) {
			Numerics->timeStep = atoi(strValue);
		} else if	(TOKEN("time")) {
			Physics->time  = atof(strValue);
		} else if	(TOKEN("dt")) {
			Physics->dt = atof(strValue);
		} else if	(TOKEN("dtAdv")) {
			Physics->dtAdv = atof(strValue);
		} else if	(TOKEN("dtT")) {
			Physics->dtT = atof(strValue);
		} else if	(TOKEN("dtPrevTimeStep")) {
			Numerics->dtPrevTimeStep = atof(strValue);
		} else if	(TOKEN("residual")) {
			Numerics->lsLastRes = atof(strValue);
		} else if	(TOKEN("Dresidual")) {
			// doesn't refer to a specific value
		} else if	(TOKEN("n_iterations")) {
			Numerics->itNonLin = atoi(strValue);
		} else if	(TOKEN("xmin")) {
			Grid->xmin = atof(strValue);
		} else if	(TOKEN("xmax")) {
			Grid->xmax = atof(strValue);
		} else if	(TOKEN("ymin")) {
			Grid->ymin = atof(strValue);
		} else if	(TOKEN("ymax")) {
			Grid->ymax = atof(strValue);
		} else if	(TOKEN("nxS")) {
			Grid->nxS = atoi(strValue);
		} else if	(TOKEN("nyS")) {
			Grid->nyS = atoi(strValue);
		} else if	(TOKEN("Char_length")) {
			Char->length = atof(strValue);
		} else if	(TOKEN("Char_time")) {
			Char->time = atof(strValue);
		} else if	(TOKEN("Char_mass")) {
			Char->mass = atof(strValue);
		} else if	(TOKEN("Char_temperature")) {
			Char->temperature = atof(strValue);
		} else if	(TOKEN("Output_counter")) {
			Output->counter = atoi(strValue);
		} else {
			printf("error in Output_readModelState: unknown token: %.*s\n", t[i].end-t[i].start, JSON_STRING + t[i].start);
			exit(0);
		}
		i+=2;
	}

	// Fill related values
	Grid->nxC = Grid->nxS-1;
	Grid->nyC = Grid->nyS-1;
}




void Breakpoint_createParticleSystem(Model* Model) {
	Particles* Particles = &(Model->Particles);
	Physics* Physics 	 = &(Model->Physics);
	Grid* Grid 			 = &(Model->Grid);
	Output* Output	 	 = &(Model->Output);


	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	// Open bin file

	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];
	sprintf(Folder_thistStep, "%sOut_%05i/", Output->breakpointFolder,Output->counter);

	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_x");
	FILE* file_x = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_x);
	fread(&Char_quantity, sizeof(double), 1, file_x);


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_y");
	FILE* file_y = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_y);
	fread(&Char_quantity, sizeof(double), 1, file_y);


#if (STORE_PARTICLE_POS_INI)
	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_xIni");
	FILE* file_xIni = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_xIni);
	fread(&Char_quantity, sizeof(double), 1, file_xIni);


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_yIni");
	FILE* file_yIni = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_yIni);
	fread(&Char_quantity, sizeof(double), 1, file_yIni);


#endif
	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_passive");
	FILE* file_passive = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_passive);
	fread(&Char_quantity, sizeof(double), 1, file_passive);


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_phase");
	FILE* file_phase = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_phase);
	fread(&Char_quantity, sizeof(double), 1, file_phase);


#if (STORE_PLASTIC_STRAIN)
	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_strain");
	FILE* file_strain = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_strain);
	fread(&Char_quantity, sizeof(double), 1, file_strain);
#endif


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_Sxx0");
	FILE* file_Sxx0 = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_Sxx0);
	fread(&Char_quantity, sizeof(double), 1, file_Sxx0);


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_Sxy0");
	FILE* file_Sxy0 = fopen(fname,"rb");
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_Sxy0);
	fread(&Char_quantity, sizeof(double), 1, file_Sxy0);


	sprintf(fname,"%s%s.bin",Folder_thistStep, "particles_timeLastPlastic");
	FILE* file_timeLP = fopen(fname,"rb"); // timeLastPlastic
	// Move the reading pointer over the header
	fread(&Particles->n , sizeof(int), 1, file_timeLP);
	fread(&Char_quantity, sizeof(double), 1, file_timeLP);
	
	// Read the header
	
	//fwrite(data, sizeof(float), iPart, fptr);
	SingleParticle *modelParticle = (SingleParticle *)malloc(sizeof(SingleParticle));
	Particles_initModelParticle(modelParticle);
	

	int iPart;
	for(iPart = 0;iPart < Particles->n;iPart++)
	{
		fread(&modelParticle->x, sizeof(double), 1, file_x);
		fread(&modelParticle->y, sizeof(double), 1, file_y);

#if (STORE_PARTICLE_POS_INI)
		fread(&modelParticle->xIni, sizeof(double), 1, file_xIni);
		fread(&modelParticle->yIni, sizeof(double), 1, file_yIni);
#endif

		fread(&modelParticle->sigma_xx_0, sizeof(double), 1, file_Sxx0);
		fread(&modelParticle->sigma_xy_0, sizeof(double), 1, file_Sxy0);
		
		fread(&modelParticle->phase, sizeof(double), 1, file_phase);
		fread(&modelParticle->passive, sizeof(double), 1, file_passive);


	#if (STORE_PLASTIC_STRAIN)
		fread(&modelParticle->strain, sizeof(double), 1, file_strain);
	#endif
	#if (STORE_TIME_LAST_PLASTIC)
		fread(&modelParticle->timeLastPlastic, sizeof(double), 1, file_timeLP);
	#endif

		Particles_findNodeForThisParticle(modelParticle, Grid);
		Particles_addSingleParticle(&Particles->linkHead[modelParticle->nodeId], modelParticle);
	}

	
	fclose(file_x);
	fclose(file_y);
#if (STORE_PARTICLE_POS_INI)
	fclose(file_xIni);
	fclose(file_yIni);
#endif
	fclose(file_passive);
	fclose(file_phase);
#if (STORE_PLASTIC_STRAIN)
	fclose(file_strain);
#endif
	fclose(file_Sxx0);
	fclose(file_Sxy0);
	fclose(file_timeLP);
	





	// Refill the sticky air
	//================================

	// Go through the grid to check for empty cells and fill them with air particles + do something about the topo cells
	int iy, ix, iNode;
	SingleParticle* thisParticle;
	Particles_initModelParticle(modelParticle);
	int nPCX = 2;
	int nPCY = 2;
	compute dxP = Grid->dx/nPCX;
	compute dyP = Grid->dy/nPCY;
	int iPx, iPy;
	srand(time(NULL));

	compute x, y;
	int* TopoNodes = (int*) malloc(Grid->nxS * sizeof(int));
	for (ix = 0; ix < Grid->nxS; ++ix) { 
		TopoNodes[ix] = 0;
		for (iy = 0; iy < Grid->nyS; ++iy) { 
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			x = Grid->X[ix];// - 0.5*Grid->dx;
			y = Grid->Y[iy];// - 0.5*Grid->dy;
			if (thisParticle==NULL) {
				for (iPy=0;iPy<nPCY;iPy++) {
					for (iPx=0;iPx<nPCX;iPx++) {
						modelParticle->x 	= x + 0.5*dxP + iPx*dxP + Particles->noiseFactor*dxP*(0.5 - (rand() % 1000)/1000.0);
						modelParticle->y 	= y + 0.5*dyP + iPy*dyP + Particles->noiseFactor*dyP*(0.5 - (rand() % 1000)/1000.0);

						if (modelParticle->x<Grid->xmin || modelParticle->x>Grid->xmax || modelParticle->y<Grid->ymin || modelParticle->y>Grid->ymax) {
							// do nothing
						} else {
#if (STORE_PARTICLE_POS_INI)
							modelParticle->xIni = modelParticle->x;
							modelParticle->yIni = modelParticle->y;
#endif
							modelParticle->phase = Physics->phaseAir;
							modelParticle->nodeId = ix + iy*Grid->nxS;
							Particles_addSingleParticle(&Particles->linkHead[modelParticle->nodeId], modelParticle);
							Particles->n++;
						}
					}
				}	
			} else {
				if (iy>TopoNodes[ix]) {
					TopoNodes[ix] = iy;
				}
			}
		}
	}

	// Loop through TopoNodes, check each quadrant and add a single air particle
	bool QuadrantIsEmpty[4];
	int iQuad;
	compute locX, locY;
	compute xMod[4] = { 1.0,-1.0,-1.0, 1.0};
	compute yMod[4] = { 1.0, 1.0,-1.0,-1.0};
	for (ix = 0; ix < Grid->nxS; ++ix) { 
		iy = TopoNodes[ix];
		iNode = ix + iy*Grid->nxS;
		thisParticle = Particles->linkHead[iNode];
		QuadrantIsEmpty[0] = true; // NE
		QuadrantIsEmpty[1] = true; // NW
		QuadrantIsEmpty[2] = true; // SW
		QuadrantIsEmpty[3] = true; // SE
		// Check the filling of quadrants
		while (thisParticle!=NULL) {
			locX = Particles_getLocX(ix, thisParticle->x,Grid);
			locY = Particles_getLocY(iy, thisParticle->y,Grid);
			if (locX>0) {
				if (locY>0) {
					QuadrantIsEmpty[0] = false;
				} else {
					QuadrantIsEmpty[2] = false;
				}
			} else {
				if (locY>0) {
					QuadrantIsEmpty[1] = false;
				} else {
					QuadrantIsEmpty[3] = false;
				}
			}
			thisParticle = thisParticle->next;
		}

		// Add air particles in the center of empty quadrants
		for(iQuad = 0;iQuad < 4;iQuad++) {
			if (QuadrantIsEmpty[iQuad]) {
				modelParticle->x = Grid->X[ix] + xMod[iQuad]*0.5*Grid->dx;
				modelParticle->xIni = modelParticle->x;
				modelParticle->y = Grid->Y[iy] + yMod[iQuad]*0.5*Grid->dy;
				modelParticle->yIni = modelParticle->y;
				modelParticle->phase = Physics->phaseAir;
				modelParticle->nodeId = ix + iy*Grid->nxS;
				Particles_addSingleParticle(&Particles->linkHead[modelParticle->nodeId], modelParticle);
				Particles->n++;
			}
		}
	}
	
	
	// Collect the data and fill the matrix
	//fwrite(MatrixToFill, sizeof(double), nxy[0]*nxy[1], fptr);

	// Close the file
	free(modelParticle);
}





