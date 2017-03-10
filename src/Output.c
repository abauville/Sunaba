/*
 * Output.c
 *
 *  Created on: 12 Sep 2016
 *      Author: abauville
 */



#include "stokes.h"

// Works only for UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


void Output_free(Output* Output) {
	//free(Output->ModelDescription); // assigned in Input_readVisu
}

void Output_writeInputCopyInOutput(Output* Output, Input* Input)
{
	// Writes a copy of the input file in the output folder
	char* InputFileString = readFile(Input->inputFile);


	// touch folder for this timestep
	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_Input[MAX_STRING_LENGTH];


	//sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_Input, "%sInput/", Output->outputFolder);


	//printf("filename: %smodelState.json\n",Folder_Input);



	struct stat st = {0};

	if (stat(Folder_Input, &st) == -1) {
		mkdir(Folder_Input, 0700);
	}


	sprintf(fname,"%sinput.json",Folder_Input);
	if ((fptr = fopen(fname,"w")) == NULL) {
		fprintf(stderr,"Failed the output file\n");
		exit(0);
	}


	fprintf(fptr,"%s", InputFileString);

	fclose(fptr);

	free(InputFileString);



}


void Output_modelState(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics)
{

	// touch folder for this timestep
	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_thistStep[MAX_STRING_LENGTH];


	//sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);


	printf("filename: %smodelState.json\n",Folder_thistStep);



	struct stat st = {0};

	if (stat(Folder_thistStep, &st) == -1) {
		mkdir(Folder_thistStep, 0700);
	}


	sprintf(fname,"%smodelState.json",Folder_thistStep);
	if ((fptr = fopen(fname,"w")) == NULL) {
		fprintf(stderr,"Failed the output file\n");
		exit(0);
	}

	fprintf(fptr,"{ \n");
	fprintf(fptr,"\t \"timeStep\"	: %i     			,\n", Numerics->timeStep);
	fprintf(fptr,"\t \"time\" 		: %f     			,\n", Physics->time);
	fprintf(fptr,"\t \"dt\" 		: %f     			,\n", Physics->dt);
	fprintf(fptr,"\t \"residual\"	: %f   				,\n", Numerics->lsLastRes);
	fprintf(fptr,"\t \"xmin\"		: %f   				,\n", Grid->xmin);
	fprintf(fptr,"\t \"xmax\"		: %f   				,\n", Grid->xmax);
	fprintf(fptr,"\t \"ymin\"		: %f   				,\n", Grid->ymin);
	fprintf(fptr,"\t \"ymax\"		: %f   				,\n", Grid->ymax);
	fprintf(fptr,"\t \"nxS\"		: %i  				,\n", Grid->nxS);
	fprintf(fptr,"\t \"nyS\"		: %i   				,\n", Grid->nyS);
	fprintf(fptr,"\t \"Char_length\"		: %f   			,\n", Char->length);
	fprintf(fptr,"\t \"Char_time\"  		: %f   			,\n", Char->mass);
	fprintf(fptr,"\t \"Char_mass\"  		: %f   			,\n", Char->time);
	fprintf(fptr,"\t \"Char_temperature\" 	: %f   			 \n", Char->temperature);
	//	fprintf(fptr,"\t \"Description\" 		: %s   			 \n", Output->ModelDescription);

	fprintf(fptr,"}");

	fclose(fptr);

	//printf("%s\n",Output->ModelDescription);






}


void Output_data(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics)
{

	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_thistStep[MAX_STRING_LENGTH];
	char Data_name[MAX_STRING_LENGTH];

	double* PointerToData;

	int iOut;

	//sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);


	//Output->nTypes = 1;
	//Output->type[0] = Out_Viscosity;
	//Output->type[1] = Out_Vy;

	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int iy, ix, iCell, iNode;
	for (iOut = 0; iOut < Output->nTypes; ++iOut) {
		compute* Data;
		printf("iOut = %i, Type = %d\n",iOut, Output->type[iOut]);
		switch (Output->type[iOut]) {
		case Out_Vx:
			nxy[0] = Grid->nxVx;
			nxy[1] = Grid->nyVx;
			xmin = Grid->xmin;
			xmax = Grid->xmax;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		case Out_Vy:
			nxy[0] = Grid->nxVy;
			nxy[1] = Grid->nyVy;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin;
			ymax = Grid->ymax;
			break;
		case Out_P:
		case Out_Pf:
		case Out_Pc:
		case Out_Viscosity:
		case Out_Porosity:
		case Out_Z:
		case Out_G:
		case Out_Khi:
		case Out_Sxx0:
		case Out_Sxx:
		case Out_StrainRate:
		case Out_SII:
		case Out_Temperature:
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		case Out_Sxy:
		case Out_Sxy0:
			nxy[0] = Grid->nxS;
			nxy[1] = Grid->nyS;
			xmin = Grid->xmin;
			xmax = Grid->xmax;
			ymin = Grid->ymin;
			ymax = Grid->ymax;
			break;

			//Out_Khi, Out_Sxx0, OutSxy0, Out_StrainRate
		default:
			printf("error: Unknown Output type");
			exit(0);
		}
		switch (Output->type[iOut]) {
		case Out_Vx:
			sprintf(Data_name,"Vx");
			PointerToData = Physics->Vx;
			Char_quantity = Char->length / Char->time;
			break;
		case Out_Vy:
			sprintf(Data_name,"Vy");
			PointerToData = Physics->Vy;
			Char_quantity = Char->length / Char->time;
			break;
		case Out_P:
			sprintf(Data_name,"P");
			PointerToData = Physics->P;
			Char_quantity = Char->stress;
			break;
		case Out_Pf:
#if (DARCY)
			sprintf(Data_name,"Pf");
			PointerToData = Physics->Pf;
			Char_quantity = Char->stress;
#endif
			break;
		case Out_Pc:
#if (DARCY)
			sprintf(Data_name,"Pc");
			PointerToData = Physics->Pc;
			Char_quantity = Char->stress;
#endif
			break;
		case Out_Viscosity:
			sprintf(Data_name,"Viscosity");
			PointerToData = Physics->eta;
			Char_quantity = Char->stress * Char->time;
			break;
		case Out_Porosity:
#if (DARCY)
			sprintf(Data_name,"Porosity");
			PointerToData = Physics->phi;
			Char_quantity = 1.0;
#endif
			break;
			//Out_Khi, Out_Sxx0, OutSxy0, Out_StrainRate
		case Out_Z:
			sprintf(Data_name,"Z");
			PointerToData = Physics->Z;
			Char_quantity = Char->stress * Char->time;
			break;
		case Out_G:
			sprintf(Data_name,"G");
			PointerToData = Physics->G;
			Char_quantity = Char->stress;
			break;
		case Out_Khi:
			sprintf(Data_name,"Khi");
			PointerToData = Physics->khi;
			Char_quantity = Char->stress * Char->time;
			break;
		case Out_Sxx0:
			sprintf(Data_name,"sigma_xx0");
			PointerToData = Physics->sigma_xx_0;
			Char_quantity = Char->stress;
			break;
		case Out_Sxy0:
			sprintf(Data_name,"sigma_xy0");
			PointerToData = Physics->sigma_xy_0;
			Char_quantity = Char->stress;
			break;
		case Out_Sxx:
			sprintf(Data_name,"sigma_xx");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					Data[iCell] = Physics->sigma_xx_0[iCell] + Physics->Dsigma_xx_0[iCell];
				}
			}
			Char_quantity = Char->stress;
			break;
		case Out_Sxy:
			sprintf(Data_name,"sigma_xy");
			Data = (compute*) malloc(Grid->nSTot * sizeof(compute));
			PointerToData = Data;
			for (iy = 0; iy < Grid->nyS; ++iy) {
				for (ix = 0; ix < Grid->nxS; ++ix) {
					iNode = ix + iy*Grid->nxS;
					Data[iNode] = Physics->sigma_xy_0[iNode] + Physics->Dsigma_xy_0[iNode];
				}
			}
			Char_quantity = Char->stress;
			break;
		case Out_SII:
			sprintf(Data_name,"sigma_II");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			compute SII;
			for (iy = 1; iy < Grid->nyEC-1; ++iy) {
				for (ix = 1; ix < Grid->nxEC-1; ++ix) {
					Physics_computeStressInvariantForOneCell(Physics, Grid, ix, iy, &SII);
					Data[ix + iy*Grid->nxEC] = SII;
				}
			}
			Physics_copyValuesToSides(Data, Grid);
			Char_quantity = Char->stress;
			break;
		case Out_StrainRate:
			sprintf(Data_name,"strainRate");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			Physics_computeStrainRateInvariant(Physics, Grid, Data);
			Physics_copyValuesToSides(Data, Grid);
			Char_quantity = 1.0 / Char->time;
			break;
		case Out_Temperature:
#if (HEAT)
			sprintf(Data_name,"Temperature");
			PointerToData = Physics->T;
			Char_quantity = Char->stress * Char->time;
			break;
#endif
		default:
			printf("error: Unknown Output type");
			exit(0);
		}


		printf("filename: %s%s.bin\n",Folder_thistStep, Data_name);



		struct stat st = {0};

		if (stat(Folder_thistStep, &st) == -1) {
			mkdir(Folder_thistStep, 0700);
		}


		sprintf(fname,"%s%s.bin",Folder_thistStep, Data_name);
		if ((fptr = fopen(fname,"w")) == NULL) {
			fprintf(stderr,"Failed the output file\n");
			exit(0);
		}


		fwrite(nxy , sizeof(int), 2, fptr);
		fwrite(&xmin, sizeof(double), 1, fptr);
		fwrite(&xmax, sizeof(double), 1, fptr);
		fwrite(&ymin, sizeof(double), 1, fptr);
		fwrite(&ymax, sizeof(double), 1, fptr);
		fwrite(&Char_quantity, sizeof(double), 1, fptr);
		fwrite(PointerToData, sizeof(double), nxy[0]*nxy[1], fptr);

		fclose(fptr);


		if (Output->type[iOut] == Out_Sxx || Output->type[iOut] == Out_Sxy || Output->type[iOut] == Out_SII || Output->type[iOut] == Out_StrainRate) {
			free(Data);
		}

	}





}







void Output_particles(Output* Output, Particles* Particles, Char* Char)
{


	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_thistStep[MAX_STRING_LENGTH];
	char Data_name[MAX_STRING_LENGTH];

	double* PointerToData;

	int iOut;

	//sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);


	//Output->nTypes = 1;
	//Output->type[0] = Out_Viscosity;
	//Output->type[1] = Out_Vy;

	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int iy, ix, iCell, iNode;
	for (iOut = 0; iOut < Output->nPartTypes; ++iOut) {
		compute* Data;
		printf("iOut = %i, Type = %d\n",iOut, Output->partType[iOut]);

	/*
// Particles
// =========================

// Single Particle storing coordinate, temp and info for a linked list

typedef struct SingleParticle SingleParticle;
struct SingleParticle {
	coord x, y;
	int phase;
	float passive; // some passive attribute used for visualization

#if (HEAT)
	compute T;
#endif

	// Old stresses
	compute sigma_xx_0;
	compute sigma_xy_0;

#if (DARCY)
	compute DeltaP0;
	compute phi;
#endif
	//bool faulted;

#if (STORE_PARTICLE_POS_INI)
	float xIni, yIni;
#endif


	// for the linked list
	int nodeId;
    SingleParticle* next;

};
*/


		int dataSize;
		float* 	dataFloat 	= (float*) 	malloc(Particles->n * sizeof(float));
		int* 	dataIntt 	= (int*) 	malloc(Particles->n * sizeof(int)  );


		switch (Output->type[iOut]) {
		case OutPart_x:
			sprintf(Data_name,"particles_x");
			Char_quantity = Char->length;
			break;
		case OutPart_y:
			sprintf(Data_name,"particles_y");
			Char_quantity = Char->length;
			break;
		case OutPart_xIni:
			sprintf(Data_name,"particles_xIni");
			Char_quantity = Char->length;
			break;
		case OutPart_yIni:
			sprintf(Data_name,"particles_yIni");
			Char_quantity = Char->length;
			break;
		case OutPart_Phase:
			sprintf(Data_name,"particles_phase");
			Char_quantity = 1;
			break;
		case OutPart_Passive:
			sprintf(Data_name,"particles_passive");
			Char_quantity = 1.0;
			break;
		case OutPart_T:
#if (HEAT)
			sprintf(Data_name,"particles_T");
			Char_quantity = Char->temperature;
#endif
			break;
		case OutPart_P0:
#if (DARCY)
			sprintf(Data_name,"particles_DeltaP0");
			Char_quantity = Char->stress;

#endif
			break;
		case OutPart_Sxx0:
			sprintf(Data_name,"particles_Sxx0");
			Char_quantity = Char->stress;
			break;
		case OutPart_Sxy0:
			sprintf(Data_name,"particles_Sxy0");
			Char_quantity = Char->stress;
			break;
		case OutPart_Phi:
			sprintf(Data_name,"particles_phi");
			Char_quantity = 1.0;
			break;

		default:
			printf("error: Unknown Particle Output type");
			exit(0);
		}



		printf("filename: %s%s.bin\n",Folder_thistStep, Data_name);



		struct stat st = {0};

		if (stat(Folder_thistStep, &st) == -1) {
			mkdir(Folder_thistStep, 0700);
		}


		sprintf(fname,"%s%s.bin",Folder_thistStep, Data_name);
		if ((fptr = fopen(fname,"w")) == NULL) {
			fprintf(stderr,"Failed the output file\n");
			exit(0);
		}


		fwrite(&Particles->n , sizeof(int), 2, fptr);
		fwrite(&Char_quantity, sizeof(double), 1, fptr);
		fwrite(PointerToData, sizeof(double), nxy[0]*nxy[1], fptr);

		fclose(fptr);


		if (Output->type[iOut] == Out_Sxx || Output->type[iOut] == Out_Sxy || Output->type[iOut] == Out_SII || Output->type[iOut] == Out_StrainRate) {
			free(Data);
		}

	}

}
































































