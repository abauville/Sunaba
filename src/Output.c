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
		case Out_SII:
			sprintf(Data_name,"sigma_II");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			int iy, ix;
			compute SII;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nyEC; ++ix) {
					Physics_computeStressInvariantForOneCell(Physics, Grid, ix, iy, &SII);
					Data[ix + iy*Grid->nxEC] = SII;
				}
			}
			Char_quantity = Char->stress;
			break;
		case Out_StrainRate:
			sprintf(Data_name,"strainRate");
			PointerToData = Physics->sigma_xy_0;
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			Physics_computeStrainRateInvariant(Physics, Grid, Data);
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


		if (Output->type[iOut] == Out_SII || Output->type[iOut] == Out_StrainRate) {
			free(Data);
		}

	}





}

