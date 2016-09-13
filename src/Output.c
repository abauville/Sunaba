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

void Output_modelState(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics)
{

	// touch folder for this timestep
	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_thistStep[MAX_STRING_LENGTH];


	sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_thistStep, "%stimeStep_%05i/", Output->outputFolder,Numerics->timeStep);


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


	fprintf(fptr,"}");

	fclose(fptr);


}


void Output_data(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics)
{
	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_thistStep[MAX_STRING_LENGTH];
	char Data_name[MAX_STRING_LENGTH];

	double* PointerToData;

	int iOut;



	sprintf(Output->outputFolder,"/Users/abauville/Work/Output_StokesFD/Test00/");
	sprintf(Folder_thistStep, "%stimeStep_%05i/", Output->outputFolder,Numerics->timeStep);


	Output->nOutputs = 1;
	Output->OutputType[0] = Out_Viscosity;
	//Output->OutputType[1] = Out_Vy;

	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	for (iOut = 0; iOut < Output->nOutputs; ++iOut) {
		switch (Output->OutputType[iOut]) {
		case Out_Vx:
			sprintf(Data_name,"Vx");
			PointerToData = Physics->Vx;
			nxy[0] = Grid->nxVx;
			nxy[1] = Grid->nyVx;
			Char_quantity = Char->length / Char->time;
			xmin = Grid->xmin;
			xmax = Grid->xmax;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		case Out_Vy:
			sprintf(Data_name,"Vy");
			PointerToData = Physics->Vy;
			nxy[0] = Grid->nxVy;
			nxy[1] = Grid->nyVy;
			Char_quantity = Char->length / Char->time;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin;
			ymax = Grid->ymax;
			break;
		case Out_P:
			sprintf(Data_name,"P");
			PointerToData = Physics->P;
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			Char_quantity = Char->stress;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		case Out_Pf:
#if (DARCY)
			sprintf(Data_name,"Pf");
			PointerToData = Physics->Pf;
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			Char_quantity = Char->stress;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
#endif
			break;
		case Out_Pc:
#if (DARCY)
			sprintf(Data_name,"Pc");
			PointerToData = Physics->Pc;
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			Char_quantity = Char->stress;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
#endif
			break;
		case Out_Viscosity:
			sprintf(Data_name,"Viscosity");
			PointerToData = Physics->eta;
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			Char_quantity = Char->stress / Char->time;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		case Out_Porosity:
#if (DARCY)
			sprintf(Data_name,"Porosity");
			PointerToData = Physics->phi;
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			Char_quantity = 1.0;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
#endif
			break;
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

	}





}

