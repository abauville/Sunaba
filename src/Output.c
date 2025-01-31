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
#include <stddef.h>


void Output_init(Output* Output) {
	Output->counter = 0;
}

void Output_free(Output* Output) {
	//free(Output->ModelDescription); // assigned in Input_readVisu
}

void Output_call(Model* Model) {
	Output* Output 			= &(Model->Output);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);
	bool writeOutput = false;
	if (Output->nTypes>0 || Output->nPartTypes>0) {
		
		if (Output->useTimeFrequency) {
			//printf("Output->counter*Output->timeFrequency = %.2e, tim = %.2e\n", Output->counter*Output->timeFrequency, Physics->time);
			if (Physics->time>Output->counter*Output->timeFrequency) {
				writeOutput = true;
			} else if (Output->saveFirstStep && Numerics->timeStep == 0) {
				writeOutput = true;
			}
		} else {
			if (((Numerics->timeStep) % Output->frequency)==0) {
				if (Numerics->timeStep>0 || Output->saveFirstStep) {
					writeOutput = true;
				}
			}
		}
		if (writeOutput) {
			printf("Write output ...\n");
			Output_modelState(Model,false);
			Output_data(Model);
			Output_particles(Model,false);
			Output->counter++;
		}
	}

	
}


void Output_writeInputCopyInOutput(Output* Output, Input* Input)
{
	// Writes a copy of the input file in the output folder
	char* InputFileString = readFile(Input->inputFile);


	// touch folder for this timestep
	FILE *fptr;
	char fname[MAX_STRING_LENGTH];
	char Folder_Input[BUFFER_STRING_LENGTH];

	sprintf(Folder_Input, "%sInput/", Output->outputFolder);

	struct stat st = {0};
	if (stat(Folder_Input, &st) == -1) {
		mkdir(Folder_Input, 0700);
	}

	sprintf(fname,"%sinput.json",Folder_Input);
	if ((fptr = fopen(fname,"w")) == NULL) {
		printf("error: Failed to copy the input file: %s\n", fname);
		exit(0);
	}

	fprintf(fptr,"%s", InputFileString);
	fclose(fptr);
	free(InputFileString);
}


void Output_modelState(Model* Model, bool restartMode)
{
	Output* Output 			= &(Model->Output);
	Breakpoint* Breakpoint 	= &(Model->Breakpoint);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);
	Numerics* Numerics 		= &(Model->Numerics);
	EqSystem* EqStokes 		= &(Model->EqStokes);
	
	// touch folder for this timestep
	FILE *fptr;
	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];


	if (restartMode) {
		sprintf(Folder_thistStep, "%sOut_%05i/", Breakpoint->breakpointFolder, Breakpoint->counter);
	} else {
		sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);
	}

	//printf("filename: %smodelState.json\n",Folder_thistStep);

	struct stat st = {0};

	if (stat(Folder_thistStep, &st) == -1) {
		mkdir(Folder_thistStep, 0700);
	}


	sprintf(fname,"%smodelState.json",Folder_thistStep);
	if ((fptr = fopen(fname,"w")) == NULL) {
		printf("error: Failed to output the modelState file: %s\n", fname);
		exit(0);
	}

	fprintf(fptr,"{ \n");
	fprintf(fptr,"\t \"timeStep\"	: %i     			,\n", Numerics->timeStep);
	fprintf(fptr,"\t \"time\" 		: %.14f     			,\n", Physics->time);
	fprintf(fptr,"\t \"dt\" 		: %.14f     			,\n", Physics->dt);
	fprintf(fptr,"\t \"dtAdv\" 		: %.14f     			,\n", Physics->dtAdv);
	fprintf(fptr,"\t \"dtT\" 		: %.14f     			,\n", Physics->dtT);
	fprintf(fptr,"\t \"dtPrevTimeStep\" 	: %.14f 		,\n", Numerics->dtPrevTimeStep);
	fprintf(fptr,"\t \"residual\"	: %.14f   				,\n", Numerics->lsLastRes);
	fprintf(fptr,"\t \"Dresidual\"	: %.14f   				,\n", fabs(EqStokes->normResidual-Numerics->oldRes));
	fprintf(fptr,"\t \"n_iterations\"		: %i   			,\n", Numerics->itNonLin);
	fprintf(fptr,"\t \"xmin\"		: %.14f   				,\n", Grid->xmin);
	fprintf(fptr,"\t \"xmax\"		: %.14f   				,\n", Grid->xmax);
	fprintf(fptr,"\t \"ymin\"		: %.14f   				,\n", Grid->ymin);
	fprintf(fptr,"\t \"ymax\"		: %.14f   				,\n", Grid->ymax);
	fprintf(fptr,"\t \"nxS\"		: %i  					,\n", Grid->nxS);
	fprintf(fptr,"\t \"nyS\"		: %i   					,\n", Grid->nyS);
	fprintf(fptr,"\t \"Char_length\"		: %.14f   			,\n", Char->length);
	fprintf(fptr,"\t \"Char_time\"  		: %.14f   			,\n", Char->time);
	fprintf(fptr,"\t \"Char_mass\"  		: %.14f   			,\n", Char->mass);
	fprintf(fptr,"\t \"Char_temperature\" 	: %.14f   			,\n", Char->temperature);
	fprintf(fptr,"\t \"Output_counter\" 	: %i     			,\n", Output->counter);
	fprintf(fptr,"\t \"Breakpoint_counter\" 	: %i     		,\n", Breakpoint->counter);
	fprintf(fptr,"\t \"realTimeSinceStart\" : %.14f    		 \n", Numerics->realTimeSinceStart);
	

	//	fprintf(fptr,"\t \"Description\" 		: %s   			 \n", Output->ModelDescription);

	fprintf(fptr,"}");

	fclose(fptr);


}


void Output_data(Model* Model)
{

	Output* Output 			= &(Model->Output);
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);
	

	FILE *fptr;
	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];
	char Data_name[BUFFER_STRING_LENGTH];

	double* PointerToData = NULL;

	int iOut;

	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);


	printf("Folder: %s\n",Folder_thistStep);

	int nxy[2];
	double Char_quantity;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int iy, ix, iCell;
	for (iOut = 0; iOut < Output->nTypes; ++iOut) {
		compute* Data = NULL;
		//printf("iOut = %i, Type = %d\n",iOut, Output->type[iOut]);
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
		case Out_Sxy0:
		case Out_Sxx:
		case Out_Sxy:
		case Out_StrainRate:
		case Out_Strain:
		case Out_SII:
		case Out_Temperature:
		case Out_Phase:
			nxy[0] = Grid->nxEC;
			nxy[1] = Grid->nyEC;
			xmin = Grid->xmin - Grid->dx/2.0;
			xmax = Grid->xmax + Grid->dx/2.0;
			ymin = Grid->ymin - Grid->dy/2.0;
			ymax = Grid->ymax + Grid->dy/2.0;
			break;
		
		case Out_Sxy_Node: // i.e. sxy or sxy0
			nxy[0] = Grid->nxS;
			nxy[1] = Grid->nyS;
			xmin = Grid->xmin;
			xmax = Grid->xmax;
			ymin = Grid->ymin;
			ymax = Grid->ymax;
			break;
		
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
			sprintf(Data_name,"viscosity");
			PointerToData = Physics->eta;
			Char_quantity = Char->stress * Char->time;
			break;
		case Out_Porosity:
#if (DARCY)
			sprintf(Data_name,"porosity");
			PointerToData = Physics->phi;
			Char_quantity = 1.0;
#endif
			break;
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
			sprintf(Data_name,"khi");
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

			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			for (iy = 1; iy < Grid->nyEC-1; ++iy) {
				for (ix = 1; ix < Grid->nxEC-1; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					Data[iCell] = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0 ,ix,iy,Grid->nxS);
				}
			}
			Physics_CellVal_SideValues_copyNeighbours_Global(Data, Grid);
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
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			
			compute sxy0, Dsxy0;
			for (iy = 1; iy < Grid->nyEC-1; ++iy) {
				for (ix = 1; ix < Grid->nxEC-1; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					sxy0 	= Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0 ,ix,iy,Grid->nxS);
					Dsxy0 	= Interp_NodeVal_Node2Cell_Local(Physics->Dsigma_xy_0,ix,iy,Grid->nxS);
					Data[iCell] = sxy0 + Dsxy0;
				}
			}
			Physics_CellVal_SideValues_copyNeighbours_Global(Data, Grid);
			Char_quantity = Char->stress;
			break;
			case Out_Sxy_Node:
				sprintf(Data_name,"sigma_xy_node");
				Data = (compute*) malloc(Grid->nSTot * sizeof(compute));
				PointerToData = Data;
				int iNode;
				for (iy = 0; iy < Grid->nyS; ++iy) {
					for (ix = 0; ix < Grid->nxS; ++ix) {
						iNode = ix + iy*Grid->nxS;
						Data[iNode] = Physics->sigma_xy_0[iNode] + Physics->Dsigma_xy_0[iNode];
					}
				}
				break;
		case Out_SII:
			sprintf(Data_name,"sigma_II");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			compute SII;
			for (iy = 1; iy < Grid->nyEC-1; ++iy) {
				for (ix = 1; ix < Grid->nxEC-1; ++ix) {
					SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);
					Data[ix + iy*Grid->nxEC] = SII;
				}
			}
			Physics_CellVal_SideValues_copyNeighbours_Global(Data, Grid);
			Char_quantity = Char->stress;
			break;
		case Out_StrainRate:
			sprintf(Data_name,"strainRate");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			for (iy = 1; iy < Grid->nyEC-1; ++iy) {
				for (ix = 1; ix < Grid->nxEC-1; ++ix) {
					Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &SII);
					Data[ix + iy*Grid->nxEC] = SII;
				}
			}
			Physics_CellVal_SideValues_copyNeighbours_Global(Data, Grid);
			Char_quantity = 1.0 / Char->time;
			break;
		case Out_Strain:
#if (STORE_PLASTIC_STRAIN)
			sprintf(Data_name,"strain");
			PointerToData = Physics->strain;
			Char_quantity = 1.0;
#endif
			break;
		case Out_Temperature:
#if (HEAT)
			sprintf(Data_name,"temperature");
			PointerToData = Physics->T;
			Char_quantity = Char->stress * Char->time;
			break;
#endif
		case Out_Phase:
			sprintf(Data_name,"phase");
			Data = (compute*) malloc(Grid->nECTot * sizeof(compute));
			PointerToData = Data;
			for (iy = 0; iy < Grid->nyEC; ++iy) {
				for (ix = 0; ix < Grid->nxEC; ++ix) {
					iCell = ix + iy*Grid->nxEC;
					Data[iCell] = (compute) Physics->phase[iCell];
				}
			}
			Char_quantity = 1.0;
			break;
		default:
			printf("error: Unknown Output type");
			exit(0);
		}


		//printf("filename: %s%s.bin\n",Folder_thistStep, Data_name);

		if (iOut==0) {
			printf("Output data: \n");
		}
		
		printf("Data: %s\n",Data_name);


		struct stat st = {0};

		if (stat(Folder_thistStep, &st) == -1) {
			mkdir(Folder_thistStep, 0700);
		}


		sprintf(fname,"%s%s.bin",Folder_thistStep, Data_name);
		if ((fptr = fopen(fname,"w")) == NULL) {
			printf("error: Failed to output the data file: %s\n", fname);
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


		if (Output->type[iOut] == Out_Sxx || Output->type[iOut] == Out_Sxy0 || Output->type[iOut] == Out_Sxy || Output->type[iOut] == Out_SII || Output->type[iOut] == Out_StrainRate || Output->type[iOut] == Out_Phase || Output->type[iOut] == Out_Sxy_Node) {
			free(Data);
		}

	}

}


void Output_particles(Model* Model, bool breakpointMode)
{
	Output* Output 			= &(Model->Output);
	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	Char* Char 				= &(Model->Char);
	Physics* Physics 		= &(Model->Physics);


	int phaseAir = Physics->phaseAir;
	

	FILE *fptr;
	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];
	char Data_name[BUFFER_STRING_LENGTH];


	int iOut;

	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);

	

	double Char_quantity;
	//int iCell;
	for (iOut = 0; iOut < Output->nPartTypes; ++iOut) {
		//printf("iOut = %i, Type = %d\n",iOut, Output->partType[iOut]);

		INIT_PARTICLE;
		int dataOffset = 0;
		float* 	data = (float*) 	malloc(Particles->n * sizeof(float));
		int thisType = -1; // 0 = double, 1 = float, 2 = int

		switch (Output->partType[iOut]) {
		case OutPart_x:
			sprintf(Data_name,"particles_x");
			Char_quantity = Char->length;
			dataOffset = offsetof(SingleParticle, x);
			thisType = 0;
			//printf("offset = %i\n",dataOffset);
			break;
		case OutPart_y:
			sprintf(Data_name,"particles_y");
			Char_quantity = Char->length;
			dataOffset = offsetof(SingleParticle, y);
			thisType = 0;
			break;
		case OutPart_xIni:
#if (STORE_PARTICLE_POS_INI)
			sprintf(Data_name,"particles_xIni");
			Char_quantity = Char->length;
			dataOffset = offsetof(SingleParticle, xIni);
			thisType = 1;
#endif
			break;
		case OutPart_yIni:
#if (STORE_PARTICLE_POS_INI)
			sprintf(Data_name,"particles_yIni");
			Char_quantity = Char->length;
			dataOffset = offsetof(SingleParticle, yIni);
			thisType = 1;
#endif
			break;
		case OutPart_Phase:
			sprintf(Data_name,"particles_phase");
			Char_quantity = 1;
			dataOffset = offsetof(SingleParticle, phase);
			thisType = 2;
			break;
		case OutPart_Passive:
			sprintf(Data_name,"particles_passive");
			Char_quantity = 1.0;
			dataOffset = offsetof(SingleParticle, passive);
			thisType = 1;
			break;
		case OutPart_T:
#if (HEAT)
			sprintf(Data_name,"particles_T");
			Char_quantity = Char->temperature;
			dataOffset = offsetof(SingleParticle, T);
			thisType = 0;
#endif
			break;
		case OutPart_DeltaP0:
#if (DARCY)
			sprintf(Data_name,"particles_DeltaP0");
			Char_quantity = Char->stress;
			dataOffset = offsetof(SingleParticle, DeltaP0);
			thisType = 0;
#endif
			break;
		case OutPart_Sxx0:
			sprintf(Data_name,"particles_Sxx0");
			Char_quantity = Char->stress;
			dataOffset = offsetof(SingleParticle, sigma_xx_0);
			thisType = 0;
			break;
		case OutPart_Sxy0:
			sprintf(Data_name,"particles_Sxy0");
			Char_quantity = Char->stress;
			dataOffset = offsetof(SingleParticle, sigma_xy_0);
			thisType = 0;
			break;
		case OutPart_Phi:
#if (DARCY)
			sprintf(Data_name,"particles_phi");
			Char_quantity = 1.0;
			dataOffset = offsetof(SingleParticle, phi);
			thisType = 0;
#endif
			break;
		case OutPart_Strain:
#if (STORE_PLASTIC_STRAIN)
			sprintf(Data_name,"particles_strain");
			Char_quantity = 1.0;
			dataOffset = offsetof(SingleParticle, strain);
			thisType = 0;
#endif
			break;
		case OutPart_Vorticity_cum:
#if (STORE_PLASTIC_STRAIN)
			sprintf(Data_name,"particles_vorticity_cumulated");
			Char_quantity = 1.0/Char->time;
			dataOffset = offsetof(SingleParticle, vorticity_cum);
			thisType = 0;
#endif
			break;
		case OutPart_TimeLastPlastic:
#if (STORE_TIME_LAST_PLASTIC)
			sprintf(Data_name,"particles_timeLastPlastic");
			Char_quantity = Char->time;
			dataOffset = offsetof(SingleParticle, timeLastPlastic);
			thisType = 0;
#endif
			break;
		default:
			thisType = -1;
			printf("error: Unknown Particle Output type");
			printf("iOut = %i, PartType = %d", iOut, Output->partType[iOut]);
			exit(0);
		}



		int iPart = 0;

		FOR_PARTICLES
			char* base = (char*) thisParticle;
			if (breakpointMode || thisParticle->phase!=phaseAir) {
				if (thisType==0) {
					compute* ptr2value = (compute*)(base+dataOffset);
					data[iPart] = (float) *ptr2value;
				} else if (thisType == 1) {
					float* ptr2value = (float*)(base+dataOffset);
					data[iPart] = (float) *ptr2value;
				} else if (thisType == 2) {
					int* ptr2value = (int*)(base+dataOffset);
					data[iPart] = (float) *ptr2value;
				} else {
					printf("error in OutputPart: unknwon thisType = %i",thisType);
					exit(0);
				}
				iPart++;
			}
		END_PARTICLES

		if (iOut==0) {
			printf("Output particles: nPart = %i, Particles->n=%i\n",iPart, Particles->n);
		}
		
		printf("Data: %s\n",Data_name);
	
		struct stat st = {0};

		if (stat(Folder_thistStep, &st) == -1) {
			mkdir(Folder_thistStep, 0700);
		}


		sprintf(fname,"%s%s.bin",Folder_thistStep, Data_name);
		if ((fptr = fopen(fname,"w")) == NULL) {
			printf("error: Failed to output the particle file: %s\n", fname);
			exit(0);
		}


		//fwrite(&Particles->n , sizeof(int), 1, fptr);
		fwrite(&iPart , sizeof(int), 1, fptr);
		fwrite(&Char_quantity, sizeof(double), 1, fptr);
		fwrite(data, sizeof(float), iPart, fptr);

		fclose(fptr);

		free(data);
		
	}
}






void Output_particleBoundaryData(Model* Model)
{

	Output* Output 			= &(Model->Output);
	Grid* Grid 				= &(Model->Grid);
	Particles* Particles 	= &(Model->Particles);
	

	FILE *fptr;
	char fname[BUFFER_STRING_LENGTH];
	char Folder_thistStep[BUFFER_STRING_LENGTH];
	


	int iOut;

	sprintf(Folder_thistStep, "%sOut_%05i/", Output->outputFolder,Output->counter);

	int nBoundPassive = Particles->boundPassiveGridRefinement  * (Grid->nyS-1) + 1;
	
	int iy, ix, iCell;

	struct stat st = {0};

	if (stat(Folder_thistStep, &st) == -1) {
		mkdir(Folder_thistStep, 0700);
	}


	sprintf(fname,"%sparticle_boundInfo.bin",Folder_thistStep);
	if ((fptr = fopen(fname,"w")) == NULL) {
		printf("error: Failed to output the particle boundary file: %s\n", fname);
		exit(0);
	}

	fwrite(&nBoundPassive, sizeof(int), 1, fptr);
	fwrite(Particles->dispAtBoundL, sizeof(compute), nBoundPassive, fptr);
	fwrite(Particles->dispAtBoundR, sizeof(compute), nBoundPassive, fptr);
	fwrite(Particles->currentPassiveAtBoundL, sizeof(int), nBoundPassive, fptr);
	fwrite(Particles->currentPassiveAtBoundR, sizeof(int), nBoundPassive, fptr);
	fclose(fptr);

}
