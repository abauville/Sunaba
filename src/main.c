/*
 ============================================================================
 Name        : main.c
 Author      : Arthur Bauville
 Version     :
 Copyright   : 
 Description :
 ============================================================================
 */

#include "stokes.h"



int main(int argc, char *argv[]) {

	printf("\n\n\n\n\n\n");
	printf("               ============================\n"
			"               ============================\n");
	printf("\n\n\n\n\n\nBeginning of the program\n");
	printf("Num procs = %i\n",omp_get_num_procs());

#if (DARCY)
	int iy, ix, iCell;
#endif

	
	// Declare structures
	// =================================
	
	Model Model;
	// Declare structures
	// =================================
	// General
	Grid* Grid 				= &(Model.Grid);
	MatProps* MatProps 		= &(Model.MatProps);
	Particles* Particles 	= &(Model.Particles);
	Physics* Physics 		= &(Model.Physics);
	Char* Char 				= &(Model.Char);
	// Stokes: Conservation of momentum + continuity (+- Darcy)
	Numbering* NumStokes 	= &(Model.NumStokes);
	BC* BCStokes 			= &(Model.BCStokes);
	EqSystem* EqStokes		= &(Model.EqStokes);
	Solver* SolverStokes 	= &(Model.SolverStokes);
#if (DARCY)
	IC* ICDarcy 			= &(Model.ICDarcy);
#endif
	// Heat conservation
	
#if (HEAT)
	IC* ICThermal 			= &(Model.ICThermal);
	BC* BCThermal 			= &(Model.BCThermal);
#endif
	
#if (HEAT)
	EqSystem* EqThermal  	= &(Model.EqThermal);
	Numbering* NumThermal 	= &(Model.NumThermal);
	Solver* SolverThermal 	= &(Model.SolverThermal);
#endif
	// Numerics
	Numerics* Numerics 		= &(Model.Numerics);
	// Visu
#if (VISU)
	Visu* Visu 				= &(Model.Visu);
#endif
	// Input/Output
	Input* Input 			= &(Model.Input);
	Output* Output 			= &(Model.Output);
	Breakpoint* Breakpoint 	= &(Model.Breakpoint);





	//INIT_TIMER
	double toc;
	double globalTic = omp_get_wtime();
	double globalToc;
	//INIT_GLOBAL_TIMER
	//GLOBAL_TIC

	strncpy(Input->currentFolder, argv[0],strlen(argv[0])-strlen("StokesFD") );
	Input->currentFolder[strlen(argv[0])-strlen("StokesFD")] = '\0';
	printf("Input->currentFolder = %s\n",Input->currentFolder);

	if (argc < 2) {
		strcpy(Input->inputFile,INPUT_FILE);
	} else {
		strcpy(Input->inputFile,argv[1]);
	}

	if (argc >= 3) {
		Breakpoint->startingNumber = atoi(argv[2]);
		Breakpoint->use = true;
	} else {
		Breakpoint->startingNumber = -1;
		Breakpoint->use = false;
	}

	// User Input, Init etc...
	if (Breakpoint->use) {
		Main_Init_restart(&Model);
	} else {
		Main_Init_start(&Model);
	}
	// Set fields from file test
	//Visu->paused = true;
	//Visu_main(&Model);
	// Other variables
	// =================================
	int iEq;


	compute stressFacIni = Numerics->dt_stressFac;

//======================================================================================================
//======================================================================================================
//
//                          				TIME LOOP
//

	

	double timeStepTic;
	compute* NonLin_x0 = (compute*) malloc(EqStokes->nEq * sizeof(compute));
	compute* NonLin_dx = (compute*) malloc(EqStokes->nEq * sizeof(compute));
	
	//printf("Numerics->maxTime = %.2e, Physics->time = %.2e\n",Numerics->maxTime,Physics->time);
	while(Numerics->timeStep!=Numerics->nTimeSteps && Physics->time <= Numerics->maxTime) {
		printf("\n\n\n          ========  Time step %i, t= %3.2e yrs  ========   \n"
					 "              ===================================== \n\n",Numerics->timeStep, Physics->time*Char->time/(3600*24*365));
		printf("time0 = %.2e\n", Physics->time*Model.Char.time);
		Numerics->dtPrevTimeStep = Physics->dt;
		Numerics->dtAlphaCorr = Numerics->dtAlphaCorrIni;

		printf("dt = %.2e yrs, maxVx = %.2e cm/yr, maxVy = %.2e cm/yr, dx/maxVx = %.2e yrs, dy/maxVy = %.2e yrs", Physics->dt*Char->time/(3600*24*365), (Physics->maxVx*Char->length/Char->time) / (0.01/(3600*24*365)), (Physics->maxVy*Char->length/Char->time) / (0.01/(3600*24*365)), (Grid->dx/Physics->maxVx) * Char->time/(3600*24*365), (Grid->dy/Physics->maxVy) * Char->time/(3600*24*365));
#if (VISU)
		timeStepTic = omp_get_wtime();;
#endif
		Numerics->itNonLin = -1;


#if (HEAT)
		// save the value from the previous time step
		if (Numerics->itNonLin == -1) {
			for (i = 0; i < Grid->nECTot; ++i) {
				Physics->T0[i] = Physics->T[i];
			}
		}
#endif

		/*
		if(Physics->time*Char->time > 5e6 * (3600*24*365.25)) {
			//Physics->g[1] = 0.0;
			BCStokes->backStrainRate = 0.0;//
		}
		*/
		
		

		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 										NON-LINEAR ITERATION											//
		Numerics->stalling = false;
		Numerics->stallingCounter = 0;
		EqStokes->normResidual = 1.0;
		Numerics->normRes0 = 1.0;
		Numerics->normResRef = 1.0;
		Numerics->cumCorrection_fac = 0.0;
		Numerics->lsLastRes = 1E15;
		Numerics->lsGlob = 1.00;
		Numerics->lsBestRes = 1e15;
		Numerics->lsBestGlob = 1.0;
#if (VISU)
		Visu->nonLinItisOver = false;
#endif

		Numerics->itNonLin = 0;

		Numerics->lsLastRes = 1E100;

		Numerics->dt_stressFac = stressFacIni;
		Numerics->oneMoreIt = true;
		//Physics_dt_update(&Model);

#if (NON_LINEAR_VISC)
		while( ( (( (EqStokes->normResidual > Numerics->absoluteTolerance ) && Numerics->itNonLin<Numerics->maxNonLinearIter ) || Numerics->itNonLin<Numerics->minNonLinearIter)  || Numerics->cumCorrection_fac<=0.999   ) || Numerics->oneMoreIt) {
#else
	while(Numerics->oneMoreIt) {
#endif
			printf("\n\n  ==== Non linear iteration %i ==== \n",Numerics->itNonLin);
			Numerics->oneMoreIt = false;
			printf("dt = %.2e yrs\n",Physics->dt * (Char->time/(3600.0*24*365)));
			// =====================================================================================//
			//																						//
			// 										COMPUTE STOKES									//
			//if (Numerics->itNonLin<=0 || Physics->dt<1e-2 || Physics->dt>1e2) {
			if (Physics->dt<0.5 || Physics->dt>2.0) {
				//printf("Rescale\n");
				//Char_rescale(&Model, NonLin_x0);
			}
			memcpy(NonLin_x0, EqStokes->x, EqStokes->nEq * sizeof(compute));
			/*
			if (Numerics->timeStep<3) {
				EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, true, Numerics);
				EqSystem_solve(EqStokes, SolverStokes, BCStokes, NumStokes, &Model);
			} else {
			*/
				pardisoSolveStokesAndUpdatePlasticity(EqStokes, SolverStokes, BCStokes, NumStokes, &Model);
			//}
			Physics_Velocity_retrieveFromSolution(&Model);
			Physics_P_retrieveFromSolution(&Model);

			// 										COMPUTE STOKES									//
			//																						//
			// =====================================================================================//


#if (VISCOSITY_TYPE==1)
			printf("/!\\ /!\\ LINEAR_VISCOUS==true, Non-linear iterations are ineffective/!\\ \n");
			Physics_Velocity_retrieveFromSolution(&Model);
			Physics_P_retrieveFromSolution(&Model);
			Physics_Rho_updateGlobal(&Model);
			Physics_Eta_updateGlobal(&Model);

			break;
#elif (VISCOSITY_TYPE==2)
			Physics_Velocity_retrieveFromSolution(&Model);
			break;
#else // VISCOSITY_TYPE==0



			// =====================================================================================//
			//																						//
			// 										LINE SEARCH										//
			// Compute dx
			for (iEq = 0; iEq < EqStokes->nEq; ++iEq) {
				NonLin_dx[iEq] = EqStokes->x[iEq] - NonLin_x0[iEq];
			}

			Numerics->minRes = 1E100;
			Numerics->lsGlob = 1.0;
			Numerics->lsState = -1;
			Numerics->oldRes = EqStokes->normResidual;

			printf("dt = %.2e yrs\n",Physics->dt * (Char->time/(3600.0*24*365)));

#if (HEAT)
			// =====================================================================================//
			//																						//
			// 										COMPUTE HEAT									//
			//TIC
			Physics_Velocity_retrieveFromSolution(&Model);
			Physics_P_retrieveFromSolution(&Model);


#if (DARCY)
			Physics_Phi_updateGlobal(&Model);
			Physics_Perm_updateGlobal(&Model);
#endif

			Physics_Rho_updateGlobal(&Model);
			Physics_Eta_updateGlobal(&Model);
			printf("Heat assembly and solve\n");
			EqSystem_assemble(EqThermal, Grid, BCThermal, Physics, NumThermal, true, Numerics);

			EqSystem_scale(EqThermal);
			EqSystem_solve(EqThermal, SolverThermal, BCThermal, NumThermal, &Model);
			EqSystem_unscale(EqThermal);
			Physics_T_retrieveFromSolution(&Model);

			//TOC
			printf("Temp Assembly+Solve+Interp: %.3f s\n", toc);


			// 										COMPUTE HEAT									//
			//																						//
			// =====================================================================================//
#endif
#if (NON_LINEAR_VISC)
			int iLS;
			while (iLS < Numerics->nLineSearch+1) {
#pragma omp parallel for private(iEq) OMP_SCHEDULE
				for (iEq = 0; iEq < EqStokes->nEq; ++iEq) {
					EqStokes->x[iEq] = NonLin_x0[iEq] + Numerics->lsGlob*(NonLin_dx[iEq]);
				}

				Physics_Velocity_retrieveFromSolution(&Model);
				Physics_P_retrieveFromSolution(&Model);


#if (DARCY)
				Physics_Phi_updateGlobal(&Model);
				Physics_Perm_updateGlobal(&Model);
#endif
				Physics_Rho_updateGlobal(&Model);
				Physics_Eta_Simple_updateGlobal(&Model);


#if (DEBUG)
				Physics_check(&Model);
#endif
				EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, false, Numerics);
				EqSystem_computeNormResidual(EqStokes);

				printf("a = %.3f,  |Delta_Res| = %.2e, |F|/|b|: %.2e\n", Numerics->lsGlob, fabs(EqStokes->normResidual-Numerics->oldRes), EqStokes->normResidual);

				if (EqStokes->normResidual<Numerics->minRes) {
					Numerics->minRes = EqStokes->normResidual;
					Numerics->lsBestGlob = Numerics->lsGlob;
				}
				iLS++;
				if (iLS<Numerics->nLineSearch) {
					Numerics->lsGlob = Numerics->lsGlob/2.0;
				} else {
					if (Numerics->lsGlob == Numerics->lsBestGlob) {
						break;
					} else {
						Numerics->lsGlob = Numerics->lsBestGlob;
					}
				}

				if (Numerics->minRes<Numerics->lsLastRes) {
					break;
				}

				if (EqStokes->normResidual>1e10) {
					break;
				}

				if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
					printf("\n\n\n\n error: Something went wrong. The norm of the residual is NaN\n");
					break;
				}



				if (Numerics->lsGoingDown || Numerics->lsGoingUp) { // if the time step is changing there is no need to check other values of lsGlob (in most cases 1.0 is the best),
					// the residual goes up only because the time step is changeing. Once it stabilizes it will go down again.
					//break;
				}


			} // end of line search
#endif

			// 		   								LINE SEARCH										//
			//																						//
			// =====================================================================================//

			Numerics->cumCorrection_fac += Numerics->lsBestGlob;
			Numerics->lsLastRes = EqStokes->normResidual;

/*
			if (Numerics->lsState == -2) {
				//printf("Break!!\n");
				break;
			}



			// anti-Numerics->stalling
			if (fabs(EqStokes->normResidual-Numerics->oldRes)<EqStokes->normResidual*Numerics->relativeTolerance) {
				break;
				Numerics->stalling = true;
				Numerics->stallingCounter++;
			} else {
				Numerics->stalling = false;
				Numerics->stallingCounter = 0;
			}
			*/
			/*
			if (fabs(EqStokes->normResidual-Numerics->oldRes)<Numerics->absoluteTolerance*1e-4) {
				break;
			}
			*/
			



#if NON_LINEAR_VISU
		Visu->update = true;
		Visu->updateGrid = false;
		Visu_main(&Model);
		if (glfwWindowShouldClose(Visu->window))
			break;
#endif


#if (!NON_LINEAR_VISC)
/*
	compute dtAdv 	= Numerics->CFL_fac_Stokes*Grid->dx/(Physics->maxVx); // note: the min(dx,dy) is the char length, so = 1
	dtAdv 	= fmin(dtAdv,  Numerics->CFL_fac_Stokes*Grid->dy/(Physics->maxVy));
	printf("dtAdv = %.2e, dt = %.2e, dt = %.2e yrs, lsGlob = %.2e\n", dtAdv, Physics->dt,  Physics->dt * (Char->time/(3600.0*24*365)), Numerics->lsGlob);
	printf("dt = %.2e yrs\n",Physics->dt * (Char->time/(3600.0*24*365)));
	if (dtAdv<Physics->dt && Physics->dt>Numerics->dtMin) {
		Physics_dt_update(&Model);
		if (Physics->dt!=Physics->dt) {
			Numerics->oneMoreIt = true;
		} else {
			Numerics->oneMoreIt = false;
		}
	} else {
		Numerics->oneMoreIt = false;
	}
*/	
#endif


			


			if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
				printf("\n\n\n\nerror: Something went wrong. The norm of the residual is NaN\n");
				break;
			}

			if (EqStokes->normResidual>1e10) {
				break;
			}



			//Numerics->oneMoreIt = false; // for some reasons it stalls sometime
#endif

			Numerics->itNonLin++;
		} // end of non-linear loop


#if (VISU)
		Visu->nonLinItisOver = true;
#endif


		if (isnan(EqStokes->normResidual) || isinf(EqStokes->normResidual)) {
			printf("\n\n\n\nerror: Something went wrong. The norm of the residual is NaN\n");
			break;
		}
		if (EqStokes->normResidual>1e10) {
			break;
		}



#if (VISU)
		double timeStepToc = omp_get_wtime();;
		toc = (timeStepToc-timeStepTic);
		printf("the timestep took: %.2f\n",toc);
#endif



		// 										NON-LINEAR ITERATION 											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//



#if (VISU)
		timeStepTic = omp_get_wtime();;
#endif

#if (VISCOSITY_TYPE==0)

		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 									INTERPOLATION FROM CELL TO PARTICLES								//

		// update stress on the particles
		// =============================
		Physics_Dsigma_updateGlobal  (&Model);
		//Physics_Sigma0_updateGlobal_fromGrid(&Model);
#if (ADV_INTERP)
		Interp_Stresses_Grid2Particles_Global(&Model);
		//Physics_Eta_computeLambda_FromParticles_updateGlobal(&Model, true);
#endif

#if (DARCY)
		Interp_Phi_Grid2Particles_Global	(&Model);
#endif


#if (HEAT)
		for (i = 0; i < Grid->nECTot; ++i) {
			Physics->DT[i] = Physics->T[i] - Physics->T0[i];
		}
		Interp_Temperature_Grid2Particles_Global(&Model);

#endif
#if (STORE_PLASTIC_STRAIN)
		Interp_Strain_Grid2Particles_Global(&Model);
#endif

		// 									INTERPOLATION FROM CELL TO PARTICLES								//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//

#endif




		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 												OUTPUT AND VISU											//

		// Output
		// =================
		Output_call(&Model);
		


#if VISU
		if (Numerics->timeStep==0 || ( !Visu->useTimeFrequency && Visu->stepsSinceLastRender==Visu->renderFrequency  ) || (  Visu->useTimeFrequency &&  Physics->time>Visu->renderCounter*Visu->renderTimeFrequency)  ) {
			// Render
			//if (Visu->renderFrequency)
			Visu->update = true;
			if (!Grid->isFixed) {
				Visu->updateGrid = true;
			}
			Visu_main(&Model);
			Visu->stepsSinceLastRender = 1;
			if (Numerics->timeStep>0) {
    			Visu->timeSinceLastRender -= Visu->renderTimeFrequency;
			}
			Visu->timeSinceLastRender += Physics->dtAdv;
			Visu->renderCounter += 1;
		} else {
			Visu->stepsSinceLastRender += 1;
    		Visu->timeSinceLastRender += Physics->dtAdv;
		}


		glfwPollEvents();
		if (glfwWindowShouldClose(Visu->window))
			break;
		

#endif


		// 												OUTPUT AND VISU											//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//











		//======================================================================================================//
		// =====================================================================================================//
		//																										//
		// 							ADVECTION AND INTERPOLATION	FROM PARTICLES TO CELL							//

#if (ADV_INTERP)	
		// Advect Particles
		// =============================
		printf("Particles: Advect\n");
		Particles_advect(Particles, Grid, Physics);

		// Update the linked list of particles
		// =================================
		printf("Particles Update Linked List\n");
		Particles_updateLinkedList(Particles, Grid, Physics);


		// Inject particles
		// =================================
		if (Grid->isFixed) {
			Particles_injectAtTheBoundaries(Particles, Grid, Physics, MatProps);
		}

		printf("Particle injection\n");
		Particles_injectOrDelete(Particles, Grid);



		// Advect the box and update Particles position if needed
		// =============================
		switch (BCStokes->SetupType) {
		case Stokes_PureShear:
		case Stokes_Sandbox:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				Grid_updatePureShear(&Model);
				Particles_teleportInsideTheDomain(Particles, Grid, Physics);
			}
			break;
		case Stokes_SimpleShear:
			Particles_Periodicize(Particles, Grid);
			break;
		case Stokes_FixedLeftWall:
			break;
		case Stokes_CornerFlow:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				printf("error: For the corner flow setup, Grid->fixedBox must be true. Correct the input file");
				exit(0);
			}
			break;
		case Stokes_WindTunnel:
			if (Grid->isFixed) {
				Particles_deleteIfOutsideTheDomain(Particles, Grid);
			} else {
				printf("error: For the WindTunnel setup, Grid->fixedBox must be true. Correct the input file");
				exit(0);
			}
			break;
		default:
			break;
		}



		Particles_switchStickyAir			(Particles, Grid, Physics, Numerics, MatProps, BCStokes);
		Particles_surfaceProcesses			(&Model);
		// Update the Phase matrix
		// =================================
		Physics_Phase_updateGlobal					(&Model);

#if (VISCOSITY_TYPE==0)
		// Update the Physics on the Cells
		// =================================
		printf("Physics: Interp from particles to grid\n");
		Interp_All_Particles2Grid_Global(&Model);
#endif
#if (INERTIA)
		if (Numerics->timeStep>0) {
			Physics_Velocity_advectEulerian(&Model);
		} else {
			#if (INERTIA)
			Physics_VelOld_POld_updateGlobal(&Model);
			#endif
		}
#endif

		Physics_Rho_updateGlobal(&Model);	
		



#if (DARCY)
		compute dx, dy;
		for (iy = 1; iy < Grid->nyEC-1; ++iy) {
			for (ix = 1; ix < Grid->nxEC-1; ++ix) {
				iCell = ix + iy*Grid->nxEC;
				dx = Grid->DXS[ix-1];
				dy = Grid->DYS[iy-1];
				Physics->divV0[iCell]  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/dx;
				Physics->divV0[iCell] += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/dy;
			}
		}
#endif

		// Update BC
		// =================================
		printf("BC: Update\n");
		BCStokes->counter = 0;
		BC_updateStokes_Vel(BCStokes, Grid, Physics, true);
#if (DARCY)
		BC_updateStokesDarcy_P(BCStokes, Grid, Physics, true);
#endif
#if (HEAT)
		BCThermal->counter = 0;
		BC_updateThermal(BCThermal, Grid, Physics, true);
#endif

#else // i.e. if (!ADV_INTERP)
		Physics_Sigma0_updateGlobal_fromGrid(&Model);
#endif // if (ADV_INTERP)



		// 							ADVECTION AND INTERPOLATION FROM PARTICLES TO CELL 							//
		//																										//
		//======================================================================================================//
		// =====================================================================================================//
		printf("timeN = %.2e\n", Physics->time*Model.Char.time);

		Physics->time += Physics->dtAdv;
		Physics->dtAdv0 = Physics->dtAdv;
		Numerics->timeStep++;

		Physics_dt_update(&Model);
		Physics_Eta_Simple_updateGlobal(&Model);

#if (VISU)
		timeStepToc = omp_get_wtime();;
		toc = (timeStepToc-timeStepTic);
		printf("interp+adv+visu timestep took: %.2f\n",toc);
#endif


		double globalToc = omp_get_wtime();;
		toc = (globalToc-globalTic);
		toc += (3600.0*24.0)*2.0 + 10.0*3600 + 35.0*60.0;
		//printf("globalToc = %.6e, globalTic = %.6e, tickDur = %.2e, toc = %.2e\n", globalToc, globalTic, tickDuration, toc);
		double daySpent = floor(toc/(3600.0*24.0));
		toc -= daySpent*(3600.0*24.0);
		double hourSpent = floor(toc/3600.0);
		toc -= hourSpent*3600.0;
		double minSpent = floor(toc/60.0);
		toc -= minSpent*60.0;
		double secSpent = toc;
		printf("time since last restart: %.0f d, %.0f h, %.0f m, %.0f s\n", daySpent, hourSpent, minSpent, secSpent );
		Breakpoint->realTimeFrequency = 30.0;
		if (Breakpoint->realTimeFrequency>0) {
			printf("Breakpoint: write data");
			if ((globalToc-globalTic)>Breakpoint->counter*Breakpoint->realTimeFrequency) {
				Breakpoint_writeData(&Model);
				Breakpoint->counter ++;
			}
		} else {
			if (Breakpoint->frequency>0 && (Output->counter % Breakpoint->frequency)==0 && (Breakpoint->counter!=Output->counter)) {
				Breakpoint->counter = Output->counter;
				Breakpoint_writeData(&Model);
			}
		}
	}

	printf("Exiting\n");

//
//                          								END OF TIME LOOP													//
//																																//
//==============================================================================================================================//
//==============================================================================================================================//







	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                                    EXIT          	                      //

	free(NonLin_x0);
	free(NonLin_dx);
	// Free memory
	printf("Free Physics->..\n");
	Physics_Memory_free(&Model);
	printf("Free NumStokes->..\n");
	Numbering_Memory_free(NumStokes);

	printf("Free EqStokes->..\n");
	EqSystem_Memory_free(EqStokes, SolverStokes);
	printf("Free BCStokes->..\n");
	BC_Memory_free(BCStokes);
#if (HEAT)
	printf("Free NumThermal->..\n");
	Numbering_Memory_free(NumThermal);
	printf("Free EqThermal->..\n");
	EqSystem_Memory_free(EqThermal,SolverThermal);
	printf("Free BCThermal->..\n");
	BC_Memory_free(BCThermal);
#endif
	printf("Free Particles->..\n");
	Particles_Memory_free(Particles, Grid);
	printf("Free Numerics->..\n");
	Numerics_Memory_free(Numerics);
	printf("Free Grid->..\n");
	Grid_Memory_free(Grid);
	printf("Free Output->..\n");
	Output_free(Output);


#if VISU
	// Quit glfw
	printf("Quit GFLW...\n");
	glfwDestroyWindow(Visu->window);
	glfwTerminate();
	printf("Free Visu->..\n");
	Visu_Memory_free(Visu);
#endif

printf("Memory freed successfully\n");

	//return EXIT_SUCCESS;
	return 1;

	//                                    EXIT          	                      //
	//                                                                            //
	//============================================================================//
	//============================================================================//

}




