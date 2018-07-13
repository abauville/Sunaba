/*
 * MainInit.c
 *
 *  Created on: 13 Jul 2018
 *      Author: abauville
 */
#include "stokes.h"

void Main_Init_start(Model* Model) {
    // Declare structures
	// =================================
	// General
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	Char* Char 				= &(Model->Char);
	// Stokes: Conservation of momentum + continuity (+- Darcy)
	Numbering* NumStokes 	= &(Model->NumStokes);
	BC* BCStokes 			= &(Model->BCStokes);
	EqSystem* EqStokes				= &(Model->EqStokes);
	Solver* SolverStokes 	= &(Model->SolverStokes);
#if (DARCY)
	IC* ICDarcy 			= &(Model->ICDarcy);
#endif
	// Heat conservation
	Numbering* NumThermal 	= &(Model->NumThermal);
#if (HEAT)
	IC* ICThermal 			= &(Model->ICThermal);
	BC* BCThermal 			= &(Model->BCThermal);
#endif
	EqSystem* EqThermal  	= &(Model->EqThermal);
#if (HEAT)
	Solver* SolverThermal 	= &(Model->SolverThermal);
#endif
	// Numerics
	Numerics* Numerics 		= &(Model->Numerics);
	// Visu
#if (VISU)
	Visu* Visu 				= &(Model->Visu);
#endif
	// Input/Output
	Input* Input 			= &(Model->Input);
	Output* Output 			= &(Model->Output);


    
	printf("Reading input\n");
	Input_read(Model);
#if (LINEAR_VISCOUS)
	if (Numerics->maxNonLinearIter>1) {
		printf("error: you requested %i non linear iterations, however, they are switched off due to LINEAR_VISCOUS==true\n", Numerics->maxNonLinearIter);
		exit(0);
	}
#endif

#if (!HEAT)
	if (Model->BCThermal.TB!=1.0 || Model->BCThermal.TT!=1.0) {
		printf("TB = %.3e, TT = %.3e\n",Model->BCThermal.TB, Model->BCThermal.TT);
		printf("error: you specified non default thermal boundary conditions, however, the heat equation is switched off due to HEAT==false ");
		exit(0);
	}
#endif


#if (VISU)
	printf("Reading Visu input\n");
	Input_readVisu(Model);
#endif

	printf("Reading input over\n");

	if (Output->nTypes>0) {
		Output_writeInputCopyInOutput(Output, Input);
	}

	Output_init(Output);
	
	//============================================================================//
	//============================================================================//
	//                                                                            //
	//                       		 USER INPUT      		   	              	  //
	//                                                                            //
	//============================================================================//
	//============================================================================//

	printf("nTimesteps = %i\n",Numerics->nTimeSteps);
	Grid->isPeriodic = false;
	if (BCStokes->SetupType==Stokes_SimpleShear) {
		Grid->isPeriodic = true;
		Grid->isFixed 	= true;
	}

	// Non-dimensionalization
	// =================================
	Char_nonDimensionalize(Model);

	
	


	if (DEBUG) {
		if (Grid->nCTot>200) {
			printf("error: The system size exceeds the maximum allowed for debugging\n");
			exit(0);
		}
	}

	if (Grid->isPeriodic && Grid->nxC%2!=0) {
		printf("error: When using the periodic boundaries nxC must be even (because of node 'coloring'\n");
		exit(0);
	}


//======================================================================================================
//======================================================================================================
//
//                          				INITIALIZATION
//
	Numerics->timeStep = -1;
	Numerics->itNonLin = -1;

	

	//Init Grid
	// =================================
	Grid_Memory_allocate(Grid);
	Grid_init(Model);

NumThermal->nSubEqSystem 	= 1;
	NumThermal->Stencil[0] = Stencil_Heat;
	EqThermal->nEqIni 		= Grid->nECTot;

#if (DARCY)
	NumStokes->nSubEqSystem 	= 4;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Darcy_Momentum_x; 	// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Darcy_Momentum_y; 	// Vy
	NumStokes->Stencil[2] 	= Stencil_Stokes_Darcy_Darcy;	   	// Pf
	NumStokes->Stencil[3] 	= Stencil_Stokes_Darcy_Continuity; 	// Pc
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot + Grid->nECTot + Grid->nECTot;
#else
#if (PENALTY_METHOD)
	NumStokes->nSubEqSystem 	= 2;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Momentum_x;		// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Momentum_y; 		// Vy
	//NumStokes->Stencil[2]	= Stencil_Stokes_Continuity;		// P
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot;// + Grid->nECTot;
#else
	NumStokes->nSubEqSystem 	= 3;
	NumStokes->Stencil[0] 	= Stencil_Stokes_Momentum_x;		// Vx
	NumStokes->Stencil[1] 	= Stencil_Stokes_Momentum_y; 		// Vy
	NumStokes->Stencil[2]	= Stencil_Stokes_Continuity;		// P
	EqStokes->nEqIni  	 	= Grid->nVxTot + Grid->nVyTot + Grid->nECTot;
#endif
#endif


	Numerics->oneMoreIt = false;

	// Init Physics
	// =================================
	printf("Init Physics: allocate memory\n");
	Physics_Memory_allocate	(Model);

	// Init Numerics
	// =================================
	printf("Init Numerics\n");
	Numerics_init(Numerics);

	// Initialize Particles
	// =================================
	printf("Particles: Init Particles\n");
	Particles_Memory_allocate	(Particles, Grid);
	Particles_initCoord			(Particles, Grid);
	Particles_updateLinkedList	(Particles, Grid, Physics); // in case a ridiculous amount of noise is put on the particle
	Input_assignPhaseToParticles(Model);
	Particles_initPassive		(Particles, Grid, Physics);

	// Initialize Physics
	// =================================
	printf("Physics: Init Physics\n");
#if (EXTRA_PART_FIELD)
	Input_setFieldsOnParticles(Model);
#endif
	Interp_All_Particles2Grid_Global	(Model);
	Physics_Rho_updateGlobal(Model);
	Physics_Phase_updateGlobal					(Model);
	
	// Set fields from file test
	


#if (HEAT)
	IC_T(Physics, Grid, ICThermal, BCThermal);
	Interp_All_Particles2Grid_Global	(Model);
#endif
#if (DARCY)
	IC_phi(Physics, Grid, Numerics, ICDarcy, MatProps, Particles);
	Interp_All_Particles2Grid_Global	(Model);
	memcpy(Physics->phi, Physics->phi0, Grid->nECTot * sizeof(compute));
#endif

	Physics_Rho_updateGlobal	(Model);

	Physics_P_initToLithostatic (Model);

	Physics_Eta_init(Model);
	//Physics_dt_update(Model);
	//Physics_Eta_init(Model);


#if (DEBUG)
	Physics_check(Model);
#endif


	// Set boundary conditions
	// =================================
	printf("BC: Set\n");
	BC_initStokes			(BCStokes , Grid, Physics, EqStokes);
#if (HEAT)
	BC_initThermal			(BCThermal, Grid, Physics, EqThermal);
#endif
	// Initialize Numbering maps without dirichlet and EqStokes->I
	// =================================
	printf("Numbering: init Stokes\n");
	EqSystem_Memory_allocateI		(EqStokes);
	printf("a\n");
	Numbering_Memory_allocate(NumStokes, EqStokes, Grid);
	printf("a\n");
	Numbering_init			(BCStokes, Grid, EqStokes, NumStokes, Physics, Numerics);
	printf("EqSystem: init Stokes\n");
	EqSystem_Memory_allocate	(EqStokes );

	printf("Number of Unknowns for Stokes: %i \n", EqStokes->nEq);

#if (HEAT)
	printf("Numbering: init Thermal\n");
	EqSystem_Memory_allocateI		(EqThermal);
	Numbering_Memory_allocate(NumThermal, EqThermal, Grid);
	Numbering_init			(BCThermal, Grid, EqThermal, NumThermal, Physics);
	printf("EqSystem: init Thermal\n");
	EqSystem_Memory_allocate	(EqThermal);

	printf("Number of Unknowns for Heat: %i \n", EqThermal->nEq);
#endif
	


#if (VISU)
	printf("Visu: Init Visu\n");

	Visu->ntri   	= 2;//Grid->nxC*Grid->nyC*2;//2;//Grid->nxC*Grid->nyC*2;
	Visu->ntrivert 	= Visu->ntri*3;
	Visu->nParticles = Particles->n+ (int) (Particles->n*0.1); // overallocate 5% of the number of particles

	Visu_Memory_allocate(Visu, Grid);
	Visu_init(Visu, Grid, Particles, Char, Input);
#endif

	printf("koko, nEq=%i, nVxTot = %i, nVyTot = %i, nECTot = %i\n", EqStokes->nEq, Grid->nVxTot, Grid->nVyTot, Grid->nECTot);
	// Init Solvers
	// =================================
	printf("EqStokes: Init Solver\n");
	EqSystem_assemble(EqStokes, Grid, BCStokes, Physics, NumStokes, false, Numerics); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (EqStokes, SolverStokes);

#if (HEAT)
	printf("EqThermal: Init Solver\n");
	EqSystem_assemble(EqThermal, Grid, BCThermal, Physics, NumThermal, false, Numerics); // dummy assembly to give the EqSystem initSolvers
	EqSystem_initSolver (EqThermal, SolverThermal);
#endif



//
//                          				INITIALIZATION
//
//======================================================================================================
//======================================================================================================




	// Update Cell Values with Part
	// =================================
	Interp_All_Particles2Grid_Global(Model);
	Physics_Rho_updateGlobal(Model);
	Physics_P_initToLithostatic 			(Model);
	
	// Update BC
	// =================================
	printf("BC: Update\n");
	BCStokes->counter = 0;
	BC_updateStokes_Vel(BCStokes, Grid, Physics, true);
#if (DARCY)
	BC_updateStokesDarcy_P(BCStokes, Grid, Physics, true);
#endif

#if (DARCY)
	Physics_Perm_updateGlobal(Model);
#endif
	Physics_Rho_updateGlobal(Model);

	


}

void Main_Init_restart(Model* Model) {

	// Do most of the Main_Init_start_stuff but not the particles

	// Read the input breakpoint files
	// then rewrite the Vx, Vy, P, Z

	// Write the particles

	// Interp_All_Particles2Grid_Global(Model);

	// Write scalar values like Physics->time, Physics->dt etc...

	// Update BC

}