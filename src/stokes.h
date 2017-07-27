/*
 * stokes.h

 *
 *  Created on: Feb 9, 2016
 *      Author: abauville
 */


//============================================================================//
//============================================================================//
//                                                                            //
//                                  INCLUDES                              	  //
//                                                                            //
//============================================================================//
//============================================================================//
#ifndef STOKES_H_
#define STOKES_H_


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>



#define DEBUG   false
#define VISU 	true
#define HEAT	false


#if (VISU)
#define NON_LINEAR_VISU false
#else
#define NON_LINEAR_VISU false
#endif
#define MULTI_VISU false


#define VISCOSITY_TYPE 0 // 0: non linear, 1: linear viscous, 2: homogeneously constant (single phase)

#define DARCY false

#define STORE_PARTICLE_POS_INI false

#define STRAIN_SOFTENING false


#define INPUT_FILE "./Setups/input.json"

#define FREE_SURFACE_STABILIZATION false

#define CRANK_NICHOLSON_VEL true
#if (CRANK_NICHOLSON_VEL)
#define CRANK_NICHOLSON_P false // BROKEN
#else
#define CRANK_NICHOLSON_P false
#endif



#define INERTIA false

// Advection of velocity and viscosity fields
#define ADVECT_VEL_AND_VISCOSITY false
#define VEL_VISC_METHOD 0
#define ADVECT_METHOD 1 // 0, from Vx, Vy nodes to particles; 1, Vx,Vy interpolated on cell centers, then, interpolated to particles




#if (VISU)
//#ifdef __APPLE__
#include <GL/glew.h>
//#endif
#include <GLFW/glfw3.h>
#include <png.h>

#endif

#include <math.h>

#include <string.h>
#include <omp.h>

#include <time.h>




#define OMP_SCHEDULE schedule(static,32)






//============================================================================//
//============================================================================//
//                                                                            //
//                                   MACROS                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
#define NB_PHASE_MAX 10
#define NB_SUBSYSTEM_MAX 8
#define UPPER_TRI true
#define TIMER false

//#define INIT_TIMER 	clock_t tic, diff;
//					float toc;
//#define TIC 	tic = clock();
//#define TOC 	diff = (clock() - tic);
//				toc = (float) (diff * 1000 /CLOCKS_PER_SEC)/1000;

#define INIT_TIMER 	double tic; \
					double toc;

#if (VISU)
#define TIC tic = glfwGetTime();
#define TOC toc = glfwGetTime(); \
            toc = toc - tic;
#else
#define TIC tic = 0;
#define TOC toc = 0; \
            toc = toc - tic;
#endif


#define PI 3.14159265358979323846


#define INIT_PARTICLE SingleParticle* thisParticle = NULL; \
						int iNode = 0;

#define FOR_PARTICLES  	  for (iNode = 0; iNode < Grid->nSTot; ++iNode) { \
								thisParticle = Particles->linkHead[iNode]; \
								while (thisParticle != NULL) {

#define END_PARTICLES  			thisParticle = thisParticle->next; \
								} \
							}


#define MAX_STRING_LENGTH 2048

#define MAX_VISU_TYPE 32

//============================================================================//
//============================================================================//
//                                                                            //
//                                 TYPE DEFS                              	  //
//                                                                            //
//============================================================================//
//============================================================================//

// Physics
// =========================
// Definitions in Physics.h
typedef struct SinglePhase SinglePhase;
typedef struct Physics Physics;

// Particles
// =========================
// Definitions in Particles.h
typedef struct SingleParticle SingleParticle;
typedef struct ParticlePointerList ParticlePointerList;
typedef struct Particles Particles;

// Visu
// =========================
// Definitions in Visu.h
typedef struct ColorMap ColorMap;
typedef struct Visu Visu;




// Define basic data types
// =========================
typedef double coord;
typedef double compute;

// Input
// =========================
typedef struct Input 
{
	char inputFile[MAX_STRING_LENGTH];
	char currentFolder[MAX_STRING_LENGTH];
} Input;


// Output
// =========================

typedef enum {OutFormat_Float, OutFormat_Double} OutFormat;
typedef enum {Out_Vx, Out_Vy, Out_P, Out_Pf, Out_Pc, Out_Viscosity, Out_Porosity, Out_Z, Out_G, Out_Khi, Out_Sxx0, Out_Sxy0, Out_Sxx, Out_Sxy, Out_SII, Out_StrainRate, Out_Temperature, Out_Phase} OutType;
typedef enum {OutPart_x, OutPart_y, OutPart_xIni, OutPart_yIni, OutPart_Phase, OutPart_Passive, OutPart_T, OutPart_DeltaP0, OutPart_Sxx0, OutPart_Sxy0, OutPart_Phi} OutPartType;
typedef struct Output {
	char outputFolder[MAX_STRING_LENGTH];
	//OutFormat OutputFormat;
	OutType type[14];
	int nTypes, nPartTypes; // number of data matrices outputted at each time step
	int counter;
	int frequency;
	compute timeFrequency;
	bool useTimeFrequency;
	bool saveFirstStep;
	//char *ModelDescription;

	OutPartType partType[11];



} Output;


// Numerics
// =========================
typedef struct Numerics 
{

	int timeStep;

	int nTimeSteps; //  negative value for infinite

	compute maxTime;

	int nLineSearch;
	int maxNonLinearIter; // should always be greater than the number of line searches
	int minNonLinearIter; // should always be greater than the number of line searches
	compute relativeTolerance; // relative tolerance to the one of this time step
	compute absoluteTolerance; // relative tolerance to the first one of the simulation
	compute maxCorrection;

	int itNonLin;


	compute CFL_fac_Stokes, CFL_fac_Darcy, CFL_fac_Thermal;
	bool use_dtMaxwellLimit;
	compute dtMin, dtMax;
	compute etaMin, etaMax;
	compute dLmin; // min grid size

	compute minRes;

	compute normRes0, normResRef;

	compute cumCorrection_fac; 	// cumulative correction factor = sum of globalization
								// for a given non linear iteration, should be at 1.0 before to pass to the next time step

	compute phiMin, phiMax, phiCrit;

	int lsState;
	int lsCounterUp, lsCounter;
	compute lsGlob, lsGlobStart; // globalization factor
	compute lsBestGlob, lsBestRes;
	compute lsLowerBound, lsUpperBound;
	compute lsLastRes;
	compute lsResTolImprovement;
	compute lsGlobMin;
	compute lsLastGlob;
	compute lsTolDiverge;
	compute lsbestRes_It;

	compute StickyAirStress;

	compute stickyAirSwitchingDepth, stickyAirTimeSwitchPassive, stickyAirTimeSinceLastPassiveSwitch;
	int stickyAirSwitchPhaseTo, stickyAirSwitchPassiveTo;

	compute dtCorr, dtPrevCorr, dtAlphaCorr, dtAlphaCorrIni;

	compute dtVep;
	compute dtMaxwellFac_EP_ov_E;  // lowest
	compute dtMaxwellFac_VP_ov_E;  // intermediate
	compute dtMaxwellFac_VP_ov_EP; // highest

	compute dtPrevTimeStep;

	bool oneMoreIt; // used by time step size routine: if the time step is not fixed then continue iterating
	bool lsGoingDown, lsGoingUp;
} Numerics;



// Characteristic physical quantities
// =========================
typedef struct Char
{
	compute time;    		// [s]
	compute length;  		// [m]
	compute mass;  			// [kg]
	compute velocity; 		// [m.s-1]
	compute density; 		// [kg.m^-3]
	compute stress;  		// [Pa] or [kg.m^-1.s^-2]
	compute viscosity; 		// [Pa.s] or [kg.m^-1.s-1]
	compute acceleration; 	// [m.s^-2]
	compute strainrate; 	// [s^-1]
	compute temperature; 	// [K]

} Char;







// Grid
// =========================
typedef struct Grid 
{
	int nxC, nyC, nCTot; 			// number of cells / cell center nodes
	int nxEC, nyEC, nECTot; 		// number of embedded cells = cells + ghost cells around, useful for interpolation
	int nxS, nyS, nSTot; 			// number of base nodes
	int nxVx, nyVx, nVxTot; 		// number of Vx nodes
	int nxVy, nyVy, nVyTot; 		// number of Vy nodes
	coord xmin, xmax, ymin, ymax; 	// grid extent
	coord xmin_ini, xmax_ini, ymin_ini, ymax_ini; 	// grid extent
	compute dx, dy; 		 	  	// grid increment / cell size
	compute *X, *Y, *DXS, *DYS, *DXEC, *DYEC;
	int nSegX, *nxPerSeg;
	int nSegY, *nyPerSeg;
	compute *xSeg, *ySeg;

	bool userDefined;
	bool isFixed;
	bool isPeriodic;
} Grid;



// Material properties
// =========================
typedef enum {TensorCorrection_None, TensorCorrection_UniAxial,TensorCorrection_SimpleShear} TensorCorrection;

typedef struct DiffCreepProps 
{
	bool isActive;
	//TensorCorrection tensorCorrection;
	compute B, E, V;
	//compute d0, p, C_OH_0, r;
} DiffCreepProps;
typedef struct DislCreepProps 
{
	bool isActive, MPa;
	TensorCorrection tensorCorrection;
	compute B, E, V, n;
	//compute n, C_OH_0, r;
} DislCreepProps;
typedef struct PeiCreepProps 
{
	bool isActive;
	compute B, E, V;
	compute tau, gamma, q;
} PeiCreepProps;



typedef struct MatProps 
{
	int nPhase;
	char name[NB_PHASE_MAX][128];

	compute rho0[NB_PHASE_MAX];
	compute alpha[NB_PHASE_MAX]; // thermal expansion
	compute beta[NB_PHASE_MAX];  // compressibility
	compute k[NB_PHASE_MAX]; 	 // thermal conductivity
	compute G[NB_PHASE_MAX]; 	 // shear modulus

	compute cohesion[NB_PHASE_MAX]; // cohesion
	compute frictionAngle[NB_PHASE_MAX]; // angle of friction


	compute perm0[NB_PHASE_MAX];
	compute perm0_eta_f[NB_PHASE_MAX]; //perm0/eta_f
	//compute eta_b[NB_PHASE_MAX];
	//compute B[NB_PHASE_MAX];

	bool isAir[NB_PHASE_MAX];
	bool isWater[NB_PHASE_MAX];
	bool use_dtMaxwellLimit[NB_PHASE_MAX];

	DiffCreepProps vDiff[NB_PHASE_MAX];
	DislCreepProps vDisl[NB_PHASE_MAX];
	PeiCreepProps  vPei [NB_PHASE_MAX];

	compute phiIni[NB_PHASE_MAX];

} MatProps;






// Boundary conditions
// ========================
typedef enum {Dirichlet, DirichletGhost, Neumann, NeumannGhost, Infinity} BCType;
typedef enum {Stokes_PureShear, Stokes_SimpleShear, Stokes_FixedLeftWall, Stokes_Sandbox, Stokes_SandboxWeakBackstop, Stokes_CornerFlow, Stokes_WindTunnel,
			  Thermal_TT_TB_LRNoFlux, Thermal_TT_TBExternal_LRNoFlux,
			  Darcy_Default} SetupType;
typedef struct BC
{
	int n;

	int nPerEqSubSys[10]; // number of BC for each subsystem of equations

	int *list;
	BCType *type;
	compute *value;
	int counter;


	bool hPeriod;
	bool mobileBox;

	//int typeL, typeR, typeT, typeB;
	SetupType SetupType;
	// Should be moved somewhere else, but I don't know where yet
	compute backStrainRate; // background strain, positive in extension
	compute refValue; // ref Value, e.g. plate velocity U for the corner flow, or a reference Temp, or reference Pressure
	compute TB, TT;

	compute specialPhase;

	compute DeltaL; // For infinity like BC


	// Temporary: chotto dirty
	compute Sandbox_TopSeg00;
	compute Sandbox_TopSeg01;
	bool Sandbox_NoSlipWall;

	compute Corner_SubductionAngle;

	bool IsFreeSlipLeft, IsFreeSlipRight, IsFreeSlipBot, IsFreeSlipTop; // Free slip info, used to enforce that Sigma_xy is 0 on the boundary
} BC;




// Initial conditions
// ========================
typedef enum {IC_HSC, IC_Gaussian} ICSetupType;
typedef struct IC 
{
	ICSetupType SetupType;
	compute data[32];
} IC;





// Equation System
// ========================
typedef struct EqSystem 
{
	//bool penaltyMethod;
	//compute penaltyFac;
	int nEqIni, nEq, nRow, nnz;

	// Stiffness matrix
	int *I; // Allocated in function Numbering_initMapAndSparseTripletIJ
	int *J;
	compute *V;

	compute *b; // right hand side
	compute *x; // solution vector;

	compute *S; // Scaling diagonal matrix (stored as a vector)

	compute normResidual;
	compute norm_b;
} EqSystem;




// Numbering
// ========================
typedef enum {Stencil_Stokes_Momentum_x, Stencil_Stokes_Momentum_y, Stencil_Stokes_Continuity, Stencil_Heat, Stencil_Stokes_Darcy_Momentum_x, Stencil_Stokes_Darcy_Momentum_y, Stencil_Stokes_Darcy_Continuity, Stencil_Stokes_Darcy_Darcy, Stencil_Poisson} StencilType;
typedef struct Numbering 
{
	int* map;
	int *IX, *IY;

	int subEqSystem0[NB_SUBSYSTEM_MAX]; // Index of the first equation of a subset (e.g. first Vx, first Vy, first P equation)// hard coded number. Assuming there will never be more than 10 subsystem of equations to solve
	int subEqSystem0Dir[NB_SUBSYSTEM_MAX];
	int nSubEqSystem;
	StencilType Stencil[NB_SUBSYSTEM_MAX];
} Numbering;




// Local Numbering Vx
// ========================
typedef struct LocalNumberingVx 
{
	int VxC, VxN, VxS, VxE, VxW;
	int VySW, VyNW, VySE, VyNE;
	int PW, PE;
	int NormalE, NormalW, ShearN, ShearS;
} LocalNumberingVx;

// Local Numbering Vy
// ========================
typedef struct LocalNumberingVy 
{
	int VxSW, VxSE, VxNW, VxNE;
	int VyS, VyW, VyC, VyE, VyN;
	int PS, PN;
	int NormalN, NormalS, ShearE, ShearW;
} LocalNumberingVy;

// Local Numbering P
// ========================
typedef struct LocalNumberingP 
{
	int VxE, VxW, VyN, VyS;
} LocalNumberingP;



// A basic linked list node struct
// ========================
typedef struct LinkedNode 
{
    int data;
    struct LinkedNode* next;
} LinkedNode;


// Pardiso solver
// ========================
typedef struct Solver 
{
	// Pardiso struct
	void    *pt[64];
	int      iparm[64];
	double   dparm[64];
	int      maxfct, mnum, msglvl, mtype, nrhs;
} Solver;




//============================================================================//
//============================================================================//
//                                                                            //
//                                PROTOTYPES                            	  //
//                                                                            //
//============================================================================//
//============================================================================//




// Char
// =========================
void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BCStokes, BC* BCThermal, IC* ICThermal, IC* ICDarcy, Numerics* Numerics, Particles* Particles, Output* Output);




// Grid
// =========================
void Grid_Memory_allocate	(Grid* Grid);
void Grid_Memory_free		(Grid* Grid);
void Grid_init				(Grid* Grid, Numerics* Numerics);
void Grid_updatePureShear	(Grid* Grid, BC* BC, Numerics* Numerics, compute dt);





// Interp
// =========================
void Interp_All_Particles2Grid_Global			(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumStokes, Numbering* NumThermal, BC* BCThermal);
void Interp_Temperature_Grid2Particles_Global	(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps, BC* BCThermal);
void Interp_Stresses_Grid2Particles_Global		(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal, MatProps* MatProps, Numerics* Numerics);
void Interp_Phi_Grid2Particles_Global			(Grid* Grid, Particles* Particles, Physics* Physics);
void Interp_Strain_Grid2Particles_Global		(Grid* Grid, Particles* Particles, Physics* Physics);
extern compute Interp_ECVal_Cell2Node_Local(compute* A, int ix, int iy, int nxEC);
extern compute Interp_ECVal_Node2Cell_Local(compute* A, int ix, int iy, int nxS);
extern compute Interp_ECVal_Cell2Particle_Local(compute* A, int ix, int iy, int nxEC, compute locX, compute locY);
extern compute Interp_ECVal_Node2Particle_Local(compute* A, int ix, int iy, int nxEC, compute locX, compute locY, int signX, int signY);




// Boundary conditions
// =========================
void BC_Memory_free			(BC* BC);
void BC_initStokes			(BC* BC, Grid* Grid, Physics* Physics, EqSystem* EqSystem);
void BC_initThermal			(BC* BC, Grid* Grid, Physics* Physics, EqSystem* EqSystem);
void BC_updateStokes_Vel	(BC* BC, Grid* Grid, Physics* Physics, bool assigning);
void BC_updateStokes_P		(BC* BC, Grid* Grid, Physics* Physics, bool assigning);
void BC_updateStokesDarcy_P	(BC* BC, Grid* Grid, Physics* Physics, bool assigning);
void BC_updateThermal		(BC* BC, Grid* Grid, Physics* Physics, bool assigning);


// Initial conditions
// =========================
#if (HEAT)
void IC_T(Physics* Physics, Grid* Grid, IC* ICThermal, BC* BCThermal);
#endif
#if (DARCY)
void IC_phi(Physics* Physics, Grid* Grid, Numerics* Numerics, IC* ICDarcy, MatProps* MatProps, Particles* Particles);
#endif



// Numbering
// =========================
void Numbering_Memory_allocate	(Numbering* Numbering, EqSystem* EqSystem, Grid* Grid);
void Numbering_Memory_free		(Numbering* Numbering);
//inline compute Interp_ECVal_Cell2Node_Local(compute* A, int ix, int iy, int nxEC) __attribute__((always_inline));
void Numbering_init				(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering, Physics* Physics, Numerics* Numerics);



// Equation system
// =========================
void EqSystem_Memory_allocateI		(EqSystem* EqSystem);
void EqSystem_Memory_allocate(EqSystem* EqSystem);
void EqSystem_Memory_free	(EqSystem* EqSystem, Solver* Solver) ;
void EqSystem_assemble		(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering, bool updateScale, Numerics* Numerics);
void EqSystem_solve			(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_check			(EqSystem* EqSystem);
void EqSystem_initSolver  	(EqSystem* EqSystem, Solver* Solver);
void pardisoSolveSymmetric	(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_computeNormResidual(EqSystem* EqSystem);
void EqSystem_scale			(EqSystem* EqSystem);
void EqSystem_unscale		(EqSystem* EqSystem);



// Local stencil
// =========================
void LocalStencil_Call(StencilType Stencil, int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic, Numerics* Numerics);
void LocalStencil_Stokes_Momentum_x	(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic, Numerics* Numerics);
void LocalStencil_Stokes_Momentum_y	(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic, Numerics* Numerics);
void LocalStencil_Stokes_Continuity	(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
#if (HEAT)
void LocalStencil_Heat				(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
#endif
#if (DARCY)
void LocalStencil_Stokes_Darcy_Momentum_x(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
void LocalStencil_Stokes_Darcy_Momentum_y(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
void LocalStencil_Stokes_Darcy_Continuity(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
void LocalStencil_Stokes_Darcy_Darcy 	 (int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
#endif


// PARDISO
// =========================
void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                  int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);




// Utils
// =========================
#if (VISU)
	void compileShaders(GLuint *ShaderProgram, const char* pVSFileName, const char* pFSFileName, const char* pGSFileName, bool useGS);
#endif
char* readFile		(const char* fileName);
void printListi		(int* list, int length);
void printListd		(double* list, int length);
int  findi 			(int* list, int length, int value);
double min			(double* List, int length);
double max			(double* List, int length);
float  minf			(float* List, int length);
float  maxf			(float* List, int length);
double absmin		(double* List, int length);
double absmax		(double* List, int length);
void addToLinkedList		(LinkedNode** pointerToHead, int x);
void freeLinkedList			(LinkedNode* head);
#if (VISU)
	int writePNGImage	(char* filename, int width, int height, unsigned char *buffer, char* title);
#endif


// Darcy
// =========================
/*
typedef enum {Air, Ocean, Solid} PhaseFlag;
void Darcy_setBC		(Grid* Grid, Physics* Physics, coord hOcean, PhaseFlag* Phase);
void Darcy_setPhaseFlag	(PhaseFlag* Phase, coord hOcean, Grid* Grid, Particles* Particles);
void Darcy_solve		(Darcy* Darcy, Grid* Grid, Physics* Physics, MatProps* MatProps, Particles* Particles);
*/

// Numerics
// ========================
void Numerics_init		(Numerics* Numerics);
void Numerics_Memory_free(Numerics* Numerics);
//int  Numerics_updateBestGlob(Numerics* Numerics, EqSystem* EqStokes, int* iLS);
void Numerics_LineSearch_chooseGlob(Numerics* Numerics, EqSystem* EqStokes);

// Input
// ========================
void Input_read(Input* Input, Grid* Grid, Numerics* Numerics, Physics* Physics, MatProps* MatProps, Particles* Particles, Char* Char, BC* BCStokes, BC* BCThermal, IC* ICThermal, IC* ICDarcy, Output* Output);
void Input_assignPhaseToParticles(Input* Input, Particles* Particles, Grid* Grid, Char* Char);

#if (VISU)
void Input_readVisu(Input* Input, Visu* Visu);
#endif


// Output
// ========================
void Output_free						(Output* Output);
void Output_writeInputCopyInOutput		(Output* Output, Input* Input);
void Output_modelState					(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics);
void Output_data 						(Output* Output, Grid* Grid, Physics* Physics, Char* Char, Numerics* Numerics);
void Output_particles					(Output* Output, Particles* Particles, Grid* Grid, Char* Char, Numerics* Numerics);


// Physics
// =========================
#include "Physics.h"

// Particles
// =========================
#include "Particles.h"

// Visu
// =========================
#if (VISU)
#include "Visu.h"
#endif



/*
// Mikito's bitmap reader
#define HEADERSIZE   54
#define PALLETSIZE 1024
#define MAXWIDTH   1000
#define MAXHEIGHT  1000


#define SWAP(x,y) {typeof(x) temp; temp=x; x=y; y=temp;}

unsigned char Bmp_headbuf[HEADERSIZE];
unsigned char Bmp_Pallet[PALLETSIZE];

char Bmp_type[2];
unsigned long Bmp_size;
unsigned int Bmp_info_header_size;
unsigned int Bmp_header_size;
long Bmp_height;
long Bmp_width;
unsigned short Bmp_planes;
unsigned short Bmp_color;
long Bmp_comp;
long Bmp_image_size;
long Bmp_xppm;
long Bmp_yppm;

typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} color;

typedef struct {
  long height;
  long width;
  color data[MAXHEIGHT][MAXWIDTH];
} img;

void ReadBmp(char *filename, img *imgp);
*/




#endif /* STOKES_H_ */




