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
#define HEAT  	false


#if (VISU)
#define NON_LINEAR_VISU true
#else
#define NON_LINEAR_VISU false
#endif
#define LINEAR_VISCOUS	false

#define DARCY false

#define STORE_PARTICLE_POS_INI false

#define STRAIN_SOFTENING true


#define INPUT_FILE "./Setups/input.json"

#define FREE_SURFACE_STABILIZATION false

#define CRANK_NICHOLSON_VEL true

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


#define MAX_STRING_LENGTH 1024

#define MAX_VISU_TYPE 32

//============================================================================//
//============================================================================//
//                                                                            //
//                                 TYPE DEFS                              	  //
//                                                                            //
//============================================================================//
//============================================================================//


// Define basic data types
// =========================
typedef double coord;
typedef double compute;

// Input
// =========================
typedef struct Input Input;
struct Input {
	char inputFile[MAX_STRING_LENGTH];
	char currentFolder[MAX_STRING_LENGTH];
};


// Output
// =========================
typedef struct Output Output;
typedef enum {OutFormat_Float, OutFormat_Double} OutFormat;
typedef enum {Out_Vx, Out_Vy, Out_P, Out_Pf, Out_Pc, Out_Viscosity, Out_Porosity, Out_Z, Out_G, Out_Khi, Out_Sxx0, Out_Sxy0, Out_Sxx, Out_Sxy, Out_SII, Out_StrainRate, Out_Temperature, Out_Phase} OutType;
typedef enum {OutPart_x, OutPart_y, OutPart_xIni, OutPart_yIni, OutPart_Phase, OutPart_Passive, OutPart_T, OutPart_DeltaP0, OutPart_Sxx0, OutPart_Sxy0, OutPart_Phi} OutPartType;

struct Output {
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



};


// Numerics
// =========================
typedef struct Numerics Numerics;
struct Numerics
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

};



// Characteristic physical quantities
// =========================
typedef struct Char Char;
struct Char
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

};





// Physics
// =========================
typedef struct SinglePhase SinglePhase;
struct SinglePhase {
	int phase;
	compute weight;
	SinglePhase* next;
};
typedef struct Physics Physics;
struct Physics
{
	compute R;
	// Physics Stokes
	compute g[2]; // gravity acceleration
	compute gFac[2]; // gravity acceleration / (norm gravity acceleration)
	compute dt, dtAdv, dtT, dtDarcy;
	compute *Vx, *Vy, *P;

#if (CRANK_NICHOLSON_VEL)
	compute *Vx0, *Vy0;
#endif

	compute maxVx, maxVy;
	compute *eta;

	//compute *eta0
	//compute *n;
	compute *rho_g;
	//compute *rho0_g; // Density*norm_g

#if (HEAT)

	compute *k;  // Thermal conductivity
	compute *T, *T0, *DT; // temperature stored on cell centers
#endif

	compute epsRef; // reference strainrate

	int *phase;

	//compute *Plitho;

#if (DARCY)


	compute *Pc, *DeltaP0, *DDeltaP; // old compaction pressure
	compute *phi, *Dphi, *phi0; // fluid phase fraction
	compute *Pf;

	compute* divV0;

	compute *perm0_eta_f, *perm_eta_f; // permeability/eta_f
	compute minPerm;
	compute *eta_b; // bulk viscosity
	//compute *B; // elastic bulk modulus

	compute eta_f, rho_f; // viscosity of the fluid
	compute rho_f_g; // rho_f*norm_g
	compute PfGrad_Air_X;
	compute PfGrad_Air_Y;

	compute y_oceanSurface;

	compute *khi_b; // sigmaII/(plastic multiplier), i.e. plastic viscosity

	compute *Zb;

#endif

	compute *khi, *khiShear; // sigmaII/(plastic multiplier), i.e. plastic viscosity

	compute *Z, *ZShear;

	compute *etaShear;

	// Stokes, elasticity related variables
	compute *sigma_xx_0, *sigma_xy_0; // old stresses
	compute *Dsigma_xx_0, *Dsigma_xy_0; // stress corrections for markers
	compute *G; // shear modulus

	// Plasticity
	compute *cohesion, *frictionAngle;

	//compute *etaVisc;


	// Physics thermal

	compute Cp; // heat capacity, taken as a single value because it varies very little between different types of rocks



	// Darcy

	compute dtMaxwellMin, dtMaxwellMax;


	compute time;


	int phaseAir;
	int phaseWater;
	int phaseRef;
	// compute stressOld

	SinglePhase **phaseListHead;
	compute *sumOfWeightsCells, *sumOfWeightsNodes;



#if (STRAIN_SOFTENING)
	compute* strain;
	compute* Dstrain;
#endif


};




// Grid
// =========================
typedef struct Grid Grid;
struct Grid
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
};



// Material properties
// =========================
typedef enum {TensorCorrection_None, TensorCorrection_UniAxial,TensorCorrection_SimpleShear} TensorCorrection;

typedef struct DiffCreepProps DiffCreepProps;
struct DiffCreepProps
{
	bool isActive;
	//TensorCorrection tensorCorrection;
	compute B, E, V;
	//compute d0, p, C_OH_0, r;
};
typedef struct DislCreepProps DislCreepProps;
struct DislCreepProps
{
	bool isActive, MPa;
	TensorCorrection tensorCorrection;
	compute B, E, V, n;
	//compute n, C_OH_0, r;
};
typedef struct PeiCreepProps PeiCreepProps;
struct PeiCreepProps
{
	bool isActive;
	compute B, E, V;
	compute tau, gamma, q;
};


typedef struct MatProps MatProps;
typedef enum {Custom} MaterialType;
struct MatProps
{
	int nPhase;

	MaterialType Material;

	char name[NB_PHASE_MAX][128];

	compute rho0[NB_PHASE_MAX];
	compute rho0_g[NB_PHASE_MAX]; // rho0 * norm_g
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

	DiffCreepProps vDiff[NB_PHASE_MAX];
	DislCreepProps vDisl[NB_PHASE_MAX];
	PeiCreepProps  vPei [NB_PHASE_MAX];

	compute phiIni[NB_PHASE_MAX];

};





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


#if (STRAIN_SOFTENING)
	compute strain;
#endif

	// for the linked list
	int nodeId;
    SingleParticle* next;

};

// Id Changed
typedef struct ParticlePointerList ParticlePointerList;
struct ParticlePointerList {
    //int data;
    SingleParticle* pointer;
    ParticlePointerList* next;
};
// Particles, i.e. info of the system of all particles
typedef struct Particles Particles;
typedef enum {PartPassive_Grid} ParticlePassiveGeom;
struct Particles
{
	int nPC, nPCX, nPCY; // number of particles per cell, tot, in x and in y
	int n; // number of particles
	compute minPartPerCellFactor, maxPartPerCellFactor;
	coord noiseFactor;
	SingleParticle **linkHead;

	ParticlePassiveGeom passiveGeom ;
	compute passiveDx, passiveDy;

	compute *dispAtBoundL, *dispAtBoundR;// ,*dispAtBoundT, *dispAtBoundB
	int *currentPassiveAtBoundL, *currentPassiveAtBoundR;

};




// Visualization
// ========================
#if (VISU)
typedef enum {Blank, Viscosity, StrainRate, Velocity, Pressure, Density, Temperature, Stress, FluidPressure, Permeability, Porosity, CompactionPressure, Phase,
			  VxRes, VyRes, PRes, PfRes, PcRes, TRes, VelocityDiv,SIIOvYield, PeOvYield, Khi, Khib, Strain} VisuType;
typedef enum {PartPhase, PartTemp,PartSigma_xx, PartSigma_xy, PartDeltaP, PartPorosity} ParticleVisuType;
typedef enum {StokesVelocity, DarcyGradient} GlyphType;
typedef enum {Triangle, ThinArrow, ThickArrow} GlyphMeshType;
typedef enum {Nearest, Linear} FilterType;

typedef struct ColorMap ColorMap;
struct ColorMap {
	//number     = number
    //type       = colormapType # "automatic would go from min to max values"
    compute colorMapRes;
    compute colorMap;
    compute scale;
    compute center; // centered value (scaled)
    compute max; // maximum value (scaled) from the center
    bool log10on;
};


typedef struct Visu Visu;
struct Visu
{

	GLFWwindow* window;

	int ntri, ntrivert;
	GLuint* elements;
	GLfloat* U;
	GLfloat* vertices;
	GLfloat scale, valueScale, valueShift;
	GLfloat colorScale[2], partColorScale[2];
	GLfloat shift[3], shiftFac[3];
	GLint log10_on;
	GLuint VAO, VBO, EBO;
	GLuint TEX;
	GLuint VAO_part, VBO_part, VBO_partMesh;
	GLuint VAO_glyph, VBO_glyph, VBO_glyphMesh;
	GLuint ShaderProgram, ParticleShaderProgram, ParticleBackgroundShaderProgram, GlyphShaderProgram;
	/*
	const char* VertexShaderFile;
	const char* FragmentShaderFile;

	const char* ParticleVertexShaderFile;
	const char* ParticleFragmentShaderFile;
	const char* ParticleGeometryShaderFile;
	const char* ParticleBackgroundVertexShaderFile;
	const char* ParticleBackgroundFragmentShaderFile;
	const char* GlyphVertexShaderFile;
	const char* GlyphFragmentShaderFile;
	*/



	char VertexShaderFile[MAX_STRING_LENGTH];
	char FragmentShaderFile[MAX_STRING_LENGTH];

	char ParticleVertexShaderFile[MAX_STRING_LENGTH];
	char ParticleFragmentShaderFile[MAX_STRING_LENGTH];
	char ParticleGeometryShaderFile[MAX_STRING_LENGTH];
	char ParticleBackgroundVertexShaderFile[MAX_STRING_LENGTH];
	char ParticleBackgroundFragmentShaderFile[MAX_STRING_LENGTH];
	char GlyphVertexShaderFile[MAX_STRING_LENGTH];
	char GlyphFragmentShaderFile[MAX_STRING_LENGTH];

	VisuType type;
	ParticleVisuType typeParticles;

	GLfloat* particles;
	int nParticles;
	bool showParticles;
	GLfloat* particleMesh;
	int particleMeshRes;
	compute particleMeshSize;
	int nGlyphs;
	int glyphSamplingRateX;
	int glyphSamplingRateY;

	GLfloat glyphScale;

	GLfloat* glyphs;
	GLfloat* glyphMesh;

	// Input variables
	bool mouse1Pressed;
	bool mouse2Pressed;
	double mouse1BeginDrag[2];
	double mouse1EndDrag[2];
	double mouse2BeginDrag[2];
	double mouse2EndDrag[2];
	double mouse2BeginDragShifted[2];

	bool paused;
	bool nonLinItisOver;
	bool initPassivePart;

	GLFWcursor* handCursor;

	unsigned char* imageBuffer; // stores the pixel data from the window to be stored in an image file
	bool writeImages;

	int retinaScale;

	char outputFolder[MAX_STRING_LENGTH];

	bool transparency;
	bool alphaOnValue;


	bool showGlyphs;
	GlyphType glyphType;
	GlyphMeshType glyphMeshType;
	int nGlyphMeshVert;

	bool update;

	int width, height;

	bool updateGrid;

	FilterType filter;


	ColorMap colorMap[MAX_VISU_TYPE];




};
#endif



// Boundary conditions
// ========================
typedef enum {Dirichlet, DirichletGhost, NeumannGhost, Infinity} BCType;
typedef enum {Stokes_PureShear, Stokes_SimpleShear, Stokes_FixedLeftWall, Stokes_Sandbox, Stokes_SandboxWeakBackstop, Stokes_CornerFlow,
			  Thermal_TT_TB_LRNoFlux, Thermal_TT_TBExternal_LRNoFlux,
			  Darcy_Default} SetupType;
typedef struct BC BC;
struct BC
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
};


// Initial conditions
// ========================
typedef enum {IC_HSC, IC_Gaussian} ICSetupType;
typedef struct IC IC;
struct IC {
	ICSetupType SetupType;
	compute data[32];
};

// Equation System
// ========================
typedef struct EqSystem EqSystem;
struct EqSystem
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
};




// Numbering
// ========================
typedef enum {Stencil_Stokes_Momentum_x, Stencil_Stokes_Momentum_y, Stencil_Stokes_Continuity, Stencil_Heat, Stencil_Stokes_Darcy_Momentum_x, Stencil_Stokes_Darcy_Momentum_y, Stencil_Stokes_Darcy_Continuity, Stencil_Stokes_Darcy_Darcy, Stencil_Poisson} StencilType;
typedef struct Numbering Numbering;
struct Numbering
{
	int* map;
	int *IX, *IY;

	int subEqSystem0[NB_SUBSYSTEM_MAX]; // Index of the first equation of a subset (e.g. first Vx, first Vy, first P equation)// hard coded number. Assuming there will never be more than 10 subsystem of equations to solve
	int subEqSystem0Dir[NB_SUBSYSTEM_MAX];
	int nSubEqSystem;
	StencilType Stencil[NB_SUBSYSTEM_MAX];
};
// Inline functions for Numbering

static inline compute shearValue(compute* A, int ix, int iy, int nxEC)
{
	// Compute a value on the shear grid from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy+1)*nxEC] + A[ix+1+(iy+1)*nxEC] + A[ix  +(iy  )*nxEC] + A[ix+1+(iy  )*nxEC])/4;
}

static inline compute centerValue(compute* A, int ix, int iy, int nxS)
{
	// Compute a value on an embedded cell center from the A Array of values defined on the shear grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy-1)*nxS] + A[ix-1+(iy-1)*nxS] + A[ix  +(iy  )*nxS] + A[ix-1+(iy  )*nxS])/4;
}




// Local Numbering Vx
// ========================
typedef struct LocalNumberingVx LocalNumberingVx;
struct LocalNumberingVx
{
	int VxC, VxN, VxS, VxE, VxW;
	int VySW, VyNW, VySE, VyNE;
	int PW, PE;
	int NormalE, NormalW, ShearN, ShearS;
};

// Local Numbering Vy
// ========================
typedef struct LocalNumberingVy LocalNumberingVy;
struct LocalNumberingVy
{
	int VxSW, VxSE, VxNW, VxNE;
	int VyS, VyW, VyC, VyE, VyN;
	int PS, PN;
	int NormalN, NormalS, ShearE, ShearW;
};

// Local Numbering P
// ========================
typedef struct LocalNumberingP LocalNumberingP;
struct LocalNumberingP
{
	int VxE, VxW, VyN, VyS;
};



// A basic linked list node struct
// ========================
typedef struct LinkedNode LinkedNode;
struct LinkedNode {
    int data;
    struct LinkedNode* next;
};


// Pardiso solver
// ========================
typedef struct Solver Solver;
struct Solver {
	// Pardiso struct
	void    *pt[64];
	int      iparm[64];
	double   dparm[64];
	int      maxfct, mnum, msglvl, mtype, nrhs;
};




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
void Grid_allocateMemory	(Grid* Grid);
void Grid_freeMemory		(Grid* Grid);
void Grid_init				(Grid* Grid, Numerics* Numerics);
void Grid_updatePureShear	(Grid* Grid, BC* BC, Numerics* Numerics, compute dt);




// Particles
// =========================
void Particles_allocateMemory 			(Particles* Particles, Grid* Grid);
void Particles_freeMemory	 			(Particles* Particles, Grid* Grid);
void Particles_initCoord				(Particles* Particles, Grid* Grid);
void Particles_initPassive				(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_updateLinkedList 		(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_injectOrDelete			(Particles* Particles, Grid* Grid);
void Particles_injectAtTheBoundaries	(Particles* Particles, Grid* Grid, Physics* Physics, MatProps* MatProps);
void Particles_advect					(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_Periodicize				(Particles* Particles, Grid* Grid);
void Particles_teleportInsideTheDomain	(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_deleteIfOutsideTheDomain	(Particles* Particles, Grid* Grid);
void Particles_switchStickyAir			(Particles* Particles, Grid* Grid, Physics* Physics, Numerics* Numerics, MatProps* MatProps);
void addToParticlePointerList 			(ParticlePointerList** pointerToHead, SingleParticle* thisParticle);
void freeParticlePointerList			(ParticlePointerList* head);
void Particles_freeAllSingleParticles	(Particles* Particles, Grid* Grid);
void addSingleParticle					(SingleParticle** pointerToHead, SingleParticle* modelParticle);




// Physics
// =========================
void Physics_allocateMemory						(Physics* Physics, Grid* Grid);
void Physics_freeMemory							(Physics* Physics, Grid* Grid);
void Physics_initPToLithostatic					(Physics* Physics, Grid* Grid);
void Physics_interpFromParticlesToCell			(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal, BC* BCThermal);
void Physics_interpFromCellToNode				(Grid* Grid, compute* CellValue, compute* NodeValue);
void Physics_interpTempFromCellsToParticle		(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes, MatProps* MatProps, BC* BCThermal);
void Physics_interpStressesFromCellsToParticle	(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal, MatProps* MatProps);
void Physics_get_VxVy_FromSolution				(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem, Numerics* Numerics);
#if (CRANK_NICHOLSON_VEL)
void Physics_updateOldVel						(Physics* Physics, Grid* Grid);
#endif
void Physics_get_P_FromSolution					(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem, Numerics* Numerics);
void Physics_get_T_FromSolution					(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem, Numerics* Numerics);
void Physics_computeStrainRateInvariant			(Physics* Physics, Grid* Grid, compute* StrainRateInvariant);
void Physics_initEta							(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics);
void Physics_computeEta							(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes, MatProps* MatProps);
void Physics_computeStressChanges				(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes);
void Physics_interpPhiFromCellsToParticle		(Grid* Grid, Particles* Particles, Physics* Physics);
void Physics_interpStrainFromCellsToParticle	(Grid* Grid, Particles* Particles, Physics* Physics);
void Physics_changePhaseOfFaults				(Physics* Physics, Grid* Grid, MatProps* MatProps, Particles* Particles);
void Physics_updateDt							(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics);
void Physics_computeStrainRateInvariantForOneCell	(Physics* Physics, Grid* Grid, int ix, int iy, compute* EII);
void Physics_computeStrainRateInvariantForOneNode	(Physics* Physics, BC* BCStokes, Grid* Grid, int ix, int iy, compute* EII);
void Physics_computeStressInvariantForOneCell	(Physics* Physics, Grid* Grid, int ix, int iy, compute* SII);
#if (DARCY)
void Physics_computePerm						(Physics* Physics, Grid* Grid, Numerics* Numerics, MatProps* MatProps);
void Physics_computePhi							(Physics* Physics, Grid* Grid, Numerics* Numerics);
#endif
void Physics_copyValuesToSides					(compute* ECValues, Grid* Grid);
void Physics_copyValuesToSidesi					(int* ECValues, Grid* Grid);
void Physics_computeRho							(Physics* Physics, Grid* Grid, MatProps* MatProps);
void Physics_get_ECVal_FromSolution 			(compute* Val, int ISub, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem);
void Physics_getPhase 							(Physics* Physics, Grid* Grid, Particles* Particles, MatProps* MatProps, BC* BCStokes);


//void Physics_computePlitho						(Physics* Physics, Grid* Grid);

void Physics_getValuesToSidesFromBC(compute* ECValues, Grid* Grid, BC* BC, Numbering* Numbering);
compute Physics_computeSideValuesFromBC_ForOneCell(compute neighValue, BC* BC, int IBC, int ix, int iy, Grid* Grid);
void Physics_reinitPhaseList(Physics* Physics, Grid* Grid) ;
void Physics_check(Physics* Physics, Grid* Grid, Char* Char) ;
compute Physics_getFromMatProps_ForOneCell(Physics* Physics, compute* ListFromMatProps, MatProps* MatProps, int iCell);

// Visualization
// =========================
#if (VISU)
	void Visu_allocateMemory	(Visu* Visu, Grid* Grid );
	void Visu_freeMemory		(Visu* Visu );
	void Visu_init				(Visu* Visu, Grid* Grid, Particles* Particles, Char* Char, Input* Input);
	void Visu_updateVertices	(Visu* Visu, Grid* Grid);
	void Visu_initWindow		(Visu* Visu);
	void error_callback			(int error, const char* description);
	void key_callback			(GLFWwindow* window, int key, int scancode, int action, int mods);

	void Visu_updateCenterValue (Visu* Visu, Grid* Grid, compute* CellValue);
	void Visu_updateCenterValuei(Visu* Visu, Grid* Grid, int* CellValue);
	void Visu_StrainRate		(Visu* Visu, Grid* Grid, Physics* Physics);
	void Visu_updateUniforms	(Visu* Visu);
	void Visu_velocity			(Visu* Visu, Grid* Grid, Physics* Physics);
	void VisudivV				(Visu* Visu, Grid* Grid, Physics* Physics);
	void Visu_stress			(Visu* Visu, Grid* Grid, Physics* Physics);
	void Visu_SIIOvYield		(Visu* Visu, Grid* Grid, Physics* Physics, Numerics* Numerics);
	void Visu_PeOvYield			(Visu* Visu, Grid* Grid, Physics* Physics, Numerics* Numerics);
	void Visu_update			(Visu* Visu, Grid* Grid, Physics* Physics, Char* Char, EqSystem* EqStokes, EqSystem* EqThermal, Numbering* NumStokes, Numbering* NumThermal, Numerics* Numerics);
	void Visu_checkInput		(Visu* Visu);
	void Visu_particles			(Visu* Visu, Particles* Particles, Grid* Grid);
	void Visu_glyphs			(Visu* Visu, Physics* Physics, Grid* Grid, Particles* Particles);
	void Visu_particleMesh		(Visu* Visu);
	void Visu_alphaValue		(Visu* Visu, Grid* Grid, Physics* Physics);
	void Visu_glyphMesh			(Visu* Visu);

	void Visu_main				(Visu* Visu, Grid* Grid, Physics* Physics, Particles* Particles, Numerics* Numerics, Char* Char, EqSystem* EqStokes, EqSystem* EqThermal, Numbering* NumStokes, Numbering* NumThermal);
	void Visu_residual			(Visu* Visu, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);
#endif



// Boundary conditions
// =========================
void BC_freeMemory			(BC* BC);
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
void Numbering_allocateMemory	(Numbering* Numbering, EqSystem* EqSystem, Grid* Grid);
void Numbering_freeMemory		(Numbering* Numbering);
//inline compute shearValue(compute* A, int ix, int iy, int nxEC) __attribute__((always_inline));
void Numbering_init				(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering, Physics* Physics);



// Equation system
// =========================
void EqSystem_allocateI		(EqSystem* EqSystem);
void EqSystem_allocateMemory(EqSystem* EqSystem);
void EqSystem_freeMemory	(EqSystem* EqSystem, Solver* Solver) ;
void EqSystem_assemble		(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering, bool updateScale);

void EqSystem_solve			(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_check			(EqSystem* EqSystem);
void EqSystem_initSolver  	(EqSystem* EqSystem, Solver* Solver);
void pardisoSolveSymmetric	(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_computePressureAndUpdateRHS(EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC);
void EqSystem_computeNormResidual(EqSystem* EqSystem);
void EqSystem_scale			(EqSystem* EqSystem);
void EqSystem_unscale		(EqSystem* EqSystem);

// Local stencil
// =========================
void LocalStencil_Call(StencilType Stencil, int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
void LocalStencil_Stokes_Momentum_x	(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
void LocalStencil_Stokes_Momentum_y	(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* Ic);
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
void Numerics_freeMemory(Numerics* Numerics);
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




