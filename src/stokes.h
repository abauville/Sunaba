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

#define DEBUG 	false
#define VISU 	true
#define HEAT	true



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
#define NXC 10
#define NYC 10
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


#define PI 		acos(-1.0)


#define INIT_PARTICLE SingleParticle* thisParticle = NULL; \
						int iNode = 0;

#define FOR_PARTICLES  	  for (iNode = 0; iNode < Grid->nSTot; ++iNode) { \
								thisParticle = Particles->linkHead[iNode]; \
								while (thisParticle != NULL) {

#define END_PARTICLES  			thisParticle = thisParticle->next; \
								} \
							}

#define FAULT_MOD 10.0/1.0

#define MAX_STRING_LENGTH 512

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
};


// Numerics
// =========================
typedef struct Numerics Numerics;
struct Numerics
{

	int timeStep;

	int nTimeSteps; //  negative value for infinite
	int nLineSearch;
	int maxNonLinearIter; // should always be greater than the number of line searches
	int minNonLinearIter; // should always be greater than the number of line searches
	compute relativeTolerance; // relative tolerance to the one of this time step
	compute absoluteTolerance; // relative tolerance to the first one of the simulation
	compute maxCorrection;

	int itNonLin;
	compute *glob; // globalization factor

	compute CFL_fac;
	compute dtmin, dtmax;
	compute etaMin, etaMax;
	compute dLmin; // min grid size

	compute minRes;

	compute normRes0, normResRef;
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
typedef struct Physics Physics;
struct Physics
{

	// Physics Stokes
	compute g[2]; // gravity acceleration
	compute dt;
	compute *Vx, *Vy, *P;
	compute maxV;
	compute *eta;

	compute *eta0, *n;
	compute *rho; // Density

#if (HEAT)
	compute *k;  // Thermal conductivity
	compute *T, *DT; // temperature stored on cell centers
#endif

	compute epsRef; // reference strainrate

	compute *psi, *Dpsi; // pressure head for Darcy stored on cell centers
	compute *kD; // Darcy diffusion term
	compute *SD; // Darcy other term


	// Stokes, elasticity related variables
	compute *sigma_xx_0, *sigma_xy_0; // old stresses
	compute *Dsigma_xx_0, *Dsigma_xy_0; // stress corrections for markers
	compute *G; // shear modulus

	// Plasticity
	compute *cohesion, *frictionAngle;

	compute *etaVisc;


	// Physics thermal

	compute Cp; // heat capacity, taken as a single value because it varies very little between different types of rocks



	// Darcy




	compute time;

	// compute stressOld
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


};



// Material properties
// =========================
typedef struct MatProps MatProps;
typedef enum {Custom} MaterialType;
struct MatProps
{
	int nPhase;

	MaterialType Material;

	char name[NB_PHASE_MAX][128];

	compute rho0[NB_PHASE_MAX], eta0[NB_PHASE_MAX], n[NB_PHASE_MAX];
	compute alpha[NB_PHASE_MAX]; // thermal expansion
	compute beta[NB_PHASE_MAX];  // compressibility
	compute k[NB_PHASE_MAX]; 	 // thermal conductivity
	compute G[NB_PHASE_MAX]; 	 // shear modulus

	compute cohesion[NB_PHASE_MAX]; // cohesion
	compute frictionAngle[NB_PHASE_MAX]; // angle of friction

	compute kD[NB_PHASE_MAX]; // Darcy diffusivity in km/yr
	compute SD[NB_PHASE_MAX]; // Darcy stuff in 1/km

};





// Particles
// =========================

// Single Particle storing coordinate, temp and info for a linked list
typedef struct SingleParticle SingleParticle;
struct SingleParticle {
	coord x, y;
	int phase;
	float passive; // some passive attribute used for visualization
	compute T;

	// Old stresses
	compute sigma_xx_0;
	compute sigma_xy_0;

	compute psi;

	bool faulted;

	// for the linked list
	int nodeId;
    struct SingleParticle* next;

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
struct Particles
{
	int nPC, nPCX, nPCY; // number of particles per cell, tot, in x and in y
	int n; // number of particles
	coord noiseFactor;
	SingleParticle **linkHead;


};





// Visualization
// ========================
#if (VISU)
typedef enum {Blank, Viscosity, StrainRate, Velocity, Pressure, Density, Temperature, Stress, WaterPressureHead, Permeability} VisuType;
typedef enum {Phase, PartTemp,PartSigma_xx, PartSigma_xy} ParticleVisuType;
typedef enum {StokesVelocity, DarcyGradient} GlyphType;
typedef enum {Triangle, ThinArrow, ThickArrow} GlyphMeshType;
typedef enum {Nearest, Linear} FilterType;
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
	const char* VertexShaderFile;
	const char* FragmentShaderFile;

	const char* ParticleVertexShaderFile;
	const char* ParticleFragmentShaderFile;
	const char* ParticleGeometryShaderFile;
	const char* ParticleBackgroundVertexShaderFile;
	const char* ParticleBackgroundFragmentShaderFile;
	const char* GlyphVertexShaderFile;
	const char* GlyphFragmentShaderFile;

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

	bool initPassivePart;

	GLFWcursor* handCursor;

	unsigned char* imageBuffer; // stores the pixel data from the window to be stored in an image file
	bool writeImages;

	int retinaScale;

	char outputFolder[1024];

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



};
#endif



// Boundary conditions
// ========================
typedef enum {Dirichlet, DirichletGhost, NeumannGhost} BCType;
typedef enum {PureShear, SimpleShearPeriodic, FixedLeftWall, Sandbox} SetupType;
typedef struct BC BC;
struct BC
{
	int n;


	int *list;
	BCType *type;
	compute *value;




	bool hPeriod;
	bool mobileBox;

	//int typeL, typeR, typeT, typeB;
	SetupType SetupType;
	// Should be moved somewhere else, but I don't know where yet
	compute backStrainRate; // background strain, positive in extension
	compute TB, TT;
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

	compute normResidual;
};




// Numbering
// ========================
typedef enum {Vx, Vy, P, T} StencilType;
typedef struct Numbering Numbering;
struct Numbering
{
	int* map;
	int *IX, *IY;

	int subEqSystem0[10]; // Index of the first equation of a subset (e.g. first Vx, first Vy, first P equation)// hard coded number. Assuming there will never be more than 10 subsystem of equations to solve
	int nSubEqSystem;
	StencilType Stencil[10];
};
// Inline functions for Numbering

static inline compute shearValue(compute* A, int ix, int iy, int nxEC)
{
	// Compute a value on the shear grid from a Array of values defined on the Embedded cell grid
	// where ix and iy refer to shear node grid
	return(A[ix  +(iy+1)*nxEC] + A[ix+1+(iy+1)*nxEC] + A[ix  +(iy  )*nxEC] + A[ix+1+(iy  )*nxEC])/4;
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



// Darcy
// ========================
typedef struct Darcy Darcy;
struct Darcy {
	coord 	hOcean;
	compute rainFlux;
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
void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BCStokes, BC* BCThermal);




// Grid
// =========================
void Grid_updatePureShear(Grid* Grid, BC* BC, compute dt);




// Particles
// =========================
void Particles_allocateMemory 	(Particles* Particles, Grid* Grid);
void Particles_freeMemory	 	(Particles* Particles, Grid* Grid);
void Particles_initCoord		(Particles* Particles, Grid* Grid);
void Particles_initPassive		(Particles* Particles, Grid* Grid);
void Particles_initPhysics		(Particles* Particles, Grid* Grid, BC* BC);
void Particles_updateLinkedList (Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_advect			(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_Periodicize		(Particles* Particles, Grid* Grid, BC* BC);
void Particles_teleportInsideTheDomain	(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_deleteIfOutsideTheDomain	(Particles* Particles, Grid* Grid);
void addToParticlePointerList 			(ParticlePointerList** pointerToHead, SingleParticle* thisParticle);
void freeParticlePointerList			(ParticlePointerList* head);
void Particles_freeAllSingleParticles	(Particles* Particles, Grid* Grid);
void addSingleParticle			(SingleParticle** pointerToHead, SingleParticle* modelParticle);




// Physics
// =========================
void Physics_allocateMemory						(Physics* Physics, Grid* Grid);
void Physics_freeMemory							(Physics* Physics);
void Physics_initPToLithostatic					(Physics* Physics, Grid* Grid);
void Physics_interpFromParticlesToCell			(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal, BC* BCThermal);
void Physics_interpFromCellToNode				(Grid* Grid, compute* CellValue, compute* NodeValue);
void Physics_interpTempFromCellsToParticle		(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal);
void Physics_interpStressesFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal, MatProps* MatProps);
void Physics_get_VxVyP_FromSolution				(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem);
void Physics_get_T_FromSolution					(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem);
void Physics_computeStrainRateInvariant			(Physics* Physics, Grid* Grid, compute* StrainRateInvariant);
void Physics_computeEta							(Physics* Physics, Grid* Grid, Numerics* Numerics);
void Physics_computeStressChanges				(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes);
void Physics_interpPsiFromCellsToParticle		(Grid* Grid, Particles* Particles, Physics* Physics);
void Physics_changePhaseOfFaults				(Physics* Physics, Grid* Grid, MatProps* MatProps, Particles* Particles);
void Physics_updateDt							(Physics* Physics, Numerics* Numerics);
void Physics_computeStrainInvariantForOneCell		(Physics* Physics, Grid* Grid, int ix, int iy, compute* EII);


// Visualization
// =========================
#if (VISU)
	void Visu_allocateMemory	(Visu* Visu, Grid* Grid );
	void Visu_freeMemory		(Visu* Visu );
	void Visu_init				(Visu* Visu, Grid* Grid, Particles* Particles);
	void Visu_updateVertices	(Visu* Visu, Grid* Grid);
	void Visu_initWindow		(Visu* Visu);
	void error_callback			(int error, const char* description);
	void key_callback			(GLFWwindow* window, int key, int scancode, int action, int mods);

	void Visu_updateCenterValue (Visu* Visu, Grid* Grid, compute* CellValue, int BCType);
	void Visu_StrainRate		(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC);
	void Visu_updateUniforms	(Visu* Visu);
	void Visu_velocity			(Visu* Visu, Grid* Grid, Physics* Physics);
	void Visu_stress			(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC);
	void Visu_update			(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC, Char* Char, MatProps* MatProps);
	void Visu_checkInput		(Visu* Visu);
	void Visu_particles			(Visu* Visu, Particles* Particles, Grid* Grid);
	void Visu_glyphs			(Visu* Visu, Physics* Physics, Grid* Grid, Particles* Particles);
	void Visu_particleMesh		(Visu* Visu);
	void Visu_alphaValue		(Visu* Visu, Grid* Grid, Particles* Particles);
	void Visu_glyphMesh			(Visu* Visu);

	void Visu_main				(Visu* Visu, Grid* Grid, Physics* Physics, Particles* Particles, Numerics* Numerics, BC* BCStokes, Char* Char, MatProps* MatProps);
#endif



// Boundary conditions
// =========================
void BC_freeMemory		(BC* BC);
void BC_initStokes		(BC* BC, Grid* Grid, EqSystem* EqSystem);
void BC_initThermal		(BC* BC, Grid* Grid, EqSystem* EqSystem);
void BC_updateStokes	(BC* BC, Grid* Grid);
void BC_updateThermal	(BC* BC, Grid* Grid);





// Numbering
// =========================
void Numbering_allocateMemory	(Numbering* Numbering, EqSystem* EqSystem, Grid* Grid);
void Numbering_freeMemory		(Numbering* Numbering);
//inline compute shearValue(compute* A, int ix, int iy, int nxEC) __attribute__((always_inline));
void Numbering_init				(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);
void Numbering_getLocalNNZ		(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, bool useNumMap, StencilType StencilType, int* sum);



// Equation system
// =========================
void EqSystem_allocateI		(EqSystem* EqSystem);
void EqSystem_allocateMemory(EqSystem* EqSystem);
void EqSystem_freeMemory	(EqSystem* EqSystem, Solver* Solver) ;
void EqSystem_assemble		(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering);

void EqSystem_solve			(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_check			(EqSystem* EqSystem);
void EqSystem_initSolver  	(EqSystem* EqSystem, Solver* Solver);
void pardisoSolveSymmetric	(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering);
void EqSystem_computePressureAndUpdateRHS(EqSystem* EqSystem, Grid* Grid, Numbering* Numbering, Physics* Physics, BC* BC);
void EqSystem_computeNormResidual(EqSystem* EqSystem);



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
typedef enum {Air, Ocean, Solid} PhaseFlag;
void Darcy_setBC		(Grid* Grid, Physics* Physics, coord hOcean, PhaseFlag* Phase);
void Darcy_setPhaseFlag	(PhaseFlag* Phase, coord hOcean, Grid* Grid, Particles* Particles);
void Darcy_solve		(Darcy* Darcy, Grid* Grid, Physics* Physics, MatProps* MatProps, Particles* Particles);


// Numerics
// ========================
void Numerics_init		(Numerics* Numerics);
void Numerics_freeMemory(Numerics* Numerics);
int  Numerics_updateBestGlob(Numerics* Numerics, EqSystem* EqStokes, int iLS);

// Input
// ========================
void Input_read(Input* Input, Grid* Grid, Numerics* Numerics, Physics* Physics, MatProps* MatProps, Particles* Particles, Char* Char, BC* BCStokes, BC* BCThermal);
void Input_assignPhaseToParticles(Input* Input, Particles* Particles, Grid* Grid, Char* Char);

#if (VISU)
void Input_readVisu(Input* Input, Visu* Visu);
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




