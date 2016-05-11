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

//#ifdef __APPLE__
	#include <GL/glew.h>
//#endif
#include <GLFW/glfw3.h>

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#include <time.h>
#include <png.h>





//============================================================================//
//============================================================================//
//                                                                            //
//                                   MACROS                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
#define DEBUG 	false
#define VISU 	true
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
#define TIC tic = glfwGetTime();
#define TOC toc = glfwGetTime(); \
            toc = toc - tic;


#define PI 		acos(-1.0)

#define WIDTH 1080
#define HEIGHT 1080

#define INIT_PARTICLE SingleParticle* thisParticle = NULL; \
						int iNode = 0;

#define FOR_PARTICLES  	  for (iNode = 0; iNode < Grid->nSTot; ++iNode) { \
								thisParticle = Particles->linkHead[iNode]; \
								while (thisParticle != NULL) {

#define END_PARTICLES  			thisParticle = thisParticle->next; \
								} \
							}

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

	compute g[2]; // gravity acceleration
	compute dt;
	compute *Vx, *Vy, *P;
	compute maxV;
	compute *eta; // Viscosity
	compute *etaShear;
	compute *eta0, *n;
	compute epsRef; // reference strainrate

	compute *rho; // Density
	compute *k;  // Thermal conductivity
	compute Cp; // heat capacity, taken as a single value because it varies very little between different types of rocks
	compute etaMin, etaMax;

	compute *T, *DT; // temperature stored on cell centers

	compute *G, *GShear; // shear modulus

	compute *sigma_xx_0, *sigma_xy_0; // old stresses
	compute *Dsigma_xx_0, *Dsigma_xy_0; // stress corrections for markers

	compute *cohesion, *frictionAngle;

	int itNonLin;
	compute glob; // globalization factor

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
	compute dx, dy; 		 	  	// grid increment / cell size

};



// Material properties
// =========================
typedef enum {LinearViscous, PowerLawViscous, ViscoElastic} FlowLaw;
typedef struct MatProps MatProps;
struct MatProps
{
	int nPhase;
	compute rho0[NB_PHASE_MAX], eta0[NB_PHASE_MAX], n[NB_PHASE_MAX];
	compute alpha[NB_PHASE_MAX]; // thermal expansion
	compute beta[NB_PHASE_MAX];  // compressibility
	compute k[NB_PHASE_MAX]; 	 // thermal conductivity
	compute G[NB_PHASE_MAX]; 	 // shear modulus
	FlowLaw flowLaw[NB_PHASE_MAX];
	compute maxwellTime[NB_PHASE_MAX]; // Mtime = eta/G
	compute cohesion[NB_PHASE_MAX]; // cohesion
	compute frictionAngle[NB_PHASE_MAX]; // angle of friction
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


	// for the linked list
	int nodeId;
    struct SingleParticle* next;
};

// Id Changed
typedef struct ParticlePointerList ParticlePointerList;
struct ParticlePointerList {
    //int data;
    SingleParticle* pointer;
    struct ParticlePointerList* next;
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
typedef enum {Blank, Viscosity, StrainRate, Velocity, Pressure, Density, Temperature, Stress} VisuType;
typedef enum {Phase, PartTemp,PartSigma_xx, PartSigma_xy} ParticleVisuType;
typedef struct Visu Visu;
struct Visu
{
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
	GLuint ShaderProgram, ParticleShaderProgram, ParticleBackgroundShaderProgram;
	const char* VertexShaderFile;
	const char* FragmentShaderFile;

	const char* ParticleVertexShaderFile;
	const char* ParticleFragmentShaderFile;
	const char* ParticleGeometryShaderFile;
	const char* ParticleBackgroundVertexShaderFile;
	const char* ParticleBackgroundFragmentShaderFile;

	VisuType type;
	ParticleVisuType typeParticles;

	GLfloat* particles;
	int nParticles;
	bool showParticles;
	GLfloat* particleMesh;
	int particleMeshRes;

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

};




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


// Memory
// =========================
void Memory_allocateMain	(Grid* Grid, Particles* Particles, Physics* Physics, EqSystem* EqStokes, Numbering* NumStokes, Numbering* NumThermal);
void Memory_freeMain		(Particles* Particles, Physics* Physics, Numbering* NumStokes, Numbering* NumThermal, BC* BC, Grid* Grid);
void addToLinkedList		(LinkedNode** pointerToHead, int x);
void freeLinkedList			(LinkedNode* head);




// Char
// =========================
void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BCStokes, BC* BCThermal);




// Grid
// =========================
void Grid_updatePureShear(Grid* Grid, BC* BC, compute dt);




// Particles
// =========================
void Particles_initCoord		(Grid* Grid, Particles* Particles);
void Particles_initPhase		(Grid* Grid, Particles* Particles);
void Particles_initPassive		(Grid* Grid, Particles* Particles);
void Particles_initPhysics		(Grid* Grid, Particles* Particles, BC* BC);
void Particles_updateLinkedList(Grid* Grid, Particles* Particles, Physics* Physics);
void Particles_advect			(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_Periodicize		(Grid* Grid, Particles* Particles, BC* BC);
void Particles_teleportInsideTheDomain(Grid* Grid, Particles* Particles, Physics* Physics);
void Particles_deleteIfOutsideTheDomain(Grid* Grid, Particles* Particles);
void addToParticlePointerList 	(ParticlePointerList** pointerToHead, SingleParticle* thisParticle);
void freeParticlePointerList	(ParticlePointerList* head);
void Particles_freeAllSingleParticles	(Particles* Particles, Grid* Grid);
void addSingleParticle(SingleParticle** pointerToHead, SingleParticle* modelParticle);




// Physics
// =========================
void Physics_interpFromParticlesToCell	(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps, BC* BCStokes, Numbering* NumThermal, BC* BCThermal);
void Physics_interpFromCellToNode		(Grid* Grid, compute* CellValue, compute* NodeValue);
void Physics_interpTempFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal);
void Physics_interpStressesFromCellsToParticle(Grid* Grid, Particles* Particles, Physics* Physics, BC* BCStokes,  BC* BCThermal, Numbering* NumThermal);
void Physics_set_VxVyP_FromSolution		(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem);
void Physics_set_T_FromSolution			(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, EqSystem* EqSystem);
void Physics_computeStrainRateInvariant	(Physics* Physics, Grid* Grid, compute* StrainRateInvariant);
void Physics_computeEta					(Physics* Physics, Grid* Grid);
void Physics_computeStressChanges		(Physics* Physics, Grid* Grid, BC* BC, Numbering* NumStokes, EqSystem* EqStokes);



// Visualization
// =========================
void Visu_allocateMemory	(Visu* Visu, Grid* Grid );
void Visu_freeMemory		(Visu* Visu );
void Visu_init				(Visu* Visu, Grid* Grid, Particles* Particles);
void Visu_updateVertices	(Visu* Visu, Grid* Grid);
void Visu_initWindow		(GLFWwindow** window, Visu* Visu);
void Visu_initOpenGL		(Visu* Visu, Grid* Grid, GLFWwindow* window);
void error_callback			(int error, const char* description);
void key_callback			(GLFWwindow* window, int key, int scancode, int action, int mods);

void Visu_updateCenterValue (Visu* Visu, Grid* Grid, compute* CellValue, int BCType);
void Visu_StrainRate		(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC);
void Visu_updateUniforms	(Visu* Visu, GLFWwindow* window);
void Visu_velocity			(Visu* Visu, Grid* Grid, Physics* Physics);
void Visu_stress			(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC);
void Visu_update			(Visu* Visu, GLFWwindow* window, Grid* Grid, Physics* Physics, BC* BC, Char* Char);
void Visu_checkInput		(Visu* Visu, GLFWwindow* window);
void Visu_particles			(Visu* Visu, Particles* Particles, Grid* Grid);
void Visu_particleMesh		(Visu* Visu);
void Visu_alphaValue		(Visu* Visu, Grid* Grid, Particles* Particles);




// Boundary conditions
// =========================
void BC_initStokes		(BC* BC, Grid* Grid, EqSystem* EqSystem);
void BC_initThermal		(BC* BC, Grid* Grid, EqSystem* EqSystem);
void BC_updateStokes	(BC* BC, Grid* Grid);
void BC_updateThermal	(BC* BC, Grid* Grid);





// Numbering
// =========================
void Numbering_init			(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);
void Numbering_getLocalNNZ	(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, bool useNumMap, StencilType StencilType, int* sum);
void Numbering_initThermal	(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);



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
void compileShaders(GLuint *ShaderProgram, const char* pVSFileName, const char* pFSFileName, const char* pGSFileName, bool useGS);
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
int writePNGImage(char* filename, int width, int height, unsigned char *buffer, char* title);

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


#endif /* STOKES_H_ */




