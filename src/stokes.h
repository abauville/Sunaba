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

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <stdbool.h>
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
#define DEBUG 	false
#define VISU 	true
#define ADVECT  false
#define NB_PHASE_MAX 10
#define NXC 10
#define NYC 10
#define UPPER_TRI true
#define TIMER false

#define INIT_TIMER 	clock_t tic, diff; \
					float toc;
#define TIC 	tic = clock();
#define TOC 	diff = (clock() - tic); \
				toc = (float) (diff * 1000 /CLOCKS_PER_SEC)/1000;

#define PI 		acos(-1.0);

#define WIDTH 1024
#define HEIGHT 1024




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
	compute velocity; 		// [m.s-1]
	compute density; 		// [kg.m^-3]
	compute stress;  		// [Pa] or [kg.m^-1.s^-2]
	compute viscosity; 		// [Pa.s] or [kg.m^-1.s-1]
	compute acceleration; 	// [m.s^-2]
	compute strainrate; 	// [s^-1]
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
	compute *rho; // Density
	// compute stressOld
};




// Grid
// =========================
typedef struct Grid Grid;
struct Grid
{
	int nxC, nyC, nCTot; 			// number of cells / cell center nodes
	int nxS, nyS, nSTot; 			// number of base nodes
	int nxVx, nyVx, nVxTot; 		// number of Vx nodes
	int nxVy, nyVy, nVyTot; 		// number of Vy nodes
	coord xmin, xmax, ymin, ymax; 	// grid extent
	compute dx, dy; 		 	  	// grid increment / cell size

};



// Material properties
// =========================
typedef struct MatProps MatProps;
struct MatProps
{
	int nPhase;
	compute rho0[NB_PHASE_MAX], eta0[NB_PHASE_MAX];
};




// Particles
// =========================
typedef struct Particles Particles;
struct Particles
{
	int nPC, nPCX, nPCY; // number of particles per cell, tot, in x and in y
	int n; // number of particles
	coord *xy;
	int *phase; // i is the index of the cell in which the particle is

	// Doubly linked list
	int *cellId;
	int *linkNext, *linkHead;
};




// Visualization
// ========================
typedef enum {Viscosity, StrainRate} VisuType;
typedef struct Visu Visu;
struct Visu
{
	int ntri, ntrivert;
	GLuint* elements;
	GLfloat* U;
	GLfloat* vertices;
	GLfloat scale, valueScale;
	GLfloat colorScale[2];
	GLfloat shift[2];
	GLint log10_on;
	GLuint VAO, VBO, CBO, EBO;
	GLuint ShaderProgram;
	const char* VertexShaderFile;
	const char* FragmentShaderFile;
	VisuType type;

	// Input variables
	bool mouse1Pressed;
	bool mouse2Pressed;
	double mouse1BeginDrag[2];
	double mouse1EndDrag[2];
	double mouse2BeginDrag[2];
	double mouse2EndDrag[2];

	GLFWcursor* handCursor;
};




// Boundary conditions
// ========================
typedef struct BC BC;
struct BC
{
	int nDir, nPDir, nNeu;
	int *listDir, *listNeu, *listNeuNeigh;
	compute *valueDir, *valueNeu, *coeffNeu, *coeffNeuNeigh;
	bool *isNeu;
	compute VxL, VxR, VxT, VxB;
	compute VyL, VyR, VyT, VyB;
	compute backStrainRate; // background strain, positive in extension
	//int typeL, typeR, typeT, typeB;
	int SetupType;
};




// Equation System
// ========================
typedef struct EqSystem EqSystem;
struct EqSystem
{
	int nEqIni, nEq, nRow, nnz;
	int VyEq0, PEq0;
	// Stiffness matrix
	int *I; // Allocated in function Numbering_initMapAndSparseTripletIJ
	int *J;
	compute *V;

	compute *b; // right hand side
	compute *x; // solution vector;

};




// Numbering
// ========================
typedef struct Numbering Numbering;
struct Numbering
{
	int* map;
	int *IX, *IY;
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




//============================================================================//
//============================================================================//
//                                                                            //
//                                PROTOTYPES                            	  //
//                                                                            //
//============================================================================//
//============================================================================//


// Memory
// =========================
void Memory_allocateMain	(Grid* Grid, Particles* Particles, Physics* Physics, EqSystem* EqSystem, Numbering* Numbering);
void Memory_freeMain		(Particles* Particles, Physics* Physics, Numbering* Numbering);
void addToLinkedList		(LinkedNode** pointerToHead, int x);
void freeLinkedList			(LinkedNode* head);




// Char
// =========================
void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BC);




// Grid
// =========================
void Grid_updatePureShear(Grid* Grid, BC* BC, compute dt);




// Particles
// =========================
void Particles_initCoord		(Grid* Grid, Particles* Particles);
void Particles_initPhase		(Grid* Grid, Particles* Particles);
void Particles_updateLinkedList (Grid* Grid, Particles* Particles);
void Particles_advect			(Particles* Particles, Grid* Grid, Physics* Physics);
void Particles_Periodicize		(Grid* Grid, Particles* Particles, BC* BC);


// Physics
// =========================
void Physics_interpFromParticlesToCell	(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps);
void Physics_interpFromCellToNode	  	(Grid* Grid, compute* CellValue, compute* NodeValue, int BCType);
void Physics_set_VxVyP_FromSolution   	(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, compute* sol);




// Visualization
// =========================
void Visu_allocateMemory	(Visu* Visu, Grid* Grid );
void Visu_freeMemory		(Visu* Visu );
void Visu_init				(Visu* Visu, Grid* Grid);
void Visu_updateVertices	(Visu* Visu, Grid* Grid);
void Visu_initWindow		(GLFWwindow** window);
void Visu_initOpenGL		(Visu* Visu, Grid* Grid);
void error_callback			(int error, const char* description);
void key_callback			(GLFWwindow* window, int key, int scancode, int action, int mods);

void Visu_updateCenterValue (Visu* Visu, Grid* Grid, compute* CellValue, int BCType);
void Visu_StrainRate		(Visu* Visu, Grid* Grid, Physics* Physics, BC* BC);
void Visu_updateUniforms	(Visu* Visu, GLFWwindow* window);

void Visu_update			(Visu* Visu, GLFWwindow* window, Grid* Grid, Physics* Physics, BC* BC, Char* Char);

// Boundary conditions
// =========================
void BC_set(BC* BC, Grid* Grid, EqSystem* EqSystem, Physics* Physics);
void BC_updateDir(BC* BC, Grid* Grid);
void BC_numberNeu(BC* BC, Grid* Grid, EqSystem* EqSystem);
void BC_updateNeuCoeff(BC* BC, Grid* Grid, Physics* Physics);


// Numbering
// =========================
void Numbering_initMapAndSparseTripletIJ(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);
void Numbering_getLocalVx(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingVx* LocVx, bool useNumMap);
void Numbering_getLocalVy(int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingVy* LocVy, bool useNumMap);
void Numbering_getLocalP (int ix, int iy, Numbering* Numbering, Grid* Grid, BC* BC, LocalNumberingP*  LocP , bool useNumMap);




// Equation system
// =========================
void EqSystem_allocateI(EqSystem* EqSystem);
void EqSystem_allocateMemory(EqSystem* EqSystem);
void EqSystem_freeMemory(EqSystem* EqSystem);
void EqSystem_assemble(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering);
void fill_J_V_local(int Type, int BCtype, int ix, int iy,int I, int iEq, int *Jsparse, compute* Vsparse, compute* RHS, int nxC, int nyC, compute dx, compute dy, int *NumMap, compute* EtaShear, compute* EtaNormal, int* BC_List_Dir, compute* BC_Value_Dir, int nDir);
void EqSystem_solve(EqSystem* EqSystem);
void EqSystem_check(EqSystem* EqSystem);
int  pardisoSolveAssymmetric(int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n);
int  pardisoSolveSymmetric(int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n);

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
void compileShaders	(GLuint *ShaderProgram, const char* pVSFileName, const char* pFSFileName);
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

#endif /* STOKES_H_ */




