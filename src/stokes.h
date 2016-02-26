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
#define DEBUG false
#define VISU true
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




// Physics
// =========================
typedef struct Physics Physics;
struct Physics
{
	compute *Vx, *Vy, *P;
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
	int nxC, nyC, nCTot;
	int nxVx, nyVx, nxVy, nyVy, nxS, nyS;
	int nVxTot, nVyTot;
	coord xmin, xmax, ymin, ymax;
	compute dx, dy;

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
	int *oldCellId, *newCellId;
	int *linkNext, *linkHead;
};




// Visualization
// ========================
typedef struct Visu Visu;
struct Visu
{
	int ntri, ntrivert;
	GLuint* elements;
	GLfloat* U;
	GLfloat* vertices;
	GLfloat scale;
	GLuint VAO, VBO, CBO, EBO;
	GLuint ShaderProgram;
	const char* VertexShaderFile;
	const char* FragmentShaderFile;
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




// Particles
// =========================
void Particles_initCoord		(Grid* Grid, Particles* Particles);
void Particles_initPhase		(Grid* Grid, Particles* Particles);
void Particles_updateLinkedList (Grid* Grid, Particles* Particles);
void Particles_getPhysicsFrom	(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps);




// Physics
// =========================
void Physics_interpFromParticlesToCell(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps);
void Physics_interpFromCellToNode	  (Grid* Grid, compute* CellValue, compute* NodeValue);
void Physics_set_VxVyP_FromSolution(Physics* Physics, Grid* Grid, BC* BC, Numbering* Numbering, compute* sol);




// Visualization
// =========================
void Visu_allocateMemory	(Visu* Visu, Grid* Grid );
void Visu_freeMemory		(Visu* Visu );
void Visu_init				(Visu* Visu, Grid* Grid);
void Visu_initWindow		(GLFWwindow** window);
void Visu_initOpenGL		(Visu* Visu, Grid* Grid);
void error_callback			(int error, const char* description);
void key_callback			(GLFWwindow* window, int key, int scancode, int action, int mods);

void Visu_updateCenterValue (Visu* Visu, Grid* Grid, compute* Value);
void Visu_StrainRate		(Visu* Visu, Grid* Grid, Physics* Physics);



// Boundary conditions
// =========================
void BC_set(BC* BC, Grid* Grid, EqSystem* EqSystem, Physics* Physics);




// Numbering
// =========================
void Numbering_initMapAndSparseTripletIJ(BC* BC, Grid* Grid, EqSystem* EqSystem, Numbering* Numbering);
void Numbering_getLocalVx(int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingVx* LocVx, bool useNumMap);
void Numbering_getLocalVy(int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingVy* LocVy, bool useNumMap);
void Numbering_getLocalP (int ix, int iy, Numbering* Numbering, Grid* Grid, LocalNumberingP*  LocP , bool useNumMap);




// Equation system
// =========================
void EqSystem_allocateI(EqSystem* EqSystem);
void EqSystem_allocateMemory(EqSystem* EqSystem);
void EqSystem_freeMemory(EqSystem* EqSystem);
void EqSystem_assemble(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering);
void fill_J_V_local(int Type, int ix, int iy,int I, int iEq, int *Jsparse, compute* Vsparse, compute* RHS, int nxC, int nyC, compute dx, compute dy, int *NumMap, compute* EtaShear, compute* EtaNormal, int* BC_List_Dir, compute* BC_Value_Dir, int nDir);
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

#endif /* STOKES_H_ */




