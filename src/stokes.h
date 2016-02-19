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
#define DEBUG FALSE
#define NB_PHASE_MAX 10
#define NXC 10
#define NYC 10










//============================================================================//
//============================================================================//
//                                                                            //
//                                 TYPE DEFS                              	  //
//                                                                            //
//============================================================================//
//============================================================================//
// Define basic types
// =========================
typedef double coord;
typedef double compute;

// Physics
// =========================
typedef struct Physics Physics;
struct Physics
{
	compute *eta; // Viscosity
	compute *rho; // Density
	// compute stressOld
};

// Grid
// =========================
typedef struct Grid Grid;
struct Grid
{
	int nxC, nyC, nCTot;
	int nxVx, nyVx, nxVy, nyVy, nxS, nyS, nxN, nyN;
	int nVxTot, nVyTot;
	coord xmin, xmax, ymin, ymax;
	coord dx, dy;

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

	// Link list
	int* oldCellId, newCellId;
	int* linkNext, linkHead;
};

// Visualization
// ========================
typedef struct Visu Visu;
struct Visu
{
	int ntri, ntrivert;
	GLuint* elements;
	compute* U;
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
void allocateMemory(Grid* Grid, Particles* Particles, Physics* Physics);
void freeMemory(Particles* Particles, Physics* Physics);

// Particles
// =========================
void Particles_initCoord(Grid* Grid, Particles* Particles);
void Particles_initPhase(Grid* Grid, Particles* Particles);
void Particles_updateLinkedList(Grid* Grid, Particles* Particles);
void getPhysicsFromParticles2Grid(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps);


// Visualization
// =========================

void allocateVisuMemory( Visu* Visu );
void freeVisuMemory( Visu* Visu );
void initVisualization(Visu* Visu, Grid* Grid);


#endif /* STOKES_H_ */




