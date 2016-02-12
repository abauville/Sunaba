/*
 * StokesFD.h
 *
 *  Created on: Feb 9, 2016
 *      Author: abauville
 */

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


// Useful macros
// =========================
#define DEBUG FALSE
#define NB_PHASE_MAX 10
#define NXC 10
#define NYC 10


// Define basic types
// =========================
typedef double coord;
typedef double compute;


// Structures
// ==========================

// Physics
typedef struct Physics Physics;
struct Physics
{
	//compute *eta;
	// compute stressOld
};

typedef struct Grid Grid;
struct Grid
{
	int nxC, nyC;
	int nxVx, nyVx, nxVy, nyVy, nxS, nyS, nxN, nyN;
	int nVxTot, nVyTot;
	coord xmin, xmax, ymin, ymax;
	coord dx, dy;
};

typedef struct MatProps MatProps;
struct MatProps
{
	int nPhase;
	compute rho[NB_PHASE_MAX], eta0[NB_PHASE_MAX];
};

typedef struct Particles Particles;
struct Particles
{
	int nPartPerCell;
	coord *x, *y;
	int *phase; // i is the index of the cell in which the particle is
};

// Prototypes
// ==========================















#endif /* STOKES_H_ */




