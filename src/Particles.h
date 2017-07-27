/*
 * Particles.h
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */

#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "stokes.h"
// Particles
// =========================

// Single Particle storing coordinate, temp and info for a linked list

//typedef struct SingleParticle SingleParticle;
struct SingleParticle
{
    coord x, y;
    int phase;
    float passive; // some passive attribute used for visualization

#if (HEAT)
    compute T;
#endif

    // Old stresses
    compute sigma_xx_0;
    compute sigma_xy_0;

#if (CRANK_NICHOLSON_VEL || INERTIA)
    compute Vx, Vy;
#endif

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
    SingleParticle *next;
};

// Id Changed
//typedef struct ParticlePointerList ParticlePointerList;
struct ParticlePointerList
{
    //int data;
    SingleParticle *pointer;
    ParticlePointerList *next;
};
// Particles, i.e. info of the system of all particles

typedef enum { PartPassive_Grid,
               PartPassive_Grid_w_Layers } ParticlePassiveGeom;
struct Particles
{
    int nPC, nPCX, nPCY; // number of particles per cell, tot, in x and in y
    int n;               // number of particles
    compute minPartPerCellFactor, maxPartPerCellFactor;
    coord noiseFactor;
    SingleParticle **linkHead;

    ParticlePassiveGeom passiveGeom;
    compute passiveDx, passiveDy;

    compute *dispAtBoundL, *dispAtBoundR; // ,*dispAtBoundT, *dispAtBoundB
    int *currentPassiveAtBoundL, *currentPassiveAtBoundR;

};

// Particles
// =========================
void Particles_Memory_allocate(Particles *Particles, Grid *Grid);
void Particles_Memory_free(Particles *Particles, Grid *Grid);
void Particles_initCoord(Particles *Particles, Grid *Grid);
void Particles_initPassive(Particles *Particles, Grid *Grid, Physics *Physics);
void Particles_updateLinkedList(Particles *Particles, Grid *Grid, Physics *Physics);
void Particles_injectOrDelete(Particles *Particles, Grid *Grid);
void Particles_injectAtTheBoundaries(Particles *Particles, Grid *Grid, Physics *Physics, MatProps *MatProps);
void Particles_advect(Particles *Particles, Grid *Grid, Physics *Physics);
void Particles_Periodicize(Particles *Particles, Grid *Grid);
void Particles_teleportInsideTheDomain(Particles *Particles, Grid *Grid, Physics *Physics);
void Particles_deleteIfOutsideTheDomain(Particles *Particles, Grid *Grid);
void Particles_switchStickyAir(Particles *Particles, Grid *Grid, Physics *Physics, Numerics *Numerics, MatProps *MatProps, BC *BCStokes);
void addToParticlePointerList(ParticlePointerList **pointerToHead, SingleParticle *thisParticle);
void freeParticlePointerList(ParticlePointerList *head);
void Particles_freeAllSingleParticles(Particles *Particles, Grid *Grid);
void addSingleParticle(SingleParticle **pointerToHead, SingleParticle *modelParticle);

#endif /* PARTICLES_H_ */
