/*
 * Physics.h
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */

#ifndef PHYSICS_H_
#define PHYSICS_H_

#include "stokes.h"

// Physics
// =========================
//typedef struct SinglePhase SinglePhase;
struct SinglePhase
{
    int phase;
    compute weight;
    SinglePhase *next;
};
struct Physics
{
    compute R;
    // Physics Stokes
    compute g[2]; // gravity acceleration
    compute dt, dtAdv, dtT, dtDarcy;
    compute *Vx, *Vy, *P;

#if (CRANK_NICHOLSON_VEL || INERTIA)
    compute *Vx0, *Vy0;
#if (CRANK_NICHOLSON_P)
    compute *P0;
#endif
#endif

    compute maxVx, maxVy;
    compute *eta;

    //compute *eta0
    //compute *n;
    compute *rho;
//compute *rho0_g; // Density*norm_g

#if (HEAT)

    compute *k;           // Thermal conductivity
    compute *T, *T0, *DT; // temperature stored on cell centers
#endif

    compute epsRef; // reference strainrate

    int *phase;

//compute *Plitho;

#if (DARCY)

    compute *Pc, *DeltaP0, *DDeltaP; // old compaction pressure
    compute *phi, *Dphi, *phi0;      // fluid phase fraction
    compute *Pf;

    compute *divV0;

    compute *perm0_eta_f, *perm_eta_f; // permeability/eta_f
    compute minPerm;
    compute *eta_b; // bulk viscosity
    //compute *B; // elastic bulk modulus

    compute eta_f, rho_f; // viscosity of the fluid
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
    compute *sigma_xx_0, *sigma_xy_0;   // old stresses
    compute *Dsigma_xx_0, *Dsigma_xy_0; // stress corrections for markers
    compute *G;                         // shear modulus

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
    compute *strain;
    compute *Dstrain;
#endif

};

// Physics
// =========================
void Physics_Memory_allocate(Physics *Physics, Grid *Grid);
void Physics_Memory_free(Physics *Physics, Grid *Grid);
void Physics_P_initToLithostatic(Physics *Physics, Grid *Grid);
void Physics_Velocity_advectEulerian(Grid *Grid, Physics *Physics, BC *BCStokes, Numbering *NumStokes);
void Physics_Velocity_retrieveFromSolution(Physics *Physics, Grid *Grid, BC *BC, Numbering *Numbering, EqSystem *EqSystem, Numerics *Numerics);
#if (CRANK_NICHOLSON_VEL || INERTIA)
void Physics_VelOld_POld_updateGlobal(Physics *Physics, Grid *Grid);
#endif
void Physics_P_retrieveFromSolution(Physics *Physics, Grid *Grid, BC *BC, Numbering *Numbering, EqSystem *EqSystem, Numerics *Numerics);
void Physics_T_retrieveFromSolution(Physics *Physics, Grid *Grid, BC *BC, Numbering *Numbering, EqSystem *EqSystem, Numerics *Numerics);
void Physics_Eta_init(Physics *Physics, Grid *Grid, MatProps *MatProps, Numerics *Numerics);
void Physics_Eta_updateGlobal(Physics *Physics, Grid *Grid, Numerics *Numerics, BC *BCStokes, MatProps *MatProps);
void Physics_Dsigma_updateGlobal(Physics *Physics, Grid *Grid, BC *BC, Numbering *NumStokes, EqSystem *EqStokes, Numerics *Numerics);
void Physics_dt_update(Physics *Physics, Grid *Grid, MatProps *MatProps, Numerics *Numerics);
void Physics_StrainRateInvariant_getLocalCell(Physics *Physics, Grid *Grid, int ix, int iy, compute *EII);
void Physics_StrainRateInvariant_getLocalNode(Physics *Physics, BC *BCStokes, Grid *Grid, int ix, int iy, compute *EII);
void Physics_StressInvariant_getLocalCell(Physics *Physics, Grid *Grid, int ix, int iy, compute *SII);
#if (DARCY)
void Physics_Perm_updateGlobal(Physics *Physics, Grid *Grid, Numerics *Numerics, MatProps *MatProps);
void Physics_Phi_updateGlobal(Physics *Physics, Grid *Grid, Numerics *Numerics);
#endif
void Physics_Rho_updateGlobal(Physics *Physics, Grid *Grid, MatProps *MatProps);
void Physics_Phase_updateGlobal(Physics *Physics, Grid *Grid, Particles *Particles, MatProps *MatProps, BC *BCStokes);

void Physics_PhaseList_reinit(Physics *Physics, Grid *Grid);

void Physics_check(Physics *Physics, Grid *Grid, Char *Char);

void Physics_CellVal_retrieveFromSolution(compute *Val, int ISub, Grid *Grid, BC *BC, Numbering *Numbering, EqSystem *EqSystem);
void Physics_CellVal_SideValues_getFromBC_Global(compute *ECValues, Grid *Grid, BC *BC, Numbering *Numbering);
compute Physics_CellVal_SideValues_getFromBC_Local(compute neighValue, BC *BC, int IBC, int ix, int iy, Grid *Grid);
void Physics_CellVal_SideValues_copyNeighbours_Global(compute *ECValues, Grid *Grid);
void Physics_CellVal_SideValues_copyNeighbours_Global_i(int *ECValues, Grid *Grid);

void Physics_Phase_addSingle(SinglePhase** pointerToHead, int phase);

#endif /* PHYSICS_H_ */
