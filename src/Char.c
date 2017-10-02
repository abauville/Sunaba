/*
 * char.c
 *
 *  Created on: Feb 29, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Char_nonDimensionalize(Model* Model)
{
	Char* Char = &(Model->Char);
	Grid* Grid = &(Model->Grid);
	Physics* Physics = &(Model->Physics);
	MatProps* MatProps = &(Model->MatProps);
	BC* BCStokes = &(Model->BCStokes);
	BC* BCThermal = &(Model->BCThermal);
	IC* ICThermal = &(Model->ICThermal);
	IC* ICDarcy = &(Model->ICDarcy);
	Numerics* Numerics = &(Model->Numerics);
	Particles* Particles = &(Model->Particles);
	Output* Output = &(Model->Output);


	// SI units
	compute s 	= Char->time;			// [s]
	compute m 	= Char->length; 		// [m]
	compute kg 	= Char->mass; 			// [kg]
	compute K 	= Char->temperature; 	// [K]

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s); 		// Watt

	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	compute mol = 1.0;

	int i;


	Char->velocity 		= m/s; 			// [m.s-1]
	Char->density  		= kg/m/m/m; 	// [kg.m^-3]
	Char->stress  		= Pa;			// [Pa] or [kg.m^-1.s^-2]
	Char->viscosity 	= Pas; 			// [Pa.s] or [kg.m^-1.s-1]
	Char->acceleration 	= m/s/s; 		// [m.s^-2]
	Char->strainrate 	= 1.0/s; 		// [s^-1]




	// Grid
	// ======================
	Grid->xmin 		/= m;
	Grid->xmax 		/= m;
	Grid->ymin 		/= m;
	Grid->ymax 		/= m;

// Physics
	// ======================
	Physics->dt		/= s;
	Physics->g[0] 	/= m/(s*s);
	Physics->g[1] 	/= m/(s*s);


	Physics->epsRef /= 1.0/s;


	Physics->Cp 	/= J/kg/K;

#if (HEAT)
	Physics->R 		= 8.3144598;
	Physics->R 		/= J/K/mol;
#else
	Physics->R 		= 1.0;
#endif

	Physics->Pback /= Pa;

#if (DARCY)
	Physics->eta_f /= Pas;
	Physics->rho_f /= kg/(m*m*m);
#endif

#if (DARCY)
	Physics->y_oceanSurface /= m;
#endif


	// Material properties
	// ======================
	for (i = 0; i < MatProps->nPhase; ++i) {
		MatProps->rho0 [i] 			/= kg/(m*m*m);
		MatProps->alpha[i]  		/= 1.0/K;
		MatProps->beta [i]  		/= 1.0/Pa;
		MatProps->k    [i]  		/= W/m/K;
		MatProps->G    [i]  		/= Pa;
		MatProps->cohesion[i] 		/= Pa;
		MatProps->frictionAngle[i]  /= 1.0;

		MatProps->perm0[i] 			/= m*m;

		MatProps->vDiff[i].B 		/= 1.0/Pas;
		MatProps->vDiff[i].E 		/= J/mol;
		MatProps->vDiff[i].V 		/= (m*m*m)/mol;

		MatProps->vDisl[i].B 		/= pow(Pa,-MatProps->vDisl[i].n) / s;
		MatProps->vDisl[i].E 		/= J/mol;
		MatProps->vDisl[i].V 		/= (m*m*m)/mol;

		MatProps->vPei [i].B 		/= 1.0/s;
		MatProps->vPei [i].E 		/= J/mol;
		MatProps->vPei [i].V 		/= (m*m*m)/mol;
		MatProps->vPei [i].tau 		/= Pa;

		//printf("MatProps->phiIni[%i] = %.2e\n", i, MatProps->phiIni[i]);

#if (DARCY)
		MatProps->perm0_eta_f[i] = MatProps->perm0[i]/Physics->eta_f;
#endif
	}

	//printf("MatProps->vDisl[0] = %.2e, MatProps->vDisl[1] = %.2e, Pa = %.2e, s = %.2e, -MatProps->vDisl[0].n = %.2e, -MatProps->vDisl[1].n = %.2e, pow(Pa,-MatProps->vDisl[0].n) = %.2e, pow(Pa,-MatProps->vDisl[1].n) = %.2e \n", MatProps->vDisl[0].B, MatProps->vDisl[1].B, Pa, s, -MatProps->vDisl[0].n, -MatProps->vDisl[1].n, pow(Pa,-MatProps->vDisl[0].n), pow(Pa,-MatProps->vDisl[1].n));
	
	BCStokes->backStrainRate /= 1.0/s;
	BCStokes->refValue 		 /= m/s;
	BCStokes->DeltaL 		 /= m;
	BCStokes->Sandbox_TopSeg00 /= m;
	BCStokes->Sandbox_TopSeg01 /= m;

	BCThermal->TT 		/= K;
	BCThermal->TB 		/= K;
	BCThermal->refValue /= K;
	BCThermal->DeltaL   /= m;


	switch (ICThermal->SetupType) {
	case IC_HSC: // Half-space cooling model
		ICThermal->data[0] /= K; // noise
		ICThermal->data[1] /= K; // Tmantle
		ICThermal->data[2] /= s; // lithosphere age (in seconds)
		break;
	case IC_Gaussian: // Gaussian cooling model
		ICThermal->data[0] /= K; // noise
		ICThermal->data[1] /= K; // background
		ICThermal->data[2] /= K; // Amplitude
		ICThermal->data[3] /= m; // xc
		ICThermal->data[4] /= m; // yc
		ICThermal->data[5] /= m; // wx
		ICThermal->data[6] /= m; // wy
		break;
	default:
		printf("error: Unknown ICThermal->SetupType %i\n", ICThermal->SetupType);
		exit(0);
		break;
	}

	switch (ICDarcy->SetupType) {
		case IC_HSC: // Half-space cooling model
			printf("error, IC_HSC is not applicable to Darcy\n");
			exit(0);
			break;
		case IC_Gaussian: // Gaussian, applied on phi
			ICDarcy->data[0] /= 1.0; // noise
			ICDarcy->data[1] /= 1.0; // background
			ICDarcy->data[2] /= 1.0; // Amplitude
			ICDarcy->data[3] /= m; // xc
			ICDarcy->data[4] /= m; // yc
			ICDarcy->data[5] /= m; // wx
			ICDarcy->data[6] /= m; // wy
			break;
		default:
			printf("error: Unknown ICDarcy->SetupType %i\n", ICDarcy->SetupType);
			exit(0);
			break;
		}


	Numerics->StickyAirStress = 1.0*MatProps->cohesion[Physics->phaseRef]/1.0;


	Particles->passiveDx /= m;
	Particles->passiveDy /= m;


	Numerics->stickyAirSwitchingDepth /= m;
	Numerics->stickyAirTimeSwitchPassive /= s;
	Numerics->stickyAirTimeSinceLastPassiveSwitch = 0.0;

	Numerics->maxTime /= s;

	Output->timeFrequency /= s;

	Numerics->dtVep /= s;
	Numerics->dtMin /= s;
	Numerics->dtMax /= s;

}



void Char_reDimensionalize(Model* Model)
{
	Char* Char = &(Model->Char);
	Grid* Grid = &(Model->Grid);
	Physics* Physics = &(Model->Physics);
	MatProps* MatProps = &(Model->MatProps);
	BC* BCStokes = &(Model->BCStokes);
	BC* BCThermal = &(Model->BCThermal);
	IC* ICThermal = &(Model->ICThermal);
	IC* ICDarcy = &(Model->ICDarcy);
	Numerics* Numerics = &(Model->Numerics);
	Particles* Particles = &(Model->Particles);
	Output* Output = &(Model->Output);


	// SI units
	compute s 	= Char->time;			// [s]
	compute m 	= Char->length; 		// [m]
	compute kg 	= Char->mass; 			// [kg]
	compute K 	= Char->temperature; 	// [K]

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s); 		// Watt

	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	compute mol = 1.0;

	int i;


	Char->velocity 		= m/s; 			// [m.s-1]
	Char->density  		= kg/m/m/m; 	// [kg.m^-3]
	Char->stress  		= Pa;			// [Pa] or [kg.m^-1.s^-2]
	Char->viscosity 	= Pas; 			// [Pa.s] or [kg.m^-1.s-1]
	Char->acceleration 	= m/s/s; 		// [m.s^-2]
	Char->strainrate 	= 1.0/s; 		// [s^-1]




	// Grid
	// ======================
	Grid->xmin 		*= m;
	Grid->xmax 		*= m;
	Grid->ymin 		*= m;
	Grid->ymax 		*= m;

// Physics
	// ======================
	Physics->dt		*= s;
	Physics->g[0] 	*= m/(s*s);
	Physics->g[1] 	*= m/(s*s);


	Physics->epsRef *= 1.0/s;


	Physics->Cp 	*= J/kg/K;

#if (HEAT)
	Physics->R 		= 8.3144598;
	Physics->R 		*= J/K/mol;
#else
	Physics->R 		= 1.0;
#endif

	Physics->Pback *= Pa;

#if (DARCY)
	Physics->eta_f *= Pas;
	Physics->rho_f *= kg/(m*m*m);
#endif

#if (DARCY)
	Physics->y_oceanSurface *= m;
#endif


	// Material properties
	// ======================
	for (i = 0; i < MatProps->nPhase; ++i) {
		MatProps->rho0 [i] 			*= kg/(m*m*m);
		MatProps->alpha[i]  		*= 1.0/K;
		MatProps->beta [i]  		*= 1.0/Pa;
		MatProps->k    [i]  		*= W/m/K;
		MatProps->G    [i]  		*= Pa;
		MatProps->cohesion[i] 		*= Pa;
		MatProps->frictionAngle[i]  *= 1.0;

		MatProps->perm0[i] 			*= m*m;

		MatProps->vDiff[i].B 		*= 1.0/Pas;
		MatProps->vDiff[i].E 		*= J/mol;
		MatProps->vDiff[i].V 		*= (m*m*m)/mol;

		MatProps->vDisl[i].B 		*= pow(Pa,-MatProps->vDisl[i].n) / s;
		MatProps->vDisl[i].E 		*= J/mol;
		MatProps->vDisl[i].V 		*= (m*m*m)/mol;

		MatProps->vPei [i].B 		*= 1.0/s;
		MatProps->vPei [i].E 		*= J/mol;
		MatProps->vPei [i].V 		*= (m*m*m)/mol;
		MatProps->vPei [i].tau 		*= Pa;


#if (DARCY)
		MatProps->perm0_eta_f[i] = MatProps->perm0[i]/Physics->eta_f;
#endif
	}

	BCStokes->backStrainRate *= 1.0/s;
	BCStokes->refValue 		 *= m/s;
	BCStokes->DeltaL 		 *= m;
	BCStokes->Sandbox_TopSeg00 *= m;
	BCStokes->Sandbox_TopSeg01 *= m;

	BCThermal->TT 		*= K;
	BCThermal->TB 		*= K;
	BCThermal->refValue *= K;
	BCThermal->DeltaL   *= m;


	switch (ICThermal->SetupType) {
	case IC_HSC: // Half-space cooling model
		ICThermal->data[0] *= K; // noise
		ICThermal->data[1] *= K; // Tmantle
		ICThermal->data[2] *= s; // lithosphere age (in seconds)
		break;
	case IC_Gaussian: // Gaussian cooling model
		ICThermal->data[0] *= K; // noise
		ICThermal->data[1] *= K; // background
		ICThermal->data[2] *= K; // Amplitude
		ICThermal->data[3] *= m; // xc
		ICThermal->data[4] *= m; // yc
		ICThermal->data[5] *= m; // wx
		ICThermal->data[6] *= m; // wy
		break;
	default:
		printf("error: Unknown ICThermal->SetupType %i\n", ICThermal->SetupType);
		exit(0);
		break;
	}

	switch (ICDarcy->SetupType) {
		case IC_HSC: // Half-space cooling model
			printf("error, IC_HSC is not applicable to Darcy\n");
			exit(0);
			break;
		case IC_Gaussian: // Gaussian, applied on phi
			ICDarcy->data[0] *= 1.0; // noise
			ICDarcy->data[1] *= 1.0; // background
			ICDarcy->data[2] *= 1.0; // Amplitude
			ICDarcy->data[3] *= m; // xc
			ICDarcy->data[4] *= m; // yc
			ICDarcy->data[5] *= m; // wx
			ICDarcy->data[6] *= m; // wy
			break;
		default:
			printf("error: Unknown ICDarcy->SetupType %i\n", ICDarcy->SetupType);
			exit(0);
			break;
		}


	Numerics->StickyAirStress = 1.0*MatProps->cohesion[Physics->phaseRef]/1.0;


	Particles->passiveDx *= m;
	Particles->passiveDy *= m;


	Numerics->stickyAirSwitchingDepth *= m;
	Numerics->stickyAirTimeSwitchPassive *= s;
	Numerics->stickyAirTimeSinceLastPassiveSwitch = 0.0;

	Numerics->maxTime *= s;

	Output->timeFrequency *= s;

	Numerics->dtVep *= s;
	Numerics->dtMin *= s;
	Numerics->dtMax *= s;

}


void Char_rescale(Model* Model) {

	
	Grid* Grid = &(Model->Grid);
	Physics* Physics = &(Model->Physics);
	MatProps* MatProps = &(Model->MatProps);
	BC* BCStokes = &(Model->BCStokes);
	BC* BCThermal = &(Model->BCThermal);
	IC* ICThermal = &(Model->ICThermal);
	IC* ICDarcy = &(Model->ICDarcy);
	Numerics* Numerics = &(Model->Numerics);
	Particles* Particles = &(Model->Particles);
	Output* Output = &(Model->Output);


	// Create a copy of the current scale
	Char Char0 = Model->Char;
	Char CharN;
	printf("koko0\n");

	// Change scale
	Char_reDimensionalize(Model);
	Model->Char.time = Physics->dt*Model->Char.time;
	Char_nonDimensionalize(Model);
	CharN = Model->Char; // shortcut for practicality


	printf("Char->time = %.2e, oldChar->time = %.2e, Char0.velocity = %.2e, Char.velocity = %.2e, Char0.viscosity = %.2e, Char.viscosity = %.2e, Physics->dt = %.2e\n", CharN.time, Char0.time, Char0.velocity, CharN.velocity, Char0.viscosity, CharN.viscosity, Physics->dt);

	
	// Rescale Vx, Vy, P, sigmaOld
	int iVx, iVy, iCell, iNode;
	
	for (iVx=0;iVx<Grid->nVxTot;++iVx) {
		Physics->Vx[iVx] = Physics->Vx[iVx]*Char0.velocity / CharN.velocity;
#if (INERTIA)
		Physics->Vx0[iVx] = Physics->Vx0[iVx]*Char0.velocity / CharN.velocity;
#endif
	}

	for (iVy=0;iVy<Grid->nVyTot;++iVy) {
		Physics->Vy[iVy] = Physics->Vy[iVy]*Char0.velocity / CharN.velocity;
#if (INERTIA)
		Physics->Vy0[iVy] = Physics->Vy0[iVy]*Char0.velocity / CharN.velocity;
#endif
	}

	
	for (iCell=0;iCell<Grid->nECTot;++iCell) {
		Physics->eta[iCell] = Physics->eta[iCell]*Char0.viscosity / CharN.viscosity;
		compute Zold = Physics->Z  [iCell] ;
		Physics->Z  [iCell] = Physics->Z  [iCell]*Char0.viscosity / CharN.viscosity;
		Physics->khi[iCell] = Physics->khi[iCell]*Char0.viscosity / CharN.viscosity;
		if (iCell == 10) {
			//printf("Zold = %.2e, Z = %.2e, Char0.viscosity =%.2e, CharN.viscosity =%.2e\n", Zold,  Physics->Z, Char0.viscosity, CharN.viscosity);
			printf("Char0.viscosity =%.2e, CharN.viscosity =%.2e\n", Char0.viscosity, CharN.viscosity);
			printf("Zold = %.2e, Z = %.2e\n", Zold,  Physics->Z[iCell]);
		}
	}

	for (iNode=0;iNode<Grid->nSTot;++iNode) {
		Physics->etaShear[iNode] = Physics->etaShear[iNode]*Char0.viscosity / CharN.viscosity;
		Physics->ZShear  [iNode] = Physics->ZShear  [iNode]*Char0.viscosity / CharN.viscosity;
		Physics->khiShear[iCell] = Physics->khiShear[iCell]*Char0.viscosity / CharN.viscosity;
	}





	//note0: Strain rates are not stored, they are going be rescaled automatically
	//note1: characteristic stress unaffected by a change in char time.

	/*
	for (iCell=0;iCell<Grid->nECTot;++iCell) {
		Physics->P[iCell] = Physics->P[iVy]*Char0.stress / CharN.stress;
#if (USE_SIGMA0_OV_G)
		Physics->sigma_xx_0_ov_G[iCell] = Physics->sigma_xx0_ov_G[iCell]*Char0.stress / CharN.stress;
#else
		Physics->sigma_xx_0[iCell] = Physics->sigma_xx0[iCell]*Char0.stress / CharN.stress;
#endif
	}

	for (iNode=0;iNode<Grid->nSTot;++iNode) {
		#if (USE_SIGMA0_OV_G)
		Physics->sigma_xy_0_ov_G[iNode] = Physics->sigma_xy0_ov_G[iNode]*Char0.stress / CharN.stress;
#else
		Physics->sigma_xy_0[iNode] = Physics->sigma_xy0[iNode]*Char0.stress / CharN.stress;
#endif
	}
	*/
	






}