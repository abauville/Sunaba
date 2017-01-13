/*
 * char.c
 *
 *  Created on: Feb 29, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Char_nonDimensionalize(Char* Char, Grid* Grid, Physics* Physics, MatProps* MatProps, BC* BCStokes, BC* BCThermal, IC* ICThermal, IC* ICDarcy, Numerics* Numerics, Particles* Particles)
{
	// SI units
	compute s 	= Char->time;			// second
	compute m 	= Char->length; 		// meter
	compute kg 	= Char->mass; 			// kilogram
	compute K 	= Char->temperature; 	// Kelvin

	// Other units
	compute J = kg*m*m/(s*s); 			// Joule
	compute W = kg*m*m/(s*s*s); 		// Watt

	compute Pa  = kg/m/s/s; 			// Pascal
	compute Pas = kg/m/s; 				// Poise, Pa.s

	compute mol = 1.0;

	int i;

	Numerics->StickyAirStress = 5e6 / Pa;

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

	compute norm_g = sqrt(Physics->g[0]*Physics->g[0] + Physics->g[1]*Physics->g[1]);

	Physics->gFac[0] 	= Physics->g[0]/norm_g;
	Physics->gFac[1] 	= Physics->g[1]/norm_g;

	Physics->epsRef /= 1.0/s;


	Physics->Cp 	/= J/kg/K;


	Physics->R 		= 8.3144598;
	Physics->R 		/= J/K/mol;

#if (DARCY)
	Physics->eta_f /= Pas;
	Physics->rho_f /= kg/(m*m*m);
	Physics->rho_f_g = Physics->rho_f*norm_g;
#endif

#if (DARCY)
	Physics->y_oceanSurface /= m;
#endif
printf("MatProps->vDisl[0] = %.2e, MatProps->vDisl[1] = %.2e\n", MatProps->vDisl[0].B, MatProps->vDisl[1].B);
printf("MatProps->vDiff[0] = %.2e, MatProps->vDiff[1] = %.2e\n", MatProps->vDiff[0].B, MatProps->vDiff[1].B);
	// Material properties
	// ======================
	for (i = 0; i < MatProps->nPhase; ++i) {
		MatProps->rho0 [i] 			/= kg/(m*m*m);
		MatProps->rho0_g [i] 		 = MatProps->rho0 [i] * norm_g;
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


#if (DARCY)
		MatProps->perm0_eta_f[i] = MatProps->perm0[i]/Physics->eta_f;
#endif
		//MatProps->eta_b[i] 	/= Pas;
		//MatProps->B	   [i] 	/= Pa;
	}

	printf("MatProps->vDisl[0] = %.2e, MatProps->vDisl[1] = %.2e, Pa = %.2e, s = %.2e, -MatProps->vDisl[0].n = %.2e, -MatProps->vDisl[1].n = %.2e, pow(Pa,-MatProps->vDisl[0].n) = %.2e, pow(Pa,-MatProps->vDisl[1].n) = %.2e \n", MatProps->vDisl[0].B, MatProps->vDisl[1].B, Pa, s, -MatProps->vDisl[0].n, -MatProps->vDisl[1].n, pow(Pa,-MatProps->vDisl[0].n), pow(Pa,-MatProps->vDisl[1].n));
	printf("MatProps->vDiff[0] = %.2e, MatProps->vDiff[1] = %.2e, Pa = %.2e, s = %.2e, -MatProps->vDisl[0].n = %.2e, -MatProps->vDisl[1].n = %.2e, pow(Pa,-MatProps->vDisl[0].n) = %.2e, pow(Pa,-MatProps->vDisl[1].n) = %.2e \n", MatProps->vDiff[0].B, MatProps->vDiff[1].B, Pa, s, -MatProps->vDisl[0].n, -MatProps->vDisl[1].n, pow(Pa,-MatProps->vDisl[0].n), pow(Pa,-MatProps->vDisl[1].n));
	BCStokes->backStrainRate /= 1.0/s;
	BCStokes->refValue 		 /= m/s;
	BCStokes->DeltaL 		 /= m;

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





	Particles->passiveDx /= m;
	Particles->passiveDy /= m;


	Numerics->stickyAirSwitchingDepth /= m;
	Numerics->stickyAirTimeSwitchPassive /= s;
	Numerics->stickyAirTimeSinceLastPassiveSwitch = 0.0;

}



