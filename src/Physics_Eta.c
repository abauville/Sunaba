/*
 * Physics_Eta.c
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */


#include "stokes.h"

#define USE_INVETA_EP false
#define COMPUTE_SHEAR_VISCOSITY false

void Physics_Eta_init(Model* Model) 
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);


	int iy, ix, iCell;
	SinglePhase* thisPhaseInfo;
	compute P, T;
	int phase;
	compute EII, weight;
	compute B, E, V, n, taup, q, s, gamma;
	compute invEtaDiff, invEtaDisl, invEtaPei;
	compute R = Physics->R;
	compute eta, G, cohesion, frictionAngle, eta_thisPhase;
	compute sumOfWeights;

	// =======================================================
	// Initial viscosity
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			//Physics->etaVisc[iCell] = Physics->eta0[iCell];
			//Physics->eta[iCell] = Physics->eta0[iCell];

			thisPhaseInfo = Physics->phaseListHead[iCell];

			EII = fabs(Physics->epsRef)/1.0;
#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];

			eta = 0.0;
			G = 0.0;
			cohesion = 0.0;
			frictionAngle = 0.0;
			compute invEta_EP = 0.0;
			while (thisPhaseInfo != NULL) {
				invEtaDiff = 0.0;
				invEtaDisl = 0.0;
				invEtaPei  = 0.0;
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				G 				+= weight/MatProps->G[phase];
				//G 				+= weight*MatProps->G[phase];
				//G 				+= log10(MatProps->G[phase])*weight;
				//G 				+= weight/log10(MatProps->G[phase]);
				cohesion 		+= MatProps->cohesion[phase] * weight;
				frictionAngle 	+= MatProps->frictionAngle[phase] * weight;
				if (MatProps->vDiff[phase].isActive) {
					B 			 = MatProps->vDiff[phase].B;
					E 			 = MatProps->vDiff[phase].E;
					V 			 = MatProps->vDiff[phase].V;
					invEtaDiff   = (2.0*(B*exp( - (E+V*P)/(R*T)   )));
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					invEtaDisl 	 = (2.0*pow(B*exp( - (E+V*P)/(R*T)   ),1.0/n)*pow(EII,-1.0/n+1.0));
				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					invEtaPei 	 = (2.0*pow(B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) ) ,1.0/s)*pow(EII,-1.0/s+1.0) );
				}
				thisPhaseInfo 	= thisPhaseInfo->next;
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
				invEta_EP += (1.0/(MatProps->G[phase]*Physics->dt) + 1.0/eta_thisPhase) * weight;
				eta += weight * eta_thisPhase;





			}
			eta = eta / sumOfWeights;
			//eta = pow(10.0,eta / sumOfWeights);
			invEta_EP /= sumOfWeights;
			/*
			if (eta>Numerics->etaMax) {
				eta = Numerics->etaMax;
			}
			if (eta<Numerics->etaMin) {
				eta = Numerics->etaMin;
			}
			*/

			Physics->eta[iCell] = eta;

			//Physics->G[iCell]  = G/Physics->sumOfWeightsCells[iCell]/G;
			Physics->G[iCell]  = Physics->sumOfWeightsCells[iCell]/G;
			//Physics->G[iCell]  = G/Physics->sumOfWeightsCells[iCell];
			//Physics->G[iCell]  = pow(10.0,G/Physics->sumOfWeightsCells[iCell]);
			//Physics->G[iCell]  = pow(10.0,Physics->sumOfWeightsCells[iCell]/G);
			Physics->khi[iCell] = 1E30;

			//Physics->Z[iCell] = 1.0/( 1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );
#if (USE_INVETA_EP)
			Physics->Z[iCell] = 1.0/( invEta_EP) ;
#else
			Physics->Z[iCell] = 1.0/( 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );
#endif

			if (Physics->Z[iCell]<Numerics->etaMin) {
				Physics->Z[iCell] = Numerics->etaMin;
			}
#if (DARCY)
			Physics->eta_b[iCell] = Physics->eta[iCell]/(Physics->phi[iCell]);
			Physics->khi_b[iCell] = 1e30;
			Physics->Zb[iCell] 	  = 1.0/( 1.0/Physics->eta_b[iCell] + 1.0/(Physics->G[iCell]/(sqrt(Physics->phi[iCell]))*Physics->dt) );
#endif

		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);
#if (DARCY)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi_b, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta_b, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Zb, Grid);
#endif


	int iNode;
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->GShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->G,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}

}






void Physics_Eta_Simple_updateGlobal(Model* Model)
{
	
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	int iNode, iy, ix;


	// ===== get EffStrainRate =====
	Physics_Eta_EffStrainRate_updateGlobal (Model);
	// ===== get EffStrainRate =====


	// ===== get the Z as a visco-elastic predictor =====
	Physics_Eta_VEpredictor_updateGlobalCell(Model);

	// ================================================================================
	// 									Shear nodes viscosity
	//#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->GShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->G,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}
	// 									Shear nodes viscosity
	// ================================================================================

}





void Physics_Eta_FromParticles_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);



	compute locX, locY;
	int ix, iy;

	INIT_PARTICLE
	int iCell;

	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->eta	[iCell] = 0.0;
		Physics->G		[iCell] = 0.0;
		Physics->khi	[iCell] = 0.0;
		Physics->Z		[iCell] = 0.0;
	}


	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		Physics->etaShear [iNode] = 0.0;
		Physics->GShear	  [iNode] = 0.0;
		Physics->khiShear [iNode] = 0.0;
		Physics->ZShear   [iNode] = 0.0;
	}

	
	compute EII;
	compute* EIICell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* SII0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exx = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* dVxdyGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* dVydxGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* Rotxy = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute dVxdy, dVydx, dVxdx, dVydy;
	compute sq_sigma_xy0, sigma_xx0;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Exx[iCell] = 0.5*(dVxdx-dVydy);

			Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);
			
			EIICell[iCell] = EII;

			sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
			sigma_xx0  = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			SII0Cell[iCell] = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx, Grid);

	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			
			dVxdyGrid[iNode] =  dVxdy;
			dVydxGrid[iNode] =  dVydx;

			Exy[iNode] = 0.5*(dVxdy+dVydx);
			Rotxy[iNode] = 0.5*(dVxdy-dVydx);
		}
	}


	compute ExyPart, ExxPart;
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			// Loop through the particles in the shifted cell
			// ======================================
			while (thisParticle!=NULL) {

				locX = Particles_getLocX(ix, thisParticle->x,Grid);
				locY = Particles_getLocY(iy, thisParticle->y,Grid);

				int IxN[4], IyN[4];
				IxN[0] =  0;  	IyN[0] =  0; // lower left
				IxN[1] =  1;	IyN[1] =  0; // lower right
				IxN[2] =  0; 	IyN[2] =  1; // upper left
				IxN[3] =  1; 	IyN[3] =  1; // upper right
				// ===== weight cells =====

				int i;
				if 		 	(locX>=0 && locY>=0) { // upper right
					i = 3;
				} else if 	(locX<0 && locY>=0) { // upper left
					// the particle is in the SE quadrant, the cell center 1 is NW (wrt to the node ix,iy)
					i = 2;
				} else if 	(locX>=0 && locY<0) { // lower right
					i = 1;
				} else if 	(locX<0 && locY<0) { // lower left
					i = 0;
				} else {
					printf("error in Interp_ECVal_Cell2Particle_Local. No case was triggered\n.");
					exit(0);
				}
				iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
				compute weightCell = fabs(locX)*fabs(locY);
				// ===== weight cells =====

				// ===== weight nodes =====
				compute weightNode = (1.0 - fabs(locX)) * (1.0 - fabs(locY));
				// ===== weight nodes =====

				int phase = thisParticle->phase;



				ExxPart = Interp_ECVal_Cell2Particle_Local(Exx, ix, iy, Grid->nxEC, locX, locY);
				ExyPart = Interp_NodeVal_Node2Particle_Local(Exy, ix, iy, Grid->nxS, Grid->nyS, locX, locY);				
				EII = sqrt(ExxPart*ExxPart + ExyPart*ExyPart);

				compute eta;

				compute T = 1.0;

				compute invEtaDiff = 0.0;
				compute invEtaDisl = 0.0;
				compute invEtaPei = 0.0;
				
				compute BDiff, BDisl, BPei;
				compute B, E, V, n, gamma, taup, q, s;
				compute R = 1.0;
				compute P = 0.0;	
				
				
				if (MatProps->vDiff[phase].isActive) {
					B 			 = MatProps->vDiff[phase].B;
					E 			 = MatProps->vDiff[phase].E;
					V 			 = MatProps->vDiff[phase].V;
					BDiff = B*exp( - (E+V*P)/(R*T)   );
					invEtaDiff   = (2.0*(BDiff));
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					BDisl = B*exp( - (E+V*P)/(R*T)   );
					invEtaDisl 	 = (2.0*pow(BDisl,1.0/n)*pow(EII,-1.0/n+1.0));
				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					BPei	 = B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) );
					invEtaPei 	 = (2.0*pow(BPei ,1.0/s)*pow(EII,-1.0/s+1.0) );
				}
				eta = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
				compute G = MatProps->G[phase];
				
				compute cohesion = MatProps->cohesion[phase];
				compute frictionAngle = MatProps->frictionAngle[phase];
				compute dt = Physics->dt;

				compute Z = 1.0/(1.0/eta + 1.0/(G*dt));
				
				compute Sxx0 = thisParticle->sigma_xx_0;
				compute Sxy0 = thisParticle->sigma_xy_0;
				//compute SII0 = sqrt(Sxx0*Sxx0 + Sxy0*Sxy0);

/*
#if (USE_UPPER_CONVECTED)
				compute RotxyPart = Interp_NodeVal_Node2Particle_Local(Rotxy, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				int nxEC = Grid->nxEC;
				int nxS = Grid->nxS;
				int nyS = Grid->nyS;

				compute sqEII_Part = ExxPart*ExxPart + ExyPart*ExyPart;
				compute sqSII0_Part = Sxx0*Sxx0 + Sxy0*Sxy0;
				


				compute dVxdyPart = Interp_NodeVal_Node2Particle_Local(dVxdyGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				compute dVydxPart = Interp_NodeVal_Node2Particle_Local(dVydxGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				//compute Eff_strainRate = 1.0/(2.0*G*dt) * sqrt(pow((2.0*ExxPart*G*dt + Sxx0 + 2.0*dt*(Sxx0*ExxPart + Sxy0*dVxdyPart)),2.0) + pow((2.0*ExyPart*G*dt - Sxx0*dt*2.0*RotxyPart+ Sxy0),2.0));
				compute Eff_strainRate = sqrt(sqEII_Part + ExxPart*Sxx0/(G*dt) + ExyPart*Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sqSII0_Part    );
#else
				compute Eff_strainRate = sqrt(EII*EII + ExxPart*Sxx0/(G*dt) + ExyPart*Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*SII0*SII0   );
#endif
*/
				compute Eff_strainRate = Interp_ECVal_Cell2Particle_Local(Physics->EII_eff,ix,iy,Grid->nxEC,locX,locY);
				

				compute sigmaII = 2.0*Z*Eff_strainRate;
				compute khi;
				compute phi = 0.0;
				compute Pe = Interp_ECVal_Cell2Particle_Local(Physics->P,ix,iy,Grid->nxEC, locX, locY);
				Pe = fmax(Pe,0.0);
				compute sigma_y = cohesion*cos(frictionAngle) + Pe*sin(frictionAngle);
				
				if (sigmaII > sigma_y) {

					khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );

					if (khi<0.0) {
						// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
						// if it happens then I consider the case where there are == , which means khi -> inf
						printf("khi = %.2e, eta = %.2e, G = %.2e, dt = %.2e, Eff_Strainrate = %.2e, 1-phi = %.2e, sigma_y = %.2e, Pe = %.2e, Pmin = %.2e\n", khi, eta, G, dt, Eff_strainRate, 1.0-phi, sigma_y, Pe, -cohesion*cos(frictionAngle)/sin(frictionAngle));
						printf("WTF!\n");
						khi = 1e30;
						//exit(0);
					}
#if (USE_INVETA_EP)
					Z 	= (1.0-phi)*1.0/(1.0/khi + invEta_EP);
#else
					Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
#endif
					sigmaII = 2.0*Z*Eff_strainRate;

				} else {
					khi = 1e30;
				}				



				Physics->eta				[iCell] += eta * weightCell;
				Physics->G					[iCell] += weightCell / G;
				Physics->khi				[iCell] += khi * weightCell;
				//Physics->khi				[iCell] += weightCell / khi;
				Physics->Z					[iCell] += Z * weightCell;
				
				Physics->etaShear 			[iNode] += eta * weightNode;
				Physics->GShear 			[iNode] += weightNode / G;
				Physics->khiShear			[iNode] += khi * weightNode;
				//Physics->khiShear			[iNode] += weightNode / khi;
				Physics->ZShear				[iNode] += Z   * weightNode;
				
				thisParticle = thisParticle->next;
			}
		}
	}


	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
			//printf("sumOfWeights[%i] = %.2e\n", iCell, Physics->sumOfWeightsCells	[iCell]);

		if (Physics->sumOfWeightsCells	[iCell]==0.0) {
			printf("error in Interp_All_Particles2Grid_Global. Cell #%i received no contribution from particles (i.e. empty cell).\n", iCell);
			exit(0);
		}


		Physics->eta	[iCell] /= Physics->sumOfWeightsCells	[iCell];
		Physics->G		[iCell]  = Physics->sumOfWeightsCells	[iCell] / Physics->G		[iCell];
		Physics->khi	[iCell] /= Physics->sumOfWeightsCells	[iCell];
		//Physics->khi	[iCell]  = Physics->sumOfWeightsCells	[iCell] / Physics->khi	[iCell];
		//Physics->Z		[iCell] /= Physics->sumOfWeightsCells	[iCell];
		compute dt = Physics->dt;
		compute phi = 0.0;
		compute eta = Physics->eta	[iCell];
		compute khi = Physics->khi	[iCell];
		compute G   = Physics->G	[iCell];
		Physics->Z		[iCell] = (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);

	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		Physics->etaShear [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic average
		Physics->GShear	  [iNode]  = Physics->sumOfWeightsNodes[iNode] / Physics->GShear	  [iNode]; // arithmetic average
		Physics->khiShear [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic average
		//Physics->khiShear [iNode]  = Physics->sumOfWeightsNodes[iNode] / Physics->khiShear [iNode]; // arithmetic average
		//Physics->ZShear   [iNode] /= Physics->sumOfWeightsNodes[iNode]; // arithmetic average
		compute dt = Physics->dt;
		compute phi = 0.0;
		compute eta = Physics->etaShear	[iNode];
		compute khi = Physics->khiShear	[iNode];
		compute G   = Physics->GShear	[iNode];
		Physics->ZShear		[iNode] = (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
	}


	//END_PARTICLES


	free(Exx);
	free(Exy);
	free(EIICell);
	free(SII0Cell);
	free(Rotxy);
	free(dVxdyGrid);
	free(dVydxGrid);
	
}



void Physics_Eta_EffStrainRate_updateGlobal(Model* Model) {
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	Numerics* Numerics 		= &(Model->Numerics);

	int ix, iy;

	compute dVxdy, dVydx, dVxdx, dVydy;

	compute* Exx_VE_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy_VE_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute dt = Physics->dt;
	int iCell;
	#pragma omp parallel for private(iy,ix, iCell, dVxdx, dVydy) OMP_SCHEDULE
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			compute G = Physics->G[iCell];
			Exx_VE_CellGlobal[iCell] = 0.5*(dVxdx-dVydy) + Physics->sigma_xx_0[iCell]/(2.0*G*dt);
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx_VE_CellGlobal, Grid);

	int iNode;
	#pragma omp parallel for private(iy,ix, iNode, dVxdy, dVydx) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;

			compute G = Physics->GShear[iNode];
			Exy_VE_NodeGlobal[iNode] = 0.5*(dVxdy+dVydx) + Physics->sigma_xy_0[iNode]/(2.0*G*dt);
		}
	}

	if (Numerics->invariantComputationType==0) {
		#pragma omp parallel for private(iy,ix, iCell) OMP_SCHEDULE
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			for (ix = 1; ix<Grid->nxEC-1; ix++) {
				iCell = ix + iy*Grid->nxEC;

				compute Exx_VE_sq = Exx_VE_CellGlobal[iCell]*Exx_VE_CellGlobal[iCell];
				//compute Exy_VE_sq = Interp_Product_NodeVal_Node2Cell_Local(Exy_VE_NodeGlobal , Exy_VE_NodeGlobal, ix, iy, Grid->nxS);
				compute Exy_VE = Interp_NodeVal_Node2Cell_Local(Exy_VE_NodeGlobal, ix, iy, Grid->nxS);
				compute Exy_VE_sq = Exy_VE * Exy_VE;

				Physics->EII_eff[iCell] = sqrt(Exx_VE_sq + Exy_VE_sq);

			}
		}
		Physics_CellVal_SideValues_copyNeighbours_Global(Physics->EII_eff, Grid);

		#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
		for (iy = 0; iy<Grid->nyS; iy++) {
			for (ix = 0; ix<Grid->nxS; ix++) {
				iNode = ix + iy*Grid->nxS;

				compute Exy_VE_sq = Exy_VE_NodeGlobal[iNode] * Exy_VE_NodeGlobal[iNode];
				//compute Exx_VE_sq = Interp_Product_ECVal_Cell2Node_Local(Exx_VE_CellGlobal,Exx_VE_CellGlobal,ix,iy,Grid->nxEC);
				compute Exx_VE = Interp_ECVal_Cell2Node_Local(Exx_VE_CellGlobal, ix, iy, Grid->nxEC);
				compute Exx_VE_sq = Exx_VE*Exx_VE;			

				Physics->EII_effShear[iNode] = sqrt(Exx_VE_sq + Exy_VE_sq);

			}
		}
	} else if (Numerics->invariantComputationType==1) {
		#pragma omp parallel for private(iy,ix, iCell) OMP_SCHEDULE
		for (iy = 1; iy<Grid->nyEC-1; iy++) {
			for (ix = 1; ix<Grid->nxEC-1; ix++) {
				iCell = ix + iy*Grid->nxEC;

				compute Exx_VE_sq = Exx_VE_CellGlobal[iCell]*Exx_VE_CellGlobal[iCell];
				compute Exy_VE_sq = Interp_Product_NodeVal_Node2Cell_Local(Exy_VE_NodeGlobal , Exy_VE_NodeGlobal, ix, iy, Grid->nxS);
				//compute Exy_VE = Interp_NodeVal_Node2Cell_Local(Exy_VE_NodeGlobal, ix, iy, Grid->nxS);
				//compute Exy_VE_sq = Exy_VE * Exy_VE;

				Physics->EII_eff[iCell] = sqrt(Exx_VE_sq + Exy_VE_sq);

			}
		}
		Physics_CellVal_SideValues_copyNeighbours_Global(Physics->EII_eff, Grid);

		#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
		for (iy = 0; iy<Grid->nyS; iy++) {
			for (ix = 0; ix<Grid->nxS; ix++) {
				iNode = ix + iy*Grid->nxS;

				compute Exy_VE_sq = Exy_VE_NodeGlobal[iNode] * Exy_VE_NodeGlobal[iNode];
				compute Exx_VE_sq = Interp_Product_ECVal_Cell2Node_Local(Exx_VE_CellGlobal,Exx_VE_CellGlobal,ix,iy,Grid->nxEC);
				//compute Exx_VE = Interp_ECVal_Cell2Node_Local(Exx_VE_CellGlobal, ix, iy, Grid->nxEC);
				//compute Exx_VE_sq = Exx_VE*Exx_VE;			

				Physics->EII_effShear[iNode] = sqrt(Exx_VE_sq + Exy_VE_sq);

			}
		}
	} else {
		printf("error: unknwon Numerics->invariantComputationType %i\n", Numerics->invariantComputationType);
		exit(0);
	}

	


	free(Exx_VE_CellGlobal);
	free(Exy_VE_NodeGlobal);

}


void Physics_Eta_VEpredictor_updateGlobalCell(Model* Model) {
	// Update Physics->Z and Physics->eta
	// according to the Visco-elastic predictor

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);


	int iy, ix, iCell;

	compute T, P, phi;
	compute alpha;

	compute eta, G;
	SinglePhase* thisPhaseInfo;

	compute invEtaDiff, invEtaDisl, invEtaPei;
	int phase;
	compute weight;
	compute B, E, V, n, taup, gamma, q, s;
	compute R = Physics->R;
	compute ZUpper, ZLower;
	compute BDiff[NB_PHASE_MAX], BDisl[NB_PHASE_MAX], BPei[NB_PHASE_MAX];
	compute EII;

	compute eta_thisPhase;


	compute dt = Physics->dt;

	compute Z, Zcorr, PrevZcorr;

	compute sigmaII;

	compute tol = 1e-7;

	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;


			Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);


#if (HEAT)
			P 	= Physics->P[iCell];
			T 	= Physics->T[iCell];
#else
			T = 1.0;
			P = 0.0;
#endif
#if (DARCY)
			phi = Physics->phi[iCell];
#else
			phi = 0.0;
#endif
			compute sumOfWeights 	= Physics->sumOfWeightsCells[iCell];


			alpha = 1.0;
			compute invEta_EP = 0.0;

			// Precompute B and viscosities using EII
			eta = 0.0;
			G = 0.0;
			compute maxInvVisc = 0.0;
			
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				invEtaDiff = 0.0;
				invEtaDisl = 0.0;
				invEtaPei  = 0.0;
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				G 				+= weight/MatProps->G[phase];

				if (MatProps->vDiff[phase].isActive) {
					B 			 = MatProps->vDiff[phase].B;
					E 			 = MatProps->vDiff[phase].E;
					V 			 = MatProps->vDiff[phase].V;
					BDiff[phase] = B*exp( - (E+V*P)/(R*T)   );
					invEtaDiff   = (2.0*(BDiff[phase]));
					maxInvVisc = fmax(invEtaDiff,maxInvVisc);
				}
				if (MatProps->vDisl[phase].isActive) {
					B 			 = MatProps->vDisl[phase].B;
					E 			 = MatProps->vDisl[phase].E;
					V 			 = MatProps->vDisl[phase].V;
					n 			 = MatProps->vDisl[phase].n;
					BDisl[phase] = B*exp( - (E+V*P)/(R*T)   );
					invEtaDisl 	 = (2.0*pow(BDisl[phase],1.0/n)*pow(EII,-1.0/n+1.0));
					maxInvVisc = fmax(invEtaDisl,maxInvVisc);
				}
				if (MatProps->vPei[phase].isActive) {
					B 			 = MatProps->vPei[phase].B;
					E 			 = MatProps->vPei[phase].E;
					V 			 = MatProps->vPei[phase].V;
					gamma 		 = MatProps->vPei[phase].gamma;
					taup  		 = MatProps->vPei[phase].tau;
					q 			 = MatProps->vPei[phase].q;
					s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
					BPei[phase]	 = B*pow(gamma*taup,-s)*exp( - (E+V*P)/(R*T) * pow((1.0-gamma),q) );
					invEtaPei 	 = (2.0*pow(BPei[phase] ,1.0/s)*pow(EII,-1.0/s+1.0) );
					maxInvVisc = fmax(invEtaPei,maxInvVisc);
				}
				
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));

				eta += weight * eta_thisPhase;
				invEta_EP += (1.0/(MatProps->G[phase]*dt)+1.0/eta_thisPhase) * weight;
				thisPhaseInfo 	= thisPhaseInfo->next;
			}
			G 				 = sumOfWeights	/ G;
			
			eta 			/= sumOfWeights;

			invEta_EP /= sumOfWeights;

			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}

			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);


			Z = 0.5*(ZUpper+ZLower);
			Zcorr = Z;



			sigmaII = 2.0*Z*Physics->EII_eff[iCell];

			// compute viscosities using sigmaII
			while (fabs(Zcorr/Z)>tol) {
				eta = 0.0;
				thisPhaseInfo = Physics->phaseListHead[iCell];
				invEta_EP = 0.0;

				while (thisPhaseInfo != NULL) {
					invEtaDiff = 0.0;
					invEtaDisl = 0.0;
					invEtaPei  = 0.0;
					phase = thisPhaseInfo->phase;
					weight = thisPhaseInfo->weight;
					if (MatProps->vDiff[phase].isActive) {
						invEtaDiff 	= (2.0*(BDiff[phase]));
					}
					if (MatProps->vDisl[phase].isActive) {
						n 			 = MatProps->vDisl[phase].n;
						invEtaDisl 	 = (2.0*BDisl[phase]*pow(sigmaII,-1.0+n));
					}
					if (MatProps->vPei[phase].isActive) {
						E 			 = MatProps->vPei[phase].E;
						V 			 = MatProps->vPei[phase].V;
						gamma 		 = MatProps->vPei[phase].gamma;
						q 			 = MatProps->vPei[phase].q;
						s   		 = (E+V*P)/(R*T)*pow((1.0-gamma),(q-1.0))*q*gamma;
						invEtaPei  	= ( 2.0*BPei[phase]*pow(sigmaII,-1.0+s) );
					}

					eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));
					eta += weight * eta_thisPhase;
					thisPhaseInfo 	= thisPhaseInfo->next;

					invEta_EP += (1.0/(MatProps->G[phase]*dt)+1.0/eta_thisPhase) * weight;


				}

				eta 			/= sumOfWeights;
				//eta = 1e30;
				//eta = pow(10.0,eta / sumOfWeights);
				//invEta_EP = pow(10.0,invEta_EP / sumOfWeights);
				invEta_EP = invEta_EP / sumOfWeights;
				PrevZcorr = Zcorr;

				//Zcorr = (1.0-phi)*(1.0/(1.0/(G*dt) + 1.0/eta)) - Z;
				Zcorr = (1.0-phi)*1.0/invEta_EP - Z;

				
				
				if (Zcorr/PrevZcorr<-0.9) {
					alpha = alpha/2.0;
				}
				Z += alpha*Zcorr;

				sigmaII = 2.0*Z*Physics->EII_eff[iCell];
			}


			Physics->Z[iCell] = Z;
			Physics->eta[iCell] = eta;
			Physics->G[iCell] = G;

		}
	}

	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);


	
}



void Physics_Eta_computeLambda_FromParticles_updateGlobal(Model* Model, bool updateStresses) {

	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	//Numerics* Numerics 		= &(Model->Numerics);

	compute locX, locY;
	int ix, iy;

	INIT_PARTICLE




	compute* Exx_Grid = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy_Grid = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute dVxdy, dVydx, dVxdx, dVydy;
	int iCell;
	#pragma omp parallel for private(iy,ix, iCell, dVxdx, dVydy) OMP_SCHEDULE
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;
			Exx_Grid[iCell] = 0.5*(dVxdx-dVydy);
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx_Grid, Grid);

#pragma omp parallel for private(iy,ix, iNode, dVxdy, dVydx) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			
			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;
			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Exy_Grid[iNode] = 0.5*(dVxdy+dVydx);
		}
	}

	compute* sumOfWeightsCells = (compute*) malloc(Grid->nECTot*sizeof(compute));
	compute* sumOfWeightsNodes = (compute*) malloc(Grid->nSTot *sizeof(compute));
	#pragma omp parallel for private(iy,ix, iCell) OMP_SCHEDULE
	for(iCell=0; iCell<Grid->nECTot;++iCell) {
		Physics->Lambda[iCell] = 0.0;
		Physics->khi[iCell] = 0.0;
		sumOfWeightsCells[iCell] = 0.0;
	}
	#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for(iNode=0; iNode<Grid->nSTot;++iNode) {
		Physics->LambdaShear[iNode] = 0.0;
		sumOfWeightsNodes[iNode] = 0.0;
	}

	

	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  1;	IyN[1] =  0; // lower right
	IxN[2] =  0; 	IyN[2] =  1; // upper left
	IxN[3] =  1; 	IyN[3] =  1; // upper right

	// Loop through particles and compute lambda
	compute dt = Physics->dt;
	

	int iColor; // indexing of the color group for nodes. Nodes of the same color don't collide with each other. i.e. similar to matrix coloring
	int ixStart[4] = {0,0,1,1};
	int iyStart[4] = {0,1,0,1};

	for (iColor = 0; iColor < 4; ++iColor) {
#pragma omp parallel for private(ix, iy, iNode, thisParticle, locX, locY) OMP_SCHEDULE
		for (iy = iyStart[iColor]; iy < Grid->nyS; iy+=2) { // Gives better result not to give contribution from the boundaries
			for (ix = ixStart[iColor]; ix < Grid->nxS; ix+=2) { // I don't get why though

	//for (iy = 0; iy < Grid->nyS; ++iy) {
	//	for (ix = 0; ix < Grid->nxS; ++ix) {
				iNode = ix  + (iy  )*Grid->nxS;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle!=NULL) {
					locX = Particles_getLocX(ix, thisParticle->x,Grid);
					locY = Particles_getLocY(iy, thisParticle->y,Grid);

					if (fabs(locX)>1.0 || fabs(locY)>1.0 ) {
						printf("Error locXY, locX = %.1f, locY = %.1f\n", locX, locY);
						exit(0);
					}

					compute Lambda, lambda, khi;
					int phase = thisParticle->phase;
					if (phase == Physics->phaseAir || phase == Physics->phaseWater) {
							// First part of the correction of stresses on the particles: add subgrid (adding remaining will be done in a second step)
							Lambda = 1.0;
							khi = 1e30;
					} else {
						compute Exx = Interp_ECVal_Cell2Particle_Local(Exx_Grid, ix, iy, Grid->nxEC, locX, locY);
						compute Exy = Interp_NodeVal_Node2Particle_Local(Exy_Grid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
						compute Z = Interp_ECVal_Cell2Particle_Local(Physics->Z, ix, iy, Grid->nxEC, locX, locY);
						compute Pe = Interp_ECVal_Cell2Particle_Local(Physics->P, ix, iy, Grid->nxEC, locX, locY);
						// Fail safe
						if (Pe<0.0) {
							Pe = 0.0;
						}
						compute G = MatProps->G[phase];
						compute cohesion = MatProps->cohesion[phase];
						compute fAngle = MatProps->frictionAngle[phase];
						
						compute Txx0 = thisParticle->sigma_xx_0;
						compute Txy0 = thisParticle->sigma_xy_0;
						
						compute Ty = cohesion*cos(fAngle) + Pe*sin(fAngle);

						compute Exx_eff = (Exx + Txx0/(2.0*G*dt));
						compute Exy_eff = (Exy + Txy0/(2.0*G*dt));
						compute EII_eff = sqrt(Exx_eff*Exx_eff + Exy_eff*Exy_eff);
						compute TII_VE = 2.0*Z*EII_eff;

						if (TII_VE>Ty) {
							//printf("ExxPart = %.2e, ExxGrid0 = %.2e, ExxGrid1 = %.2e, ExxGrid2 = %.2e, ExxGrid3 = %.2e\n", Exx, Exx_Grid[ix+(iy)*Grid->nxEC], Exx_Grid[ix+1+(iy)*Grid->nxEC], Exx_Grid[ix+(iy+1)*Grid->nxEC], Exx_Grid[ix+1+(iy+1)*Grid->nxEC]);
							Lambda = Ty/TII_VE;
							lambda = 2.0*EII_eff*(1.0-Lambda);
							khi = Ty/lambda;
						} else {
							Lambda = 1.0;
							khi = 1e30;
						}

						if (updateStresses) {
							thisParticle->sigma_xx_0 = 2.0*Z*(Exx_eff) * Lambda;
							thisParticle->sigma_xy_0 = 2.0*Z*(Exy_eff) * Lambda;
						}
					}
					int signX, signY;
					if (locX<0.0) {
						signX = -1;
					} else {
						signX = 1;
					}
					if (locY<0.0) {
						signY = -1;
					} else {
						signY = 1;
					}
					int i;
					if 		 	(signX>=0 && signY>=0) { // upper right
						i = 3;
					} else if 	(signX<0 && signY>=0) { // upper left
						// the particle is in the SE quadrant, the cell center 1 is NW (wrt to the node ix,iy)
						i = 2;
					} else if 	(signX>=0 && signY<0) { // lower right
						i = 1;
					} else if 	(signX<0 && signY<0) { // lower left
						i = 0;
					} else {
						printf("error in Interp_ECVal_Cell2Particle_Local. No case was triggered\n.");
						exit(0);
					}
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);
					compute weight = fabs(locX)*fabs(locY);
					sumOfWeightsCells[iCell] += weight;
					Physics->Lambda[iCell] += Lambda*weight;
					Physics->khi[iCell] += khi*weight;

					weight = (1.0 - fabs(locX)) * (1.0 - fabs(locY));
					sumOfWeightsNodes[iNode] += weight;
					Physics->LambdaShear[iNode] += Lambda*weight;
					
					thisParticle = thisParticle->next;
				} // while particles
			} // ix
		} // iy 
	} // iColor
	
	// Finished to compute the averages
	#pragma omp parallel for private(iy,ix, iCell) OMP_SCHEDULE
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			if (sumOfWeightsCells[iCell]==0.0) {
				printf("error: zero sum on cell ix = %i, iy = %i", ix, iy);
				exit(0);
			} 
			Physics->Lambda [iCell] /= sumOfWeightsCells[iCell];
			//Physics->Eps_pxx[iCell] /= sumOfWeightsCells[iCell];

			Physics->khi[iCell] /= sumOfWeightsCells[iCell];
		}	
	}
	//Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Eps_pxx, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Lambda, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->LambdaShear[iNode] /= sumOfWeightsNodes[iNode];
			//Physics->Eps_pxy[iNode] /= sumOfWeightsNodes[iNode];
		}
	}
	
	/*
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			compute lambda = Physics->lambda [iCell];

			if (lambda>0.0) {}
				compute Z = Physics->Z[iCell];
				compute G = Physics->G[iCell];
				compute dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;
				compute dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

				compute Exx = 0.5*(dVxdx-dVydy);
				//compute Exy = Interp_NodeVal_Node2Cell_Local(Eps_xy_NodeGlobal, ix, iy, nxS);

				compute Txx0 = Physics->sigma_xx_0[iCell];
				

				compute Txx_VE = 2.0 * Z*(Exx + Txx0/(2.0*G*dt));
				compute Txy_VE = 2.0 * Z*(Exx + Txy0/(2.0*G*dt));
				
				compute TII_VE = sqrt(Txx_VE*Txx_VE + Txy_VE*Txy_VE);
			}
		}	
	}
	*/
	


	free(Exx_Grid);
	free(Exy_Grid);
	free(sumOfWeightsCells);
	free(sumOfWeightsNodes);

}


void Physics_Eta_ZandLambda_updateGlobal(Model* Model) {

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics		= &(Model->Physics);
	MatProps* MatProps		= &(Model->MatProps);
	Numerics* Numerics		= &(Model->Numerics);

	Physics_Eta_EffStrainRate_updateGlobal(Model);


	compute* Ty_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));


	int ix, iy, iCell, iNode;

	int Method = Numerics->yieldComputationType;
	SinglePhase* thisPhaseInfo;
					// ===== Plastic stress corrector =====
#pragma omp parallel for private(iy,ix, iCell, thisPhaseInfo) OMP_SCHEDULE
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			// update cohesion and friction angle
			compute sumOfWeights 	= Physics->sumOfWeightsCells[iCell];
			int phase;
			compute weight;
			compute cohesion, frictionAngle;
			cohesion = 0.0;
			frictionAngle = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			compute staticPfFac = 0.0;
			compute staticPfFacWeakFac = 0.0;
			compute frictionAngleWeakFac = 0.0;
			compute cohesionWeakFac = 0.0;
			compute strainWeakStart = 0.0;
			compute strainWeakEnd = 0.0;
			while (thisPhaseInfo != NULL) {
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				cohesion 		+= MatProps->cohesion[phase] * weight;
				frictionAngle 	+= MatProps->frictionAngle[phase] * weight;

				staticPfFac 		+= MatProps->staticPfFac[phase] * weight;
				staticPfFacWeakFac 	+= MatProps->staticPfFacWeakFac[phase] * weight;
				frictionAngleWeakFac+= MatProps->frictionAngleWeakFac[phase] * weight;
				cohesionWeakFac		+= MatProps->cohesionWeakFac[phase] * weight;

				strainWeakStart 	+= MatProps->strainWeakStart[phase] * weight;
				strainWeakEnd 		+= MatProps->strainWeakEnd[phase] * weight;

				thisPhaseInfo = thisPhaseInfo->next;
			}
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;

			staticPfFac 		/= sumOfWeights;
			staticPfFacWeakFac 	/= sumOfWeights;
			frictionAngleWeakFac/= sumOfWeights;
			cohesionWeakFac		/= sumOfWeights;

			strainWeakStart 	/= sumOfWeights;
			strainWeakEnd 		/= sumOfWeights;
			
			
			// Strain weakening
			compute CriticalStrain0= strainWeakStart;
			compute CriticalStrain1 = strainWeakEnd;

			compute Cini = cohesion;
			compute Cend = cohesion*(1.0-cohesionWeakFac);

			compute fricIni = frictionAngle;
			compute fricEnd = frictionAngle*(1.0-frictionAngleWeakFac);

			compute staticPfFacIni = staticPfFac;
			//compute staticPfFacEnd = staticPfFac*(1.0-staticPfFacWeakFac);
			compute staticPfFacEnd = (1.0-staticPfFacWeakFac)*(  staticPfFac ) + staticPfFacWeakFac; // such that the weakening factor in front of the pressure is (1.0-staticPfFacWeakFac)*(1.0-Pf)
			staticPfFacEnd = fmin(staticPfFacEnd,0.99);
			//staticPfFacEnd = fmax(staticPfFacEnd,0.0);
			compute preFac = 0.05;
			compute Fac;
			compute plasticStrain = Physics->strain[iCell] + Physics->Dstrain[iCell];
			if (plasticStrain<CriticalStrain0) {
				Fac = (1.0 - preFac *  (plasticStrain)/(CriticalStrain0));
			} else {
				Fac = 1.0 -  preFac  - (plasticStrain-CriticalStrain0)/(CriticalStrain1-CriticalStrain0);
			}
			
			Fac = fmin(Fac,1.0);
			Fac = fmax(Fac,0.0);

			if (Physics->phase[iCell] == Physics->phaseAir || Physics->phase[iCell] == Physics->phaseWater) {
				Fac = 1.0;
			}

				
			/*
			compute staticPfFac = 0.0;
			
			if (Physics->phase[iCell]==1 || iy<=1) {
				FluidPFac = 0.8;
				//cohesion = 0.0;
			} else if (Physics->phase[iCell]==3) {
				FluidPFac = 0.5;
			} else if (Physics->phase[iCell]==4) {
				FluidPFac = 0.95;
			}
			*/
			

			if (ix>=Grid->nxEC-3) {
				frictionAngle = fricIni;
				cohesion = Cini;
			} else {
				cohesion = Cini*Fac + (1.0-Fac)*Cend;
				frictionAngle = fricIni*Fac + (1.0-Fac)*fricEnd;
				staticPfFac = staticPfFacIni*Fac + (1.0-Fac)*staticPfFacEnd;
			}


			
			if (iy<=1) {
				//cohesion = 0.0;
				//frictionAngle = fmin(frictionAngle,30.0*PI/180.0);
				//FluidPFac = 0.5;
			}
			
			compute Z_VE = 1.0/(1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );

			compute Pe = (1.0-staticPfFac) * Physics->P[iCell];
			
			if (Pe<0.0) {
				Pe = 0.0;
			}


			compute TII_VE;
			if (Method==0) {
				Physics->Lambda[iCell] = 1.0;
				TII_VE = Physics_StressInvariant_getLocalCell(Model, ix, iy);

			} else if (Method==1) {
				Physics->Lambda[iCell] = 1.0;
				//compute TII_VE = Physics_StressInvariant_getLocalCell(Model, ix, iy);
				compute EII_eff = Physics->EII_eff[iCell];
				TII_VE = 2.0 * Z_VE * EII_eff;

			} else {
				printf("error unknwon yieldComputationType %i, should be 0 or 1\n",Numerics->yieldComputationType);

			}


			compute Ty = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);

			Ty_CellGlobal[iCell] = Ty;

			if (TII_VE>Ty) {


				compute Lambda = Ty/TII_VE;
				compute lambda = 2.0*Physics->EII_eff[iCell]*(1.0-Lambda);

				if (Method==0) {
					Physics->Lambda[iCell] = Lambda;
					Physics->khi[iCell] = Ty/lambda;
				} else if (Method==1) {
					Physics->Z[iCell] = Z_VE * Lambda;
					Physics->Lambda[iCell] = 1.0;//Lambda;
					Physics->khi[iCell] = Ty/lambda;
				}

			} else {
				Physics->khi[iCell] = 1e30;
				if (Method==1) {
					Physics->Z[iCell] = Z_VE;
				}
				Physics->Lambda[iCell] = 1.0;
			}

			if (Physics->khi[iCell]<1e30 && Physics->phase[iCell] != Physics->phaseAir && Physics->phase[iCell] != Physics->phaseWater) {
				compute SII = Physics_StressInvariant_getLocalCell(Model, ix, iy);// //(Physics, Grid, ix, iy, &SII);
				Physics->Dstrain[iCell] = SII/(2.0*Physics->khi[iCell])*Physics->dtAdv; // Recovering the incremental plastic strain
			} else {
				Physics->Dstrain[iCell] = 0.0;
			}

		}
	}
	// ===== Plastic stress corrector =====
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Lambda, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Ty_CellGlobal, Grid);

	//int iNode;
#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->LambdaShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Lambda, ix, iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z, ix, iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi, ix, iy, Grid->nxEC);
			
			compute Z_VE = 1.0/(1.0/Physics->etaShear[iNode] + 1.0/(Physics->GShear[iNode]*Physics->dt) );

			compute TII_VE = 0;
			if (Method==0) {
				Physics->LambdaShear[iNode] = 1.0;
				TII_VE = Physics_StressInvariant_getLocalNode(Model, ix, iy);

			} else if (Method==1) {
				Physics->LambdaShear[iNode] = 1.0;
				//compute TII_VE = Physics_StressInvariant_getLocalCell(Model, ix, iy);
				compute EII_eff = Physics->EII_effShear[iNode];
				TII_VE = 2.0 * Z_VE * EII_eff;

			} else {
				printf("error in Physics_Eta_ZandLambda_updateGlobal Method can be only 0 or 1.\n");
			}


			compute Ty = Interp_ECVal_Cell2Node_Local(Ty_CellGlobal, ix, iy, Grid->nxEC);

			if (TII_VE>Ty) {
				compute Lambda = Ty/TII_VE;
				compute lambda = 2.0*Physics->EII_effShear[iNode]*(1.0-Lambda);

				if (Method==0) {
					Physics->LambdaShear[iNode] = Lambda;
					Physics->khiShear[iNode] = Ty/lambda;
				} else if (Method==1) {
					Physics->ZShear[iNode] = Z_VE * Lambda;
					Physics->LambdaShear[iNode] = 1.0;//Lambda;
					Physics->khiShear[iNode] = Ty/lambda;
				}

			} else {
				Physics->khiShear[iNode] = 1e30;
				if (Method==1) {
					Physics->ZShear[iNode] = Z_VE;
				}
				Physics->LambdaShear[iNode] = 1.0;
			}
			if ((Physics->ZShear[iNode]==0.0)) {
				printf("Zshear=0, Z_VE = %.2e, Physics->LambdaShear[iNode] = %.2e,Ty =%.2e, TII_VE =%.2e\n", Z_VE, Physics->LambdaShear[iNode], Ty, TII_VE);
			}


		}




	}

	free(Ty_CellGlobal);

}
