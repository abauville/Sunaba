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
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}

}




void Physics_Eta_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);
	BC* BCStokes 			= &(Model->BCStokes);

	int iCell, iy, ix;

	// local copies of values
	compute eta, cohesion, frictionAngle;
	compute sigma_y, sigmaII;
	compute EII;
	compute phi = 0.0;


	compute dt = Physics->dt;

	compute Pe;
#if (DARCY)
	compute eta_b;

	compute Rad = 2.0; // radius of the griffith curve
	//compute phiMin = Numerics->phiMin;
	compute phiCrit = Numerics->phiCrit;
	compute sigmaT;//, PeSwitch;
	compute khi_b, Zb, Py;
	compute Bulk, divV, DeltaP0;
	compute DeltaP;
	//compute tol;
#endif

	compute sigmaII0;
	compute Z;
	compute sigma_xx0, sq_sigma_xy0;




	compute G;

	compute khi;



	compute Eff_strainRate;

	SinglePhase* thisPhaseInfo;
	int phase;
	compute weight, sumOfWeights;

	compute B, E, V, n, taup, gamma, q, s;


	compute P, T;


	compute invEtaDiff, invEtaDisl, invEtaPei;


	compute R = Physics->R;
	compute ZUpper, ZLower;


	compute alpha = 0.25;


	compute BDiff[NB_PHASE_MAX], BDisl[NB_PHASE_MAX], BPei[NB_PHASE_MAX];

	compute maxInvVisc;

	compute tol = 1e-7;
	compute PrevZcorr, Zcorr;
	compute eta_thisPhase;

	compute* ZprePlasticity = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* sigma_y_Stored = (compute*) malloc(Grid->nECTot * sizeof(compute));

#if (STRAIN_SOFTENING)
	compute strainReductionFac = .95; // 0.0 stays the same 1.0 = 100% reduction
#endif


	//compute sigma_y;
#if (!DARCY)
#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi) OMP_SCHEDULE collapse(2)
#else
#pragma omp parallel for private(iy,ix, iCell, sq_sigma_xy0, sigma_xx0, sigmaII0, EII, sumOfWeights, P, T, phi, alpha, eta, eta_thisPhase, G, maxInvVisc, cohesion, frictionAngle, thisPhaseInfo, phase, weight, B, E, V, n, gamma, taup, q, s, BDiff, BDisl, BPei,invEtaDiff, invEtaDisl, invEtaPei, ZUpper, ZLower, Z, Zcorr, Eff_strainRate, sigmaII, PrevZcorr, Pe, sigma_y, khi, sigmaT, Bulk, khi_b, eta_b, divV, DeltaP0, Zb, DeltaP, Py) OMP_SCHEDULE collapse(2)
#endif
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;

			//  Compute sigmaII0
			sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
			sigma_xx0  = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

			Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];



			compute dVxdy, dVydx, dVxdx, dVydy, Eps_xx;
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Eps_xx = 0.5*(dVxdx-dVydy);


			// Anton's trick
			dVxdy = 0;
			dVydx = 0;
			compute Exy_x_Sxy0 = 0.0;
			compute Exy_x_Sxy0_ov_G = 0.0;
			compute Exy = 0.0;
			compute dVxdy_av = 0.0;
			compute dVydx_av = 0.0;
			int iNode, Ix, Iy;
			int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
			int IyMod[4] = {0,0,1,1};
			for (iNode = 0; iNode < 4; ++iNode) {
				Ix = (ix-1)+IxMod[iNode];
				Iy = (iy-1)+IyMod[iNode];

				dVxdy = ( Physics->Vx[(Ix  )+(Iy+1)*Grid->nxVx]
									  - Physics->Vx[(Ix  )+(Iy  )*Grid->nxVx] )/Grid->dy;


				dVydx = ( Physics->Vy[(Ix+1)+(Iy  )*Grid->nxVy]
									  - Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy] )/Grid->dx;

				dVxdy_av += 0.25*dVxdy;
				dVydx_av += 0.25*dVydx;

				Exy += 0.25*(0.5*(dVxdy+dVydx));
#if (USE_SIGMA0_OV_G)
				Exy_x_Sxy0_ov_G += 0.25*(0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0_ov_G[Ix+Iy*Grid->nxS];
#else 
				Exy_x_Sxy0 += 0.25*(0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0[Ix+Iy*Grid->nxS];
#endif

			}















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
			alpha = 1.0;


			compute invEta_EP = 0.0;

			// Precompute B and viscosities using EII
			eta = 0.0;
			G = 0.0;
			maxInvVisc = 0.0;
			cohesion = 0.0;
			frictionAngle = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
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
				thisPhaseInfo 	= thisPhaseInfo->next;
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));


				eta += weight * eta_thisPhase;
				
				//invEta_EP += log10(1.0/(G*dt)+1.0/eta_thisPhase) * weight;
				//invEta_EP += (1.0/(G*dt)+1.0/eta_thisPhase) * weight;
				invEta_EP += (1.0/(MatProps->G[phase]*dt)+1.0/eta_thisPhase) * weight;


			}
			G 				 = sumOfWeights	/ G;
			//G 				 /= sumOfWeights	;
			//G  = pow(10.0,G/Physics->sumOfWeightsCells[iCell]);
			//G  = pow(10.0,Physics->sumOfWeightsCells[iCell]/G);
			eta 			/= sumOfWeights;
			//eta = pow(10.0,eta / sumOfWeights);
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;

			//invEta_EP = pow(10.0,invEta_EP / sumOfWeights);
			invEta_EP = invEta_EP / sumOfWeights;

#if (STRAIN_SOFTENING)
			compute strainLimit = 1.0;
			compute coeff = (Physics->strain[iCell]/strainLimit);
			coeff = fmin(1.0,coeff);
			
			//frictionAngle *= (1.0-coeff*strainReductionFac);
			/*
			if (ix>=Grid->nxEC-2) {
				frictionAngle = 15.0*PI/180.0;
			}
			*/
			if (Physics->phase[iCell]!=Physics->phaseAir) {
				cohesion *= (1.0-coeff*strainReductionFac);
				compute minCohesion = .1*1e6 / (Model->Char.mass/Model->Char.length/Model->Char.time/Model->Char.time);
				cohesion = fmax(cohesion, minCohesion);
				frictionAngle *= (1.0-coeff*.1);
			}
#endif



			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}
#if (USE_INVETA_EP)
			ZLower = 1.0/invEta_EP;
#else
			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);
#endif
			
			

			Z = 0.5*((1.0-phi)*ZUpper+(1.0-phi)*ZLower);
			Zcorr = Z;

			
#if (USE_SIGMA0_OV_G)
			compute sigma_xx0_ov_G  = Physics->sigma_xx_0_ov_G[iCell];// + Physics->Dsigma_xx_0[iCell];
			compute sq_sigma_xy0_ov_G;
			sq_sigma_xy0_ov_G  = Physics->sigma_xy_0_ov_G[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0_ov_G[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0_ov_G += Physics->sigma_xy_0_ov_G[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0_ov_G[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0_ov_G += Physics->sigma_xy_0_ov_G[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0_ov_G[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0_ov_G += Physics->sigma_xy_0_ov_G[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0_ov_G[ix  +(iy  )*Grid->nxS];
			
			compute sigmaII0_ov_G = sqrt((sigma_xx0_ov_G)*(sigma_xx0_ov_G)    + 0.25*(sq_sigma_xy0_ov_G));

			Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0_ov_G/dt + Exy_x_Sxy0_ov_G/(dt) + (1.0/(2.0*dt))*(1.0/(2.0*dt))*sigmaII0_ov_G*sigmaII0_ov_G   );
#else

#if (USE_UPPER_CONVECTED)
			
			// effective strain rate including the upper convected correction of stresses
			compute Exx = Eps_xx;
			compute Txx0 = Physics->sigma_xx_0[iCell];
			compute Txy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0,ix,iy,Grid->nxS);
			dVxdy = dVxdy_av;
			dVydx = dVydx_av;
			Eff_strainRate = 1.0/(2.0*G*dt) * sqrt(pow((2.0*Exx*G*dt + Txx0 + 2.0*dt*(Txx0*dVxdx + Txy0*dVxdy)),2.0) + pow((2.0*Exy*G*dt - Txx0*dt*(dVxdy - dVydx) + Txy0),2.0));
			//Eff_strainRate = 1.0/(2.0*G*dt) * sqrt(pow((2.0*Exx*G*dt + Txx0 + 2.0*dt*(Txx0*Exx + Txy0*Exy)),2.0) + pow((2.0*Exy*G*dt - Txx0*dt*(dVxdy - dVydx) + Txy0),2.0));
#else
			Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(G*dt) + Exy_x_Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
#endif
#endif
			sigmaII = 2.0*Z*Eff_strainRate;

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
					//invEta_EP += log10(1.0/(G*dt)+1.0/eta_thisPhase) * weight;
					//invEta_EP += (1.0/(G*dt)+1.0/eta_thisPhase) * weight;
					invEta_EP += (1.0/(MatProps->G[phase]*dt)+1.0/eta_thisPhase) * weight;


				}

				eta 			/= sumOfWeights;
				//eta = pow(10.0,eta / sumOfWeights);
				//invEta_EP = pow(10.0,invEta_EP / sumOfWeights);
				invEta_EP = invEta_EP / sumOfWeights;
				PrevZcorr = Zcorr;
#if (USE_INVETA_EP)
				Zcorr = (1.0-phi)*(1.0/(invEta_EP)) - Z;
#else
				Zcorr = (1.0-phi)*(1.0/(1.0/(G*dt) + 1.0/eta)) - Z;
#endif

				
				
				if (Zcorr/PrevZcorr<-0.9) {
					alpha = alpha/2.0;
				}
				Z += alpha*Zcorr;

				sigmaII = 2.0*Z*Eff_strainRate;
			}


			// Compute the effective Pressure Pe
#if (DARCY)
			// Limit the effective pressure
			if (phi>=phiCrit) {
				Bulk = G/sqrt(phi);
				khi_b = 1E30;
				eta_b = eta/phi;

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];

				Zb 	= (1.0-phi)*1.0/(1.0/eta_b + 1.0/(Bulk*dt));
				Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

			} else {
				Pe 		= Physics->P [iCell];
				khi_b = 1E30;
				eta_b = eta/phi;

				divV  = (  Physics->Vx[ix+iy*Grid->nxVx] - Physics->Vx[ix-1+ iy   *Grid->nxVx]  )/Grid->dx;
				divV += (  Physics->Vy[ix+iy*Grid->nxVy] - Physics->Vy[ix  +(iy-1)*Grid->nxVy]  )/Grid->dy;
				DeltaP0 = Physics->DeltaP0[iCell];

				Bulk = G/sqrt(phi);
				Zb 	= (1.0-phi)* 1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
			}


#else
			compute Pf = 0.0;
			/*
			compute rho = 1000.0/(Model->Char.mass/(Model->Char.length*Model->Char.length*Model->Char.length));
			compute g = fabs(Physics->g[1]);
			compute hSurf = Grid->ymax;//8000.0/Model->Char.length;
			compute h = (hSurf-Grid->dy*iy);
			compute Pf = rho*g*h;
			*/
			Pe 		= Physics->P [iCell] - Pf;

#endif

#if (USE_INVETA_EP)
			Z 	= (1.0-phi)*1.0/(invEta_EP);
#else
			Z 	= (1.0-phi)*1.0/(1.0/eta + 1.0/(G*dt));
#endif
			


			sigmaII = 2.0*Z*Eff_strainRate;

			ZprePlasticity[iCell] = Z/(1.0-phi);


			khi = 1e30;
#if (DARCY)
			khi_b = 1e30;
#endif
			alpha = 1.0;



			sigma_y = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);


#if (DARCY)
			sigmaT = (cohesion*cos(frictionAngle))/Rad;

			if (Pe<-sigmaT) {

				Py = -sigmaT;
				khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
				Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
				Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc

				sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle
			}
#else
			// Since there is no griffiths handling for negative pressure for the non darcy case yet
			// here I assume a flat Mohr Coulomb when Pe <0
			if (Pe<=0) {
				sigma_y = cohesion * cos(frictionAngle);
			}
#endif

			sigma_y_Stored[iCell] = sigma_y;

			compute Pe0;

			sigmaII0 = sigmaII;
			Pe0 = Pe;
#if (DARCY)
			Py = sigmaII - sigmaT;
			khi_b = 1e30;
#endif
			compute yieldTol = 1e-4;
			int iDum = 0;

			do {
				sigmaII = sigmaII0;

				iDum++;


				if (sigmaII > sigma_y) {
					
#if (USE_INVETA_EP)
					khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - invEta_EP    );
#else
					khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );
#endif
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



#if (DARCY)
				Py = sigmaII - sigmaT;

#if (1)
				if (phi>=phiCrit) {

					Pe0 = Pe;
					if (Pe < Py && Pe!=0) {

						compute khi_bOld = khi_b;
						khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b );
						Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
						Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc
						if (isnan(Pe)!=0) {
							printf("error in computeEta: Pe  = %.2e\n", Pe);
							exit(0);
						}

					} else {
						khi_b = 1e30;
					}

					// sigma_y chunk
					// =================
					sigma_y = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);
					sigmaT = (cohesion*cos(frictionAngle))/Rad;
					if (Pe<-sigmaT) {
						Py = -sigmaT;
						khi_b = 1.0/((1.0-phi)/Py * (- divV + DeltaP0/(Bulk*dt))   - 1.0/(Bulk*dt) - 1.0/eta_b    );
						Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt));
						Pe = Zb * ( - divV + DeltaP0/(Bulk*dt) ); // Pc


						sigma_y = (sigmaT)/2.0; // arbitrary limit on the minimum mohr circle
					}

					//printf("koko\n");
					// =================
				} else {
					break;
				}
#endif


				if (fabs(1.0-Pe/Py)<yieldTol || Pe>Py-yieldTol) {
					if (sigmaII0<sigma_y) {
						break;
					}
					if (fabs(1.0-sigmaII/sigma_y)<yieldTol) {
						break;
					}
				}
#else
				break;
#endif

			} while (iDum<1) ;

			// Copy updated values back
#if (USE_INVETA_EP)
					Z 	= (1.0-phi)*1.0/(1.0/khi + invEta_EP);
#else
					Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
#endif
			if (Z<Numerics->etaMin) {
				Z = Numerics->etaMin;
			}


			Physics->eta[iCell] = eta;
			Physics->khi[iCell] = khi;
			Physics->G	[iCell] = G;
			Physics->Z	[iCell] = Z;
#if (DARCY)
			Zb 	= (1.0-phi)*1.0/(1.0/khi_b + 1.0/eta_b + 1.0/(Bulk*dt)); // this might not be needed, but it's an extra security

			if (fabs(Zb)<Numerics->etaMin) {
				Zb = Zb/fabs(Zb) * Numerics->etaMin;
			}
			if (fabs(Zb)>Numerics->etaMax) {
				Zb = Zb/fabs(Zb) * Numerics->etaMax;
			}


			Physics->eta_b[iCell] = eta_b;
			Physics->khi_b[iCell] = khi_b;
			Physics->Zb[iCell] = Zb;
#endif

		}
	}


	//Physics_CellVal_SideValues_copyNeighbours_Global(Physics->etaVisc, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);

	Physics_CellVal_SideValues_copyNeighbours_Global(sigma_y_Stored, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(ZprePlasticity, Grid);

#if (DARCY)
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi_b, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta_b, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Zb, Grid);
#endif





	// ================================================================================
	// 									Shear nodes viscosity
	compute sq_sigma_xx0;
	compute sigma_xy0;
	int iNode;
	//#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
#if (DARCY)
			phi = Interp_ECVal_Cell2Node_Local(Physics->phi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
#else
# if (COMPUTE_SHEAR_VISCOSITY)
			phi = 0.0;
			if (ix == 0 || iy == 0 || ix == Grid->nxS-1 || iy == Grid->nyS-1) {
				Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
			} else {
			eta = Interp_ECVal_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			G   = Interp_ECVal_Cell2Node_Local(Physics->G,  ix   , iy, Grid->nxEC);

			sigma_y = Interp_ECVal_Cell2Node_Local(sigma_y_Stored,  ix   , iy, Grid->nxEC);
			sq_sigma_xx0  = Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC];
			sigma_xy0  	  = Physics->sigma_xy_0[iNode];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xy0)*(sigma_xy0)    + 0.25*(sq_sigma_xx0));


			
			
			compute dVxdy, dVydx;
			compute dVxdx, dVydy;
			dVxdy = ( Physics->Vx[(ix  )+(iy+1)*Grid->nxVx]
					- Physics->Vx[(ix  )+(iy  )*Grid->nxVx] )/Grid->dy;


			dVydx = ( Physics->Vy[(ix+1)+(iy  )*Grid->nxVy]
					- Physics->Vy[(ix  )+(iy  )*Grid->nxVy] )/Grid->dx;

			compute Eps_xy = .5*(dVxdy+dVydx);
			
			// Anton's trick
			dVxdx = 0.0;
			dVydy = 0.0;
			compute Exx_x_Sxx0 = 0.0;
			compute Exx_x_Sxx0_ov_G = 0.0;
			compute sq_Exx = 0.0;
			
			int Ix, Iy;
			int IxMod[4] = {0,1,1,0}; // lower left, lower right, upper right, upper left
			int IyMod[4] = {0,0,1,1};
			for (iCell = 0; iCell < 4; ++iCell) {
				Ix = (ix)+IxMod[iCell];
				Iy = (iy)+IyMod[iCell];

				dVxdx = ( Physics->Vx[(Ix  )+(Iy)*Grid->nxVx]
						- Physics->Vx[(Ix-1)+(Iy)*Grid->nxVx] )/Grid->dx;


				dVydy = ( Physics->Vy[(Ix  )+(Iy  )*Grid->nxVy]
						- Physics->Vy[(Ix  )+(Iy-1)*Grid->nxVy] )/Grid->dy;

				
				sq_Exx += 0.25*0.5*(dVxdx-dVydy)*0.5*(dVxdx-dVydy);
#if (USE_SIGMA0_OV_G)
				Exx_x_Sxx0_ov_G += 0.25*(0.5*(dVxdx-dVydy)) * Physics->sigma_xx_0_ov_G[Ix+Iy*Grid->nxEC];
#else 
				Exx_x_Sxx0 += 0.25*(0.5*(dVxdx-dVydy)) * Physics->sigma_xx_0[Ix+Iy*Grid->nxEC];
#endif

			}


			EII = sqrt(sq_Exx + Eps_xy*Eps_xy);


#if (USE_SIGMA0_OV_G)
			printf("USE_SIGMA0_OV_G not implemented for computation on shear nodes\n");
			exit(0);
			
#else
			Eff_strainRate = sqrt(EII*EII + Eps_xy*sigma_xy0/(G*dt) + Exx_x_Sxx0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
#endif
			Z = Interp_ECVal_Cell2Node_Local(ZprePlasticity,  ix   , iy, Grid->nxEC);
			Z = Z*(1.0-phi);

			//printf("Z = %.2e, Zoth = %.2e, eta = %.2e, etaGrid = %.2e, G = %.2e, GGrid = %.2e\n",Z, 1.0/(1.0/eta + 1/(G*dt)), eta, Physics->eta[ix + iy*Grid->nxEC], G, Physics->G[ix + iy*Grid->nxEC]);

			sigmaII = 2.0*Z*Eff_strainRate;

			//printf("Z = %.2e, sigmaII = %.2e, G[0] = %.2e, sigmay = %.2e\n",Z,sigmaII,Physics->G[0],sigma_y);
			khi = 1e30;
			if (sigmaII > sigma_y) {
				//printf("iNode = %i\n", iNode);
				//khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );
				khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );


				//printf("khi = %.2e, khiCorr = %.2e\n",khi, khiCorr);
				if (khi<0.0) {
					// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
					// if it happens then I consider the case where there are == , which means khi -> inf
					printf("WTF!\n");
					khi = 1e30;
					exit(0);
				}


				Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
			}


			if (Z<Numerics->etaMin) {
				Z = Numerics->etaMin;
			}

				Physics->ZShear[iNode] = Z;
				
			
				
			}
			
			Physics->khiShear[iNode] = khi;
#else
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
#endif // end COMPUTE_SHEAR_VISCOSITY
#endif // end DARCY

		}

	}

	// 									Shear nodes viscosity
	// ================================================================================

	free(ZprePlasticity);
	free(sigma_y_Stored);

}




void Physics_Eta_Simple_updateGlobal(Model* Model)
{
	
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);
	BC* BCStokes 			= &(Model->BCStokes);



	int iCell, iy, ix;
	SinglePhase* thisPhaseInfo;
	// ===== get G =====
//#pragma omp parallel for private(iCell, thisPhaseInfo) schedule(dynamic,16)

	//for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			Physics->G[iCell] = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				Physics->G[iCell] += thisPhaseInfo->weight/MatProps->G[thisPhaseInfo->phase];
				thisPhaseInfo = thisPhaseInfo->next;
			}
			Physics->G[iCell] = Physics->sumOfWeightsCells[iCell]	/ Physics->G[iCell] ;
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->G, Grid);
	// ===== get G =====

	// ===== get EffStrainRate =====
	compute* EffStrainRate_CellGlobal = (compute*) malloc(Grid->nECTot*sizeof(compute));
	Physics_Eta_EffStrainRate_getGlobalCell(Model, EffStrainRate_CellGlobal);
	// ===== get EffStrainRate =====


	// ===== get the Z as a visco-elastic predictor =====
	void Physics_Eta_VEpredictor_getGlobalCell(Model, EffStrainRate_CellGlobal);
	// ===== get the Z as a visco-elastic predictor =====

	
	
	
	compute sumOfWeights;
	compute phi, khi, Pe, sigmaII, Z;
	compute EffStrainRate;
	compute sigma_y;
	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;

			
			sumOfWeights 	= Physics->sumOfWeightsCells[iCell];

			phi = 0.0;


			Pe 		= Physics->P [iCell];
			EffStrainRate = EffStrainRate_CellGlobal[iCell];
			Z 	= Physics->Z[iCell];
			sigmaII = 2.0*Z*EffStrainRate;
			khi = 1e30;

			int phase;
			compute weight;
			compute cohesion, frictionAngle;
			cohesion = 0.0;
			frictionAngle = 0.0;
			thisPhaseInfo = Physics->phaseListHead[iCell];
			while (thisPhaseInfo != NULL) {
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				cohesion 		+= MatProps->cohesion[phase] * weight;
				frictionAngle 	+= MatProps->frictionAngle[phase] * weight;
				thisPhaseInfo = thisPhaseInfo->next;
			}
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;



			sigma_y = cohesion * cos(frictionAngle)   +  Pe * sin(frictionAngle);

			// Since there is no griffiths handling for negative pressure for the non darcy case yet
			// here I assume a flat Mohr Coulomb when Pe <0
			if (Pe<=0) {
				sigma_y = cohesion * cos(frictionAngle);
			}


			compute G = Physics->G[iCell];
			compute dt = Physics->dt;
			compute eta = Physics->eta[iCell];


			if (sigmaII > sigma_y) {

				khi = 1.0/((1.0-phi)/sigma_y * (2.0*EffStrainRate)   - 1.0/(G*dt) - 1.0/eta    );

				if (khi<0.0) {
					// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
					// if it happens then I consider the case where there are == , which means khi -> inf
					printf("khi = %.2e, eta = %.2e, G = %.2e, dt = %.2e, Eff_Strainrate = %.2e, 1-phi = %.2e, sigma_y = %.2e, Pe = %.2e, Pmin = %.2e\n", khi, eta, G, dt, EffStrainRate, 1.0-phi, sigma_y, Pe, -cohesion*cos(frictionAngle)/sin(frictionAngle));
					printf("WTF!\n");
					
					exit(0);
				}

				Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));

			} else {
				khi = 1e30;
			}			

			if (Z<Numerics->etaMin) {
				Z = Numerics->etaMin;
			}


			
			Physics->khi[iCell] = khi;
			Physics->Z	[iCell] = Z;

		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->khi, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);
	


	// ================================================================================
	// 									Shear nodes viscosity
	compute sq_sigma_xx0;
	compute sigma_xy0;
	int iNode;
	//#pragma omp parallel for private(iy,ix, iNode) OMP_SCHEDULE
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			Physics->etaShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_ECVal_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}
	// 									Shear nodes viscosity
	// ================================================================================


}





void Physics_Eta_smoothGlobal (Model* Model) 
{
	/*
	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	//Numerics* Numerics 		= &(Model->Numerics);
	//BC* BCStokes 			= &(Model->BCStokes);
	printf("in smooth\n");
	int iy, ix, iN, iCell, iCellN; // where N stands for Neighbours
	//int IxN[4] = {-1,1,0,0};
	//int IyN[4] = {0, 0, -1, 1};

	int IxN[8] = {-1,-1,-1, 1, 1, 1, 0, 0};
	int IyN[8] = {-1, 0, 1,-1, 0, 1,-1, 1};

	compute EtaRatio;
	compute EtaRatio_Tol = 5.0;

	compute GRatio;
	compute GRatio_Tol = 100.0;
	int iCount = 0;

	compute* Zcopy = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Gcopy = (compute*) malloc(Grid->nECTot * sizeof(compute));
	int i;
	for(i=0;i<Grid->nECTot;++i) {
		Zcopy[i] = Physics->Z[i];
		Gcopy[i] = Physics->G[i];
	}
	compute contrib; 
	compute weight;
	
	//int flags = (int*) malloc(Grid->nECTot * sizeof(compute));
	for (iCount = 0; iCount < 100; ++iCount)
	{
		for (iy = 1; iy < Grid->nyEC - 1; iy++)
		{
			for (ix = 1; ix < Grid->nxEC - 1; ix++)
			{
				iCell = ix + iy * Grid->nxEC;

				// Loop over neighbours
				for (iN = 0; iN < 8; ++iN)
				{
					iCellN = ix + IxN[iN] + (iy + IyN[iN]) * Grid->nxEC;
					EtaRatio = Physics->Z[iCell] / Physics->Z[iCellN];
					GRatio = Physics->G[iCell] / Physics->G[iCellN];
					
					contrib = 0.0;
					weight = 1.0;
					if (fabs(EtaRatio) > EtaRatio_Tol)
					{ // Do something only if Eta > EtaN
						Zcopy[iCell] = (Physics->Z[iCell] + Physics->Z[iCellN]) / 2.0;
						//contrib += Physics->Z[iCellN];
						//weight += 1.0;
						//break;
					}
					//Zcopy[iCell] = (Physics->Z[iCell] + contrib)/weight;
					
					if (fabs(GRatio) > GRatio_Tol)
					{ // Do something only if Eta > EtaN
						//Gcopy[iCell] = 2.0/(1.0/Physics->G[iCell] + 1.0/Physics->G[iCellN]);
						Gcopy[iCell] = pow(10.0,(log10(Physics->G[iCell]) + log10(Physics->G[iCellN]))/2.0);
						
						
						//Gcopy[iCell] = (Physics->G[iCell] + Physics->G[iCellN]) / 2.0;
						//printf("iCell = %i, Gcopy = %.2e, Gl = %.2e, Gr = %.2e\n", iCell, Gcopy[iCell], Physics->G[iCell], Physics->G[iCellN]);
					}
					
				}
			}
		}
		// copy the changes to Z
		for(i=0;i<Grid->nECTot;++i) {
			Physics->Z[i] = Zcopy[i];
			//Physics->G[i] = Gcopy[i];
			//Physics->Z[i] = 1.0/(1.0/Physics->eta[i] + 1.0/(Physics->G[i]*Physics->dt) + 1.0/Physics->khi[i]);
		}
	}
	

	free(Zcopy);
	free(Gcopy);
	*/

}


void Physics_Eta_FromParticles_updateGlobal(Model* Model)
{
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	BC* BCThermal 			= &(Model->BCThermal);
	Numbering* NumThermal 	= &(Model->NumThermal);
	Numerics* Numerics 		= &(Model->Numerics);

	int signX, signY;
	compute locX, locY;
	int ix, iy;

	INIT_PARTICLE

	compute Dsigma_xx_0_Grid;
	compute Dsigma_xy_0_Grid;

	compute sigma_xx_0_Grid;
	compute sigma_xy_0_Grid;

	
	compute EII;
	compute* EIICell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* SII0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exx = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* dVxdyGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* dVydxGrid = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* Rotxy = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute dVxdy, dVydx, dVxdx, dVydy;
	int iCell;
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

			//dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx] - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			//dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy] - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;


			//Exx[iCell]  = 0.5*(dVxdx-dVydy);
		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx, Grid);

	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			//dVxdy = (Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx] - Physics->Vx[(ix  ) + (iy  )*Grid->nxVx])/Grid->dy;
			//dVydx = (Physics->Vy[(ix+1) + (iy  )*Grid->nxVy] - Physics->Vy[(ix  ) + (iy  )*Grid->nxVy])/Grid->dx;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Exy[iNode] = 0.5*(dVxdy+dVydx);
			dVxdyGrid[iNode] =  dVxdy;
			dVydxGrid[iNode] =  dVydx;

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

				/*
				// ===== weight cells =====
				int signX, signY;
				int i;
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
				iCell = (ix+IxN[i] + (iy+IyN[i]) * nxEC);
				weightCell = fabs(locX)*fabs(locY);
				// ===== weight cells =====


				// ===== weight nodes =====

				locX = fabs(locX);
				locY = fabs(locY);


				weight = (1.0 - locX) * (1.0 - locY);
				iNodeNeigh = iNode;
				// ===== weight nodes =====
				*/
				//Physics->sigma_xy_0_ov_G 		[iNodeNeigh] += (thisParticle->sigma_xy_0 / MatProps->G[phase]) * weight;












				//Physics->eta				[iCell] += eta * weight;




















				compute EII;
				ExxPart = Interp_ECVal_Cell2Particle_Local(Exx, ix, iy, Grid->nxEC, locX, locY);
				sigma_xx_0_Grid = Interp_ECVal_Cell2Particle_Local(Physics->sigma_xx_0, ix, iy, Grid->nxEC, locX, locY);
				ExyPart = Interp_NodeVal_Node2Particle_Local(Exy, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				sigma_xy_0_Grid = Interp_NodeVal_Node2Particle_Local(Physics->sigma_xy_0, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				//EII = Interp_ECVal_Cell2Particle_Local(EIICell, ix, iy, Grid->nxEC, locX, locY);
				
				EII = sqrt(ExxPart*ExxPart + ExyPart*ExyPart);

				compute eta;
				int phase = thisParticle->phase;
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
				compute SII0 = sqrt(Sxx0*Sxx0 + Sxy0*Sxy0);

				

				compute RotxyPart = Interp_NodeVal_Node2Particle_Local(Rotxy, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
#if (USE_UPPER_CONVECTED) 
				int nxEC = Grid->nxEC;
				int nxS = Grid->nxS;
				int nyS = Grid->nyS;

				compute sqEII_Part = ExxPart*ExxPart + ExyPart*ExyPart;
				compute sqSII0_Part = Sxx0*Sxx0 + Sxy0*Sxy0;
				


				compute dVxdyPart = Interp_NodeVal_Node2Particle_Local(dVxdyGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				compute dVydxPart = Interp_NodeVal_Node2Particle_Local(dVydxGrid, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
				compute Eff_strainRate = 1.0/(2.0*G*dt) * sqrt(pow((2.0*ExxPart*G*dt + Sxx0 + 2.0*dt*(Sxx0*ExxPart + Sxy0*dVxdyPart)),2.0) + pow((2.0*ExyPart*G*dt - Sxx0*dt*2.0*RotxyPart+ Sxy0),2.0));

#else
				compute Eff_strainRate = sqrt(EII*EII + ExxPart*Sxx0/(G*dt) + ExyPart*Sxy0/(G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*SII0*SII0   );
#endif


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

				
				

				thisParticle = thisParticle->next;
			}
		}
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



void Physics_Eta_EffStrainRate_getGlobalCell(Model* Model, compute* EffStrainRate) {
	Grid* Grid 				= &(Model->Grid);
	MatProps* MatProps 		= &(Model->MatProps);
	Particles* Particles 	= &(Model->Particles);
	Physics* Physics 		= &(Model->Physics);
	BC* BCStokes 			= &(Model->BCStokes);
	BC* BCThermal 			= &(Model->BCThermal);
	Numbering* NumThermal 	= &(Model->NumThermal);
	Numerics* Numerics 		= &(Model->Numerics);

	int signX, signY;
	compute locX, locY;
	int ix, iy;


	compute Dsigma_xx_0_Grid;
	compute Dsigma_xy_0_Grid;

	compute sigma_xx_0_Grid;
	compute sigma_xy_0_Grid;

	int Mode = 1; // 0: stress based, 1: strain rate based
	
	compute EII;
	compute* EII_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* SII0_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exx_CellGlobal = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Exy_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* dVxdy_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute* dVydx_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));

	compute* Rotxy_NodeGlobal = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute dVxdy, dVydx, dVxdx, dVydy;
	int iCell;

	compute sq_sigma_xy0, sigma_xx0;

	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			
			dVxdx = (Physics->Vx[(ix) + (iy)*Grid->nxVx]
						 - Physics->Vx[(ix-1) + (iy)*Grid->nxVx])/Grid->dx;

			dVydy = (Physics->Vy[(ix) + (iy)*Grid->nxVy]
						 - Physics->Vy[(ix) + (iy-1)*Grid->nxVy])/Grid->dy;

			Exx_CellGlobal[iCell] = 0.5*(dVxdx-dVydy);

			Physics_StrainRateInvariant_getLocalCell(Model, ix, iy, &EII);
			

			sq_sigma_xy0  = Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy-1)*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix-1+(iy  )*Grid->nxS];
			sq_sigma_xy0 += Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS] * Physics->sigma_xy_0[ix  +(iy  )*Grid->nxS];
			sigma_xx0  = Physics->sigma_xx_0[iCell];// + Physics->Dsigma_xx_0[iCell];
			SII0_CellGlobal[iCell] = sqrt((sigma_xx0)*(sigma_xx0)    + 0.25*(sq_sigma_xy0));

		}
	}
	Physics_CellVal_SideValues_copyNeighbours_Global(Exx_CellGlobal, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(SII0_CellGlobal, Grid);

	int iNode;
	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			//dVxdy = (Physics->Vx[(ix  ) + (iy+1)*Grid->nxVx] - Physics->Vx[(ix  ) + (iy  )*Grid->nxVx])/Grid->dy;
			//dVydx = (Physics->Vy[(ix+1) + (iy  )*Grid->nxVy] - Physics->Vy[(ix  ) + (iy  )*Grid->nxVy])/Grid->dx;

			dVxdy = ( Physics->Vx[ix  + (iy+1)*Grid->nxVx]  - Physics->Vx[ix  + (iy  )*Grid->nxVx] )/Grid->dy;

			dVydx = ( Physics->Vy[ix+1+ iy*Grid->nxVy]	  - Physics->Vy[ix  + iy*Grid->nxVy] )/Grid->dx;
			Exy_NodeGlobal[iNode] = 0.5*(dVxdy+dVydx);
			dVxdy_NodeGlobal[iNode] =  dVxdy;
			dVydx_NodeGlobal[iNode] =  dVydx;

			Rotxy_NodeGlobal[iNode] = 0.5*(dVxdy-dVydx);

		}
	}


	for (iy = 1; iy<Grid->nyEC-1; iy++) {
		for (ix = 1; ix<Grid->nxEC-1; ix++) {
			iCell = ix + iy*Grid->nxEC;
			compute Exx = Exx_CellGlobal[iCell];
			compute dVxdx = Exx_CellGlobal[iCell];
			
			compute Txx0 = Physics->sigma_xx_0[iCell];
			compute Txy0 = Interp_NodeVal_Node2Cell_Local(Physics->sigma_xy_0,ix,iy,Grid->nxS);
			dVxdy = Interp_NodeVal_Node2Cell_Local(dVxdy_NodeGlobal,ix,iy,Grid->nxS);
			dVydx = Interp_NodeVal_Node2Cell_Local(dVydx_NodeGlobal,ix,iy,Grid->nxS);
			compute Exy =Interp_NodeVal_Node2Cell_Local(Exy_NodeGlobal,ix,iy,Grid->nxS);
			compute G = Physics->G[iCell];
			compute dt = Physics->dt;
			EffStrainRate[iCell] = 1.0/(2.0*G*dt) * sqrt(pow((2.0*Exx*G*dt + Txx0 + 2.0*dt*(Txx0*dVxdx + Txy0*dVxdy)),2.0) + pow((2.0*Exy*G*dt - Txx0*dt*(dVxdy - dVydx) + Txy0),2.0));
		}
	}

	Physics_CellVal_SideValues_copyNeighbours_Global(EffStrainRate, Grid);

}


void Physics_Eta_VEpredictor_getGlobalCell(Model* Model, compute* EffStrainRate) {
	// Update Physics->Z and Physics->eta
	// according to the Visco-elastic predictor

	Grid* Grid 				= &(Model->Grid);
	Physics* Physics 		= &(Model->Physics);
	MatProps* MatProps 		= &(Model->MatProps);
	Numerics* Numerics 		= &(Model->Numerics);
	BC* BCStokes 			= &(Model->BCStokes);

	int iy, ix, iCell;

	compute T, P, phi;
	compute alpha;

	compute eta, G;
	SinglePhase* thisPhaseInfo;

	compute maxInvVisc = 0.0;
	compute invEtaDiff, invEtaDisl, invEtaPei;
	int phase;
	compute weight;
	compute B, E, V, n, taup, gamma, q, s;
	compute R = Physics->R;
	compute ZUpper, ZLower;
	compute BDiff[NB_PHASE_MAX], BDisl[NB_PHASE_MAX], BPei[NB_PHASE_MAX];
	compute EII;

	compute eta_thisPhase;
	compute sumOfWeights 	= Physics->sumOfWeightsCells[iCell];

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
				thisPhaseInfo 	= thisPhaseInfo->next;
				eta_thisPhase = (1.0 / (invEtaDiff + invEtaDisl + invEtaPei));


				eta += weight * eta_thisPhase;

				invEta_EP += (1.0/(MatProps->G[phase]*dt)+1.0/eta_thisPhase) * weight;


			}
			G 				 = sumOfWeights	/ G;
			
			eta 			/= sumOfWeights;

			invEta_EP /= sumOfWeights;

			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}
#if (USE_INVETA_EP)
			ZLower = 1.0/invEta_EP;
#else
			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);
#endif

			Z = 0.5*((1.0-phi)*ZUpper+(1.0-phi)*ZLower);
			Zcorr = Z;



			sigmaII = 2.0*Z*EffStrainRate[iCell];

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
				//eta = pow(10.0,eta / sumOfWeights);
				//invEta_EP = pow(10.0,invEta_EP / sumOfWeights);
				invEta_EP = invEta_EP / sumOfWeights;
				PrevZcorr = Zcorr;
#if (USE_INVETA_EP)
				Zcorr = (1.0-phi)*(1.0/(invEta_EP)) - Z;
#else
				Zcorr = (1.0-phi)*(1.0/(1.0/(G*dt) + 1.0/eta)) - Z;
#endif

				
				
				if (Zcorr/PrevZcorr<-0.9) {
					alpha = alpha/2.0;
				}
				Z += alpha*Zcorr;

				sigmaII = 2.0*Z*EffStrainRate[iCell];
			}


			Physics->Z[iCell] = Z;
			Physics->eta[iCell] = eta;

		}
	}

	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->eta, Grid);
	Physics_CellVal_SideValues_copyNeighbours_Global(Physics->Z, Grid);


	
}