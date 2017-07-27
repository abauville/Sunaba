/*
 * Physics_Eta.c
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */


#include "stokes.h"

void Physics_Eta_init(Physics* Physics, Grid* Grid, MatProps* MatProps, Numerics* Numerics) {

	int iy, ix, iCell;
	SinglePhase* thisPhaseInfo;
	compute P, T;
	int phase;
	compute EII, weight;
	compute B, E, V, Binc, n, taup, q, s, gamma;
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
			while (thisPhaseInfo != NULL) {
				invEtaDiff = 0.0;
				invEtaDisl = 0.0;
				invEtaPei  = 0.0;
				phase = thisPhaseInfo->phase;
				weight = thisPhaseInfo->weight;
				G 				+= weight/MatProps->G[phase];
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

				eta += weight * eta_thisPhase;





			}
			eta = eta / sumOfWeights;
			if (eta>Numerics->etaMax) {
				eta = Numerics->etaMax;
			}
			if (eta<Numerics->etaMin) {
				eta = Numerics->etaMin;
			}

			Physics->eta[iCell] = eta;

			Physics->G[iCell]  = Physics->sumOfWeightsCells[iCell]/G;
			Physics->khi[iCell] = 1E30;

			Physics->Z[iCell] = 1.0/( 1.0/Physics->khi[iCell] + 1.0/Physics->eta[iCell] + 1.0/(Physics->G[iCell]*Physics->dt) );

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
			Physics->etaShear[iNode] = Interp_Any_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_Any_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_Any_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}

}




void Physics_Eta_updateGlobal(Physics* Physics, Grid* Grid, Numerics* Numerics, BC* BCStokes,MatProps* MatProps)
{
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
	compute strainReductionFac = 0.9; // 1.0 stays the same
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

			Physics_StrainRateInvariant_getLocalCell(Physics, Grid, ix, iy, &EII);
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

				Exy_x_Sxy0 += (0.5*(dVxdy+dVydx)) * Physics->sigma_xy_0[Ix+Iy*Grid->nxS];
			}
			Exy_x_Sxy0 /= 4.0; // Eps_xy*sigma_xy0














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





			}
			G 				 = sumOfWeights	/ G;
			eta 			/= sumOfWeights;
			cohesion 		/= sumOfWeights;
			frictionAngle 	/= sumOfWeights;


#if (STRAIN_SOFTENING)
			compute strainLimit = 1.0;
			compute coeff = (1.0-Physics->strain[iCell]/strainLimit);
			if (coeff<strainReductionFac){
				coeff = strainReductionFac;
			}
			frictionAngle *= coeff;
#endif



			maxInvVisc = fmax(1.0/(G*dt),maxInvVisc);
			ZUpper = 1.0/maxInvVisc;
			if (ZUpper>1e10) {
				ZUpper = 1e10;
			}
			ZLower = 1.0/(1.0/(G*dt) + 1.0/eta);

			Z = 0.5*((1.0-phi)*ZUpper+(1.0-phi)*ZLower);
			Zcorr = Z;

			Eff_strainRate = sqrt(EII*EII + Eps_xx*sigma_xx0/(2.0*G*dt) + Exy_x_Sxy0/(2.0*G*dt) + (1.0/(2.0*G*dt))*(1.0/(2.0*G*dt))*sigmaII0*sigmaII0   );
			sigmaII = 2.0*Z*Eff_strainRate;

			// compute viscosities using sigmaII
			while (fabs(Zcorr/Z)>tol) {
				eta = 0.0;
				thisPhaseInfo = Physics->phaseListHead[iCell];

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



				}

				eta 			/= sumOfWeights;

				PrevZcorr = Zcorr;
				Zcorr = (1.0-phi)*(1.0/(1.0/(G*dt) + 1.0/eta)) - Z;
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
			Pe 		= Physics->P [iCell] - Pf;

#endif

			Z 	= (1.0-phi)*1.0/(1.0/eta + 1.0/(G*dt));


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
			if (Pe<0) {
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
					khi = 1.0/((1.0-phi)/sigma_y * (2.0*Eff_strainRate)   - 1.0/(G*dt) - 1.0/eta    );
					if (khi<0.0) {
						// quite rare case where (1.0-phi)/sigma_y * (2.0*Eff_strainRate) <  - 1.0/(G*dt) - 1.0/eta
						// if it happens then I consider the case where there are == , which means khi -> inf
						printf("khi = %.2e, eta = %.2e, G = %.2e, dt = %.2e, Eff_Strainrate = %.2e, 1-phi = %.2e, sigma_y = %.2e, Pe = %.2e, Pmin = %.2e\n", khi, eta, G, dt, Eff_strainRate, 1.0-phi, sigma_y, Pe, -cohesion*cos(frictionAngle)/sin(frictionAngle));
						printf("WTF!\n");
						khi = 1e30;
						//exit(0);
					}
					Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt));
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
			Z 	= (1.0-phi)*1.0/(1.0/khi + 1.0/eta + 1.0/(G*dt)); // this might not be needed, but it's an extra security
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
			Physics->etaShear[iNode] = Interp_Any_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			Physics->khiShear[iNode] = Interp_Any_Cell2Node_Local(Physics->khi,  ix   , iy, Grid->nxEC);
#if (DARCY)
			phi = Interp_Any_Cell2Node_Local(Physics->phi,  ix   , iy, Grid->nxEC);
			Physics->ZShear[iNode] = Interp_Any_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
#else
			/*
			phi = 0.0;

			eta = Interp_Any_Cell2Node_Local(Physics->eta,  ix   , iy, Grid->nxEC);
			G = Interp_Any_Cell2Node_Local(Physics->G,  ix   , iy, Grid->nxEC);

			sigma_y = Interp_Any_Cell2Node_Local(sigma_y_Stored,  ix   , iy, Grid->nxEC);
			sq_sigma_xx0  = Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy+1)*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix+1+(iy  )*Grid->nxEC];
			sq_sigma_xx0 += Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC] * Physics->sigma_xx_0[ix  +(iy  )*Grid->nxEC];
			sigma_xy0  	  = Physics->sigma_xy_0[iNode];// + Physics->Dsigma_xx_0[iCell];
			sigmaII0 = sqrt((sigma_xy0)*(sigma_xy0)    + 0.25*(sq_sigma_xx0));


			Physics_StrainRateInvariant_getLocalNode(Physics,BCStokes,Grid,ix,iy,&EII);


			Eff_strainRate = EII + (1.0/(2.0*G*dt))*sigmaII0;

			Z = Interp_Any_Cell2Node_Local(ZprePlasticity,  ix   , iy, Grid->nxEC);
			Z = Z*(1.0-phi);

			//printf("Z = %.2e, Zoth = %.2e, eta = %.2e, etaGrid = %.2e, G = %.2e, GGrid = %.2e\n",Z, 1.0/(1.0/eta + 1/(G*dt)), eta, Physics->eta[ix + iy*Grid->nxEC], G, Physics->G[ix + iy*Grid->nxEC]);

			sigmaII = 2.0*Z*Eff_strainRate;

			//printf("Z = %.2e, sigmaII = %.2e, G[0] = %.2e, sigmay = %.2e\n",Z,sigmaII,Physics->G[0],sigma_y);

			if (sigmaII > sigma_y) {
				//printf("iCell = %i, C = %i\n", iCell, C);
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
			if (ix == 0 || iy == 0 || ix == Grid->nxS-1 || iy == Grid->nyS-1) {
				Physics->ZShear[iNode] = Interp_Any_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
			}
			*/

			Physics->ZShear[iNode] = Interp_Any_Cell2Node_Local(Physics->Z,  ix   , iy, Grid->nxEC);
#endif

		}

	}

	// 									Shear nodes viscosity
	// ================================================================================

	free(ZprePlasticity);
	free(sigma_y_Stored);

}



