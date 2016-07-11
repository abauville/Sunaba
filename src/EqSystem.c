/*
 * eqSystem.c
 *
 *  Created on: Feb 25, 2016
 *      Author: abauville
 */


#include "stokes.h"

static void Static_LocalStencilVx(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift);
static void Static_LocalStencilVy(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift);
static void Static_LocalStencilP(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift);

#if (HEAT)
static void Static_LocalStencilT(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift);
#endif


void EqSystem_allocateI (EqSystem* EqSystem)
{
	EqSystem->I 			= (int*) malloc((EqSystem->nEq+1)  * sizeof(int));
}

void EqSystem_allocateMemory(EqSystem* EqSystem)
{
	EqSystem->J = (int*)     malloc(EqSystem->nnz * sizeof(int));
	EqSystem->V = (compute*) malloc(EqSystem->nnz * sizeof(compute));
	EqSystem->b = (compute*) malloc( EqSystem->nEq * sizeof(compute));
	EqSystem->x = (compute*) malloc( EqSystem->nEq * sizeof(compute));

}

void EqSystem_freeMemory(EqSystem* EqSystem, Solver* Solver)
{
	//Free Pardiso

	int error = 0;
	int idum = 0;
	double ddum = 0;
	int phase = -1;                 // Release internal memory.
	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, &ddum, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);

	free(EqSystem->I);
	free(EqSystem->J);
	free(EqSystem->V);
	free(EqSystem->b);
	free(EqSystem->x);

}

void EqSystem_assemble(EqSystem* EqSystem, Grid* Grid, BC* BC, Physics* Physics, Numbering* Numbering)
{

	//==========================================================================
	//
	//                          BUILD SPARSE TRIPLET
	//
	//==========================================================================


	// Init variables
	// ===============================
	int nLoc;
	int i, J;



	int SetupType = BC->SetupType;





	int Jloc[11];
	compute Vloc[11];
	compute bloc;

	int shift = 0;




	// Reinitialize b
	/*
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->J[i] = -1;
	}

	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->b[i] = 0.0;
	}
*/




	// Fill J, V and EqSystem->b for Free and Dirichlet nodes
	// ===============================================
	if (DEBUG) {
		printf("Start Filling loop\n");
		printf("nEq = %i, EqSystem->nRow = %i\n", EqSystem->nEq, EqSystem->nRow);
	}
	int iEq, ix, iy, I;
	int IC = 0;
	int Iloc, IBC;
	StencilType Stencil;
	//int INumMap;

	int order[11] = {0,1,2,3,4,5,6,7,8,9,10};

#pragma omp parallel for private(iEq, I, ix, iy, i, Stencil, order, nLoc, IC, Jloc, Vloc, bloc, shift, J,  Iloc, IBC) schedule(static,32)
	for (iEq=0; iEq<EqSystem->nEq; iEq++) {

		I = EqSystem->I[iEq];
		ix = Numbering->IX[iEq];
		iy = Numbering->IY[iEq];

		i = 1;
		while (iEq>=Numbering->subEqSystem0[i]) {
			i++;
		}
		Stencil = Numbering->Stencil[i-1];


		// Reinitialize order
		for (i=0;i<11;i++) {
			order[i] = i;
		}




		if (Stencil==Vx)		{
			nLoc = 11;
			IC = 2;
			Static_LocalStencilVx(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift);
		}
		else if (Stencil==Vy) 	{
			nLoc = 11;
			IC = 6;
			Static_LocalStencilVy(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift);
		}
		else if (Stencil==P) 	{
			nLoc = 4;
			Static_LocalStencilP(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift);
		}
#if (HEAT)
		else if (Stencil==T) 	{
			nLoc = 5;
			IC = 2;
			Static_LocalStencilT(order, Jloc, Vloc, &bloc, ix, iy, Grid, Physics, SetupType, &shift);
		}
#endif



		// ===========================================
		// Fill the right hand side and apply BC
		// ===========================================

		// Fill right hand side with the local right hand side
		EqSystem->b[iEq] = bloc;

		for (i=0; i<nLoc; i++) {
			//printf ("%i ",Jloc[order[i]]);
			Iloc = Numbering->map[Jloc[order[i]]];

			if (Iloc < 0) { // if Boundary node
				IBC = abs(Iloc) - 1;
				if (BC->type[IBC]==Dirichlet) { // Dirichlet on normal node
					//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC];
					EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC];
				}
				else  if (BC->type[IBC]==DirichletGhost) { // Dirichlet
					Vloc[order[IC]]  += -Vloc[order[i]]; // +1 to VxC
					//EqSystem->b[iEq] += -Vloc[i] * 2*BC->value[IBC];
					EqSystem->b[iEq] += -Vloc[order[i]] * 2*BC->value[IBC];

				}
				else if (BC->type[IBC]==NeumannGhost) { // NeumannGhost
					Vloc[order[IC]] += Vloc[order[i]]; // +1 to VxC
					if (Stencil==Vx) {
						if 		(i==0) { // VxS
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[0];
						}
						else if (i==4) { // VxN
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dy;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DYEC[Grid->nyS-1];
						}
					}


					else if (Stencil==Vy) {
						if 		(i==5) { // VyW
							//EqSystem->b[iEq] += -Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += -Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[0];
						}
						else if (i==7) { // VyE
							//EqSystem->b[iEq] += +Vloc[i] * BC->value[IBC] * dx;
							EqSystem->b[iEq] += +Vloc[order[i]] * BC->value[IBC] * Grid->DXEC[Grid->nxS-1];
						}
					}

					else if (Stencil==Vy) {
						// For the moment only 0 gradient is implement
						// This section should be filled to account for a given gradient
					}

				}
				else {
					printf("error: unknown boundary type\n");
					exit(0);
				}
			}
		}




		// ===========================================
		// Fill J and V global
		// ===========================================

		// Jloc gets the numbering with Dirichlet;
		// When values are inputted in Jsparse they are transformed to NumMap[Jloc[i]]
		// i.e. numbering without dirichlet
		J = 0;
		for (i=0; i<nLoc; i++) {
			if (Numbering->map[Jloc[i]] >=0) { // if free
				if (i>=shift) {
					EqSystem->J[I+J] = Numbering->map[Jloc[i]];
					EqSystem->V[I+J] = Vloc[i];
					J++;
				}
			}

		}



	} // end of the equation loop




	// Explicitly add zeros in the diagonal for the pressure equations (required for compatibility with Pardiso, i.e. to make the matrix square)
	if (UPPER_TRI) {
		for (i=EqSystem->nRow; i<EqSystem->nEq; i++) {
			EqSystem->J[EqSystem->I[i]] = i;
			EqSystem->V[EqSystem->I[i]] = 0.0;
		}
	}



}






static void Static_LocalStencilVx(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift)
{


	int VxPeriod = 0;
	int PPeriod  = 0;
	int	NormalPeriod = 0;



	int NormalE, NormalW, ShearN, ShearS;
	// Define size variables
	// ===============================
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y

	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;

	int nxN = Grid->nxC;
	int nxEC = Grid->nxEC;

	int nxS = Grid->nxS;

	compute EtaN, EtaS, EtaE, EtaW;
	compute ZN, ZS, ZW, ZE; // visco-elasticity factor


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);

	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXS[0]+Grid->DXS[nxVx-1]);
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVx-1) {
			dxW = Grid->DXS[ix-1];
			dxE = 0.5*(Grid->DXS[0]+Grid->DXS[nxVx-1]);
			dxC = 0.5*(dxW+dxE);

		} else {
			dxW = Grid->DXS[ix-1];
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);
		}

	} else {
		dxW = Grid->DXS[ix-1];
		dxE = Grid->DXS[ix];
		dxC = 0.5*(dxW+dxE);

	}

	compute dt = Physics->dt;
	compute sigma_xx_0_E, sigma_xx_0_W, sigma_xy_0_N, sigma_xy_0_S;

	compute GShearN, GShearS;


	//printf("dxW = %.2f, dyS = %.2f, Grid->dx = %.2f, ix = %i, iy = %i,Grid->nxVx = %i, Grid->nxS = %i\n",dxW, dyS,Grid->dx,ix,iy,Grid->nxVx, Grid->nxS);

	if (UPPER_TRI) {
		*shift = 2;
	}
	else {
		*shift = 0;
	}

	// Special case for periodic BC
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 1;
			}
			VxPeriod = Grid->nxVx-1;
			PPeriod  = nxN;

			NormalPeriod = nxN;

			order[ 0] =  0; // VxS
			order[ 1] =  3; // VxW
			order[ 2] =  1; // VxC
			order[ 3] =  2; // VxE
			order[ 4] =  4; // VxN
			order[ 5] =  5; // VySW
			order[ 6] =  6; // VySE
			order[ 7] =  7; // VyNW
			order[ 8] =  8; // VyNE
			order[ 9] = 10; // PW
			order[10] =  9; // PE
		}
		if (ix==Grid->nxVx-2) {
			if (UPPER_TRI) {
				*shift = 3;
			}
			order[ 0] =  0; // VxS
			order[ 1] =  2; // VxW
			order[ 2] =  3; // VxC
			order[ 3] =  1; // VxE
			order[ 4] =  4; // VxN
			order[ 5] =  6; // VySW
			order[ 6] =  5; // VySE
			order[ 7] =  8; // VyNW
			order[ 8] =  7; // VyNE
			order[ 9] =  9; // PW
			order[10] = 10; // PE
		}
	}

	// =====================================================================
	//                           Fill Jloc, Vloc, bloc
	// =====================================================================

	Jloc[order[ 0]]  =   ix      + iy*nxVx     - nxVx          ; // VxS
	Jloc[order[ 1]]  =   ix      + iy*nxVx     - 1      + VxPeriod; // VxW
	Jloc[order[ 2]]  =   ix      + iy*nxVx                     ; // VxC
	Jloc[order[ 3]]  =   ix      + iy*nxVx     + 1             ; // VxE
	Jloc[order[ 4]]  =   ix      + iy*nxVx     + nxVx          ; // VxN
	Jloc[order[ 5]] =   ix+0    + (iy-1)*nxVy + nVxTot        ; // VySW
	Jloc[order[ 6]] =   ix+1    + (iy-1)*nxVy + nVxTot        ; // VySE
	Jloc[order[ 7]] =   ix+0    + iy*nxVy     + nVxTot        ; // VyNW
	Jloc[order[ 8]] =   ix+1    + iy*nxVy     + nVxTot        ; // VyNE
	Jloc[order[ 9]]   =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot + PPeriod; // PW
	Jloc[order[10]]   =   ix      + (iy-1)*nxN  + nVxTot+nVyTot ; // PE

	NormalE = ix  +1    + (iy-1+1)*nxEC;
	NormalW = ix-1+1    + (iy-1+1)*nxEC + NormalPeriod;
	ShearN = ix      + iy*nxS;
	ShearS = ix      + (iy-1)*nxS;

	EtaN    = shearValue(Physics->eta, ix,  iy   , nxEC); // Shear N
	EtaS    = shearValue(Physics->eta, ix, (iy-1), nxEC); // ShearS
	EtaE    = Physics->eta[ NormalE ]; // NormalE
	EtaW    = Physics->eta[ NormalW ]; // NormalW

	GShearN = shearValue(Physics->G, ix,  iy   , nxEC);
	GShearS = shearValue(Physics->G, ix, (iy-1), nxEC);

	ZN = (dt*GShearN ) / (dt*GShearN  + EtaN);
	ZS = (dt*GShearS) / (dt*GShearS + EtaS);
	ZE = (dt*Physics->G     [NormalE]) / (dt*Physics->G     [NormalE] + EtaE);
	ZW = (dt*Physics->G     [NormalW]) / (dt*Physics->G     [NormalW] + EtaW);

	sigma_xx_0_E =  Physics->sigma_xx_0[NormalE];
	sigma_xx_0_W =  Physics->sigma_xx_0[NormalW];
	sigma_xy_0_N =  Physics->sigma_xy_0[ShearN ];
	sigma_xy_0_S =  Physics->sigma_xy_0[ShearS ];

	// Fill Vloc: list of coefficients
	// ================================
	Vloc[order[ 0]] =  EtaS*ZS/dyS/dyC;
	Vloc[order[ 1]] =  2.0 * EtaW*ZW/dxW/dxC;
	Vloc[order[ 2]] = -2.0 * EtaE*ZE/dxE/dxC   -2.0 * EtaW*ZW/dxW/dxC   -1.0 * EtaN*ZN/dyN/dyC   -1.0 * EtaS*ZS/dyS/dyC;
	Vloc[order[ 3]] =  2.0 * EtaE*ZE/dxE/dxC;
	Vloc[order[ 4]] =  EtaN*ZN/dyN/dyC;
	Vloc[order[ 5]] =  EtaS*ZS/dxW/dyS;
	Vloc[order[ 6]] = -EtaS*ZS/dxE/dyS;
	Vloc[order[ 7]] = -EtaN*ZN/dxW/dyN;
	Vloc[order[ 8]] =  EtaN*ZN/dxW/dyN;
	Vloc[order[ 9]] =  1.0/dxW;
	Vloc[order[10]] = -1.0/dxE;

	*bloc = - Physics->g[0] * 0.5 * ( Physics->rho[NormalE] + Physics->rho[NormalW] );

	// add contributions of old stresses
	*bloc += - ( sigma_xx_0_E*(1-ZE)  -   sigma_xx_0_W*(1-ZW))/dxC  -  (sigma_xy_0_N*(1-ZN)  -  sigma_xy_0_S*(1-ZS))/dyC;
}









static void Static_LocalStencilVy(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift)
{


	int VxPeriod = 0;
	int VyPeriod = 0;
	int PPeriod  = 0;
	int ShearPeriod = 0;



	int NormalN, NormalS, ShearE, ShearW;
	// Define size variables
	// ===============================
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y

	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;

	int nxN = Grid->nxC;
	int nxEC = Grid->nxEC;

	int nxS = Grid->nxS;

	compute EtaN, EtaS, EtaE, EtaW;
	compute ZN, ZS, ZE, ZW;
	compute sigma_yy_0_N, sigma_yy_0_S, sigma_xy_0_E, sigma_xy_0_W;
	compute GShearE, GShearW;
	compute dxW = Grid->DXEC[ix-1];
	compute dxE = Grid->DXEC[ix  ];
	compute dxC = 0.5*(dxW+dxE);
	compute dyS = Grid->DYS [iy-1];
	compute dyN = Grid->DYS [iy  ];
	compute dyC = 0.5*(dyS+dyN);
	compute dt = Physics->dt;




	if (UPPER_TRI) {
		*shift = 6;
	}
	else {
		*shift = 0;
	}
	// =====================================================================
	//                          Local numbering
	// =====================================================================
	NormalN = ix-1+1    + (iy  +1)*nxEC ;
	NormalS = ix-1+1    + (iy-1+1)*nxEC ;
	ShearE  = ix      + iy*nxS    ;
	ShearW  = ix-1    + iy*nxS    ;

	VxPeriod = 0		; // VxSW
	VyPeriod = 0 		; // VyW
	PPeriod  = 0   		; // PS
	ShearPeriod = 0;

	// Special cases for periodic BC
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 5;
			}
			VxPeriod = nxVx-1		; // VxSW
			VyPeriod  = nxVy-2  	; // VyW
			PPeriod   = nxN    		; // PS
			ShearPeriod = nxS-1;


			NormalN += nxN			;
			NormalS += nxN			;
			ShearW  += nxS-1 		;

			order[ 0] =  1; // VxSW
			order[ 1] =  0; // VxSE
			order[ 2] =  3; // VxNW
			order[ 3] =  2; // VxNE
			order[ 4] =  4; // VyS
			order[ 5] =  7; // VyW
			order[ 6] =  5; // VyC
			order[ 7] =  6; // VyE
			order[ 8] =  8; // VyN
			order[ 9] =  9; // PS
			order[10] = 10; // PN
		}
		else if (ix==nxVy-3) {
			if (UPPER_TRI) {
				*shift = 7;
			}
			order[ 0] =  0; // VxSW
			order[ 1] =  1; // VxSE
			order[ 2] =  2; // VxNW
			order[ 3] =  3; // VxNE
			order[ 4] =  4; // VyS
			order[ 5] =  6; // VyW
			order[ 6] =  7; // VyC
			order[ 7] =  5; // VyE
			order[ 8] =  8; // VyN
			order[ 9] =  9; // PS
			order[10] = 10; // PN
		}
	}


	// =====================================================================
	//                        fill Jloc, Vloc, bloc
	// =====================================================================
	// Get Viscosities
	// ================
	EtaN    = Physics->eta[NormalN];
	EtaS    = Physics->eta[NormalS];
	EtaE    = shearValue(Physics->eta,  ix   , iy, nxEC);
	EtaW    = shearValue(Physics->eta, (ix-1)+ShearPeriod, iy, nxEC);

	GShearE = shearValue(Physics->G,  ix   , iy, nxEC);
	GShearW = shearValue(Physics->G, (ix-1), iy, nxEC);


	ZN = (dt*Physics->G     [NormalN]) / (dt*Physics->G     [NormalN] + EtaN);
	ZS = (dt*Physics->G     [NormalS]) / (dt*Physics->G     [NormalS] + EtaS);
	ZE = (dt*GShearE) / (dt*GShearE + EtaE);
	ZW = (dt*GShearW) / (dt*GShearW + EtaW);

	sigma_yy_0_N = -Physics->sigma_xx_0[NormalN];
	sigma_yy_0_S = -Physics->sigma_xx_0[NormalS];
	sigma_xy_0_E =  Physics->sigma_xy_0[ShearE ];
	sigma_xy_0_W =  Physics->sigma_xy_0[ShearW ];



	Jloc[order[ 0]] =   ix      + (iy  )*nxVx - 1     + VxPeriod              ; // VxSW
	Jloc[order[ 1]] =   ix      + (iy  )*nxVx                       ; // VxSE
	Jloc[order[ 2]] =   ix      + (iy+1)*nxVx - 1     + VxPeriod              ; // VxNW
	Jloc[order[ 3]] =   ix      + (iy+1)*nxVx                       ; // VxNE
	Jloc[order[ 4]] =   ix      + iy*nxVy     + nVxTot    - nxVy    ; // VyS
	Jloc[order[ 5]] =   ix      + iy*nxVy     + nVxTot   - 1     + VyPeriod   ; // VyW
	Jloc[order[ 6]] =   ix      + iy*nxVy     + nVxTot              ; // VyC
	Jloc[order[ 7]] =   ix      + iy*nxVy     + nVxTot    + 1       ; // VyE
	Jloc[order[ 8]] =   ix      + iy*nxVy     + nVxTot    + nxVy    ; // VyN
	Jloc[order[ 9]] =   ix-1    + (iy-1)*nxN + nVxTot+nVyTot   + PPeriod     ; // PS
	Jloc[order[10]] =   ix-1    + (iy  )*nxN + nVxTot+nVyTot   + PPeriod     ; // PN



	// Fill Vloc: list of coefficients
	// ================================
	Vloc[order[ 0]] =  EtaW*ZW/dxW/dyS; // VxSW
	Vloc[order[ 1]] = -EtaE*ZE/dxE/dyS; // VxSE
	Vloc[order[ 2]] = -EtaW*ZW/dxW/dyN; // VxNW
	Vloc[order[ 3]] =  EtaE*ZE/dxE/dyN; // VxNE
	Vloc[order[ 4]] =  2.0 * EtaS*ZS/dyS/dyC; // VyS
	Vloc[order[ 5]] =  EtaW*ZW/dxW/dxC; 		 //VyW
	Vloc[order[ 6]] = -2.0 * EtaN*ZN/dyN/dyC   -2.0 * EtaS*ZS/dyS/dyC   -1.0 * EtaE*ZE/dxE/dxC   -1.0 * EtaW*ZW/dxW/dxC; // VyC
	Vloc[order[ 7]] =  EtaE*ZE/dxE/dxC; // VyE
	Vloc[order[ 8]] =  2.0 * EtaN*ZN/dyN/dyC; //VyN
	Vloc[order[ 9]] =  1.0/dyS; // PS
	Vloc[order[10]] = -1.0/dyN; // PN

	*bloc = - Physics->g[1] * 0.5 * ( Physics->rho[NormalN] + Physics->rho[NormalS] );

	// add contributions of old stresses
	*bloc += - (sigma_yy_0_N*(1-ZN) - sigma_yy_0_S*(1-ZS))/dyC  -  (sigma_xy_0_E*(1-ZE) - sigma_xy_0_W*(1-ZW))/dxC;

}







static void Static_LocalStencilP(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift)
{


	// Define size variables
	// ===============================
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x

	int nVxTot = nxVx*nyVx;

	int nxN = Grid->nxC;

	compute dx = Grid->DXS[ix];
	compute dy = Grid->DYS[iy];



	// Maximum number of non zeros for Stokes on the staggered grid
	if (UPPER_TRI) {
		*shift = 4;
	}
	else {
		*shift = 0;
	}

	// =====================================================================
	//                               locJ
	// =====================================================================

	// Fill Jloc: list of all J indices (including Dirichlet)
	// ================================================================
	if (SetupType==SimpleShearPeriodic) {

		if (ix==nxN-1) {
			if (UPPER_TRI) {
				*shift = 4;
			}
			order[0] = 1;
			order[1] = 0;
			order[2] = 2;
			order[3] = 3;
		}
	}


	Jloc[order[0]] = ix   + (iy+1)*nxVx              ; // VxW
	Jloc[order[1]] = ix+1 + (iy+1)*nxVx              ; // VxE
	Jloc[order[2]] = ix+1 + iy*(nxVy)      + nVxTot  ; // VyS
	Jloc[order[3]] = ix+1 + (iy+1)*(nxVy)  + nVxTot  ; // VyN

	// =====================================================================
	//                               locV
	// =====================================================================
	// Fill Vloc: list of coefficients
	// ================================
	Vloc[order[0]] = -1.0/dx;
	Vloc[order[1]] =  1.0/dx;
	Vloc[order[2]] = -1.0/dy;
	Vloc[order[3]] =  1.0/dy;



	*bloc = 0;
}

#if (HEAT)
static void Static_LocalStencilT(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift)
{



	int TPeriod = 0;

	int TS, TW, TC, TE, TN;
	// Define size variables
	// ===============================
	int nxEC = Grid->nxEC;


	compute kN, kS, kW, kE;
	compute dxW = Grid->DXEC[ix-1];
	compute dxE = Grid->DXEC[ix];
	compute dxC = 0.5*(dxW+dxE);
	compute dyS = Grid->DYEC[iy-1];
	compute dyN = Grid->DYEC[iy];
	compute dyC = 0.5*(dyS+dyN);



	compute dt = Physics->dtT;

	if (UPPER_TRI) {
		*shift = 2;
	}
	else {
		*shift = 0;
	}


	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 1;
			}
			TPeriod  = (nxEC)-2  ; // VyW

			order[ 0] =  0; // TS
			order[ 1] =  3; // TW
			order[ 2] =  1; // TC
			order[ 3] =  2; // TE
			order[ 4] =  4; // TN

		}
		else if (ix==nxEC-3) {
			if (UPPER_TRI) {
				*shift = 3;
			}
			order[ 0] =  0; // TS
			order[ 1] =  2; // TW
			order[ 2] =  3; // TC
			order[ 3] =  1; // TE
			order[ 4] =  4; // TN
		}
	}


	TS =  ix 	+ (iy-1)*(nxEC);
	TW = (ix-1) +  iy   *(nxEC) + TPeriod;
	TC =  ix 	+  iy   *(nxEC);
	TE = (ix+1) +  iy   *(nxEC) ;
	TN =  ix 	+ (iy+1)*(nxEC);


	Jloc[order[0]] = TS;
	Jloc[order[1]] = TW;
	Jloc[order[2]] = TC;
	Jloc[order[3]] = TE;
	Jloc[order[4]] = TN;


	kN = (2*Physics->k[TN]*Physics->k[TC])/(Physics->k[TN]+Physics->k[TC]); // harmonic average
	kS = (2*Physics->k[TS]*Physics->k[TC])/(Physics->k[TS]+Physics->k[TC]);
	kW = (2*Physics->k[TW]*Physics->k[TC])/(Physics->k[TW]+Physics->k[TC]);
	kE = (2*Physics->k[TE]*Physics->k[TC])/(Physics->k[TE]+Physics->k[TC]);

	Vloc[order[0]] =  -kS/dyS/dyC; // TS
	Vloc[order[1]] =  -kW/dxW/dxC; // TW
	Vloc[order[2]] =  -(-kW/dxW/dxC -kE/dxE/dxC -kN/dyN/dyC -kS/dyS/dyC) + Physics->rho[TC]*Physics->Cp/dt; // TC
	Vloc[order[3]] =  -kE/dxE/dxC; // TE
	Vloc[order[4]] =  -kN/dyN/dyC; // TN


	*bloc = + Physics->rho[TC]*Physics->Cp*Physics->T[TC]/dt;
}
#endif






















void EqSystem_check(EqSystem* EqSystem)
{
	int i, j, I;
	// Check

	printf(" ===== Isparse =====\n");
	for (i=0;i<EqSystem->nEq+1;i++) {
		printf("%i  ", EqSystem->I[i]);
	}
	printf(" \n");


	printf(" ===== Jsparse =====\n");
	for (i=0;i<EqSystem->nnz;i++) {
		printf("%i  ", EqSystem->J[i]);
	}
	printf(" \n");

	printf(" ===== Vsparse =====\n");
	printf("%.2e  \n", EqSystem->V[0]);
	for (i=0;i<EqSystem->nnz;i++) {
		printf("%.2e  ", EqSystem->V[i]);
	}
	printf(" \n");

	printf(" ===== b=====\n");
	for (i=0;i<EqSystem->nEq;i++) {
		printf("%.2e  ", EqSystem->b[i]);
	}
	printf(" \n");


	printf("===== SPY =====\n");
	printf("   ");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("%*i ",2,i);
	}
	printf("\n");
	int padding;
	for (i=0; i<EqSystem->nRow; i++) {
		I = EqSystem->I[i];
		printf("%*i ",2,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			if (j==0){
				padding = EqSystem->J[I+j]+1;
			}
			else{
				padding = EqSystem->J[I+j]-EqSystem->J[I+j-1];
			}
			printf("%*s",3*padding,"x ");
		}
		printf("\n");
	}

	// List J per row
	printf("===== J per row =====\n");
	//for (i=0; i<EqSystem->nEq; i++) {
	for (i=0; i<EqSystem->nRow; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*i",7,EqSystem->J[I+j]);
		}
		printf("\n");
	}
	printf("\n");


	printf("===== V per row=====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		I = EqSystem->I[i];
		printf("row #%*i :",3,i);
		for (j=0; j<EqSystem->I[i+1]-EqSystem->I[i]; j++) {
			printf("%*.2f",10,EqSystem->V[I+j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("===== RHS =====\n");
	for (i=0; i<EqSystem->nEq; i++) {
		printf("RHS[%i] = %.2f\n", i, EqSystem->b[i]);
	}
	printf("\n");

}






void EqSystem_solve(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering)
{
	//int i;
	INIT_TIMER
	TIC


	// Reinitialize Pressure
	int iCell;
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		Physics->P[iCell] = 0;
	}



	if (UPPER_TRI) {
		pardisoSolveSymmetric(EqSystem, Solver, Grid, Physics, BC, Numbering);
	}
	else {
		printf("No solver function for assymmetric matrices\n");
		exit(0);
	}


	TOC

	printf("Direct solve: %.2f s\n", toc);



}










void EqSystem_initSolver (EqSystem* EqSystem, Solver* Solver)
{

	//int *ia ,int *ja ,compute *a ,compute *x ,compute *b, int n
	printf("===== Init Solver =====\n");
	INIT_TIMER
	TIC
	int i;


	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0;
	 	EqSystem->b[i] = 0;
	}
	for (i=0; i<EqSystem->nnz; i++) {
		EqSystem->V[i] = 0;
	}



	Solver->mtype = -2;        /* Real symmetric matrix */

	Solver->nrhs = 1;          /* Number of right hand sides. */

	/* Internal solver memory pointer pt,                  */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
	/* or void *pt[64] should be OK on both architectures  */
	//void    *pt[64];

	/* Pardiso control parameters. */
	//int      iparm[64];
	//double   dparm[64];
	//int      maxfct, mnum, phase, error, msglvl, solver;

	int phase;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters.                                */
	/* -------------------------------------------------------------------- */

	int error = 1;
	int solver = 0;// use sparse direct solver
	Solver->msglvl = 0;         // Print statistical information
	pardisoinit (Solver->pt,  &Solver->mtype, &solver, Solver->iparm, Solver->dparm, &error);

	if (error != 0)
	{
		if (error == -10 )
			printf("No license file found \n");
		if (error == -11 )
			printf("License is expired \n");
		if (error == -12 )
			printf("Wrong username or hostname \n");
		exit(0);
	}
	else {
		printf("[PARDISO]: License check was successful ... \n");
	}





	// Numbers of processors, value of OMP_NUM_THREADS
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1");
		exit(1);
	}
	printf("Number of procs: %i\n", num_procs);


	// Solver options
	Solver->iparm[2]  = num_procs;

	Solver->iparm[32] = 0; /* compute determinant */
	Solver->iparm[7]  = 1; /* Max numbers of iterative refinement steps. */

	Solver->iparm[27] = 1; // 0: sequential reordering, 1: parallel reordering in METIS
	//Solver->iparm[50] = 0; // 0: openMP, 1:MPI
	//Solver->iparm[51] = 0; // number of compute nodes for MPI


	Solver->iparm[29] = 100; //size of supernodes, default 80

	Solver->maxfct = 1;		// Maximum number of numerical factorizations.
	Solver->mnum   = 1;     // Which factorization to use.

	Solver->msglvl = 0;     // Print statistical information
	error  = 0;         	// Initialize error flag



	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */

	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] += 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] += 1;
	}




	/* -------------------------------------------------------------------- */
	/*  .. pardiso_chk_matrix(...)                                          */
	/*     Checks the consistency of the given matrix.                      */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */

	if  (DEBUG) {
		printf("--  chkmatrix\n");
		pardiso_chkmatrix  (&Solver->mtype, &EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &error);
		if (error != 0) {
			printf("\nERROR in consistency of matrix: %d", error);
			exit(1);
		}
	}

	/* -------------------------------------------------------------------- */
	/* ..  pardiso_chkvec(...)                                              */
	/*     Checks the given vectors for infinite and NaN values             */
	/*     Input parameters (see PARDISO user manual for a description):    */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */

	if  (DEBUG) {
		printf("--  chkvec\n");
		pardiso_chkvec (&EqSystem->nEq, &Solver->nrhs, EqSystem->b, &error);
		if (error != 0) {
			printf("\nERROR  in right hand side: %d", error);
			exit(1);
		}
	}


	/* -------------------------------------------------------------------- */
	/* .. pardiso_printstats(...)                                           */
	/*    prints information on the matrix to STDOUT.                       */
	/*    Use this functionality only for debugging purposes                */
	/* -------------------------------------------------------------------- */

	if  (DEBUG) {
		printf("--  printstats\n");
		pardiso_printstats (&Solver->mtype, &EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &Solver->nrhs, EqSystem->b, &error);
		if (error != 0) {
			printf("\nERROR right hand side: %d", error);
			exit(1);
		}
	}

	/* -------------------------------------------------------------------- */
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */

	phase = 11;

	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error, Solver->dparm);


	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors  = %d", Solver->iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d\n", Solver->iparm[18]);









	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] -= 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] -= 1;
	}

	TOC
	printf("Solver initialization: %.3f s\n", toc);



}


void pardisoSolveSymmetric(EqSystem* EqSystem, Solver* Solver, Grid* Grid, Physics* Physics, BC* BC, Numbering* Numbering)
{



	INIT_TIMER
	int i, phase;
	double   	ddum;              // Double dummy
	int      	idum;              // Integer dummy.
	int 		error;


	for (i=0; i<EqSystem->nEq; i++) {
		EqSystem->x[i] = 0;
	}

	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */


	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] += 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] += 1;
	}



	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */

	if (TIMER) {
		TIC
	}

	phase = 22;

	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, &ddum, &ddum, &error,  Solver->dparm);

	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//printf("Factorization completed ...\n ");

	if (TIMER) {
		TOC
		printf("Phase 22 - Numerical factorization: %.3f s\n", toc);
	}



	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	if (TIMER) {
		TIC
	}


	phase = 33;

	// Solve full system Vx, Vy, P


	pardiso (Solver->pt, &Solver->maxfct, &Solver->mnum, &Solver->mtype, &phase,
			&EqSystem->nEq, EqSystem->V, EqSystem->I, EqSystem->J, &idum, &Solver->nrhs,
			Solver->iparm, &Solver->msglvl, EqSystem->b, EqSystem->x, &error,  Solver->dparm);




	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}







	if (TIMER) {
		TOC
		printf("Phase 33 - Back substitution: %.3f s\n", toc);
	}
	if  (DEBUG) {
		printf("\nThe solution of the system is: \n");

		for (i = 0; i < EqSystem->nEq; i++) {
			printf(" x [%d] = % f\n", i, EqSystem->x[i] );
		}
	}






	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < EqSystem->nEq+1; i++) {
		EqSystem->I[i] -= 1;
	}
	for (i = 0; i < EqSystem->nnz; i++) {
		EqSystem->J[i] -= 1;
	}







}




void EqSystem_computeNormResidual(EqSystem* EqSystem)
{
	compute* Residual = (compute*) malloc(EqSystem->nEq * sizeof(compute));

	int iEq;
	int J,i;
	EqSystem->normResidual = 0;


#pragma omp parallel for private(iEq, i, J) schedule(static,32)
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		Residual[iEq] = EqSystem->b[iEq];
		for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
			J = EqSystem->J[i];
			Residual[iEq] += - (EqSystem->V[i]*EqSystem->x[J]);
			/*
			if (UPPER_TRI) {
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}*/
		}
	}

	if (UPPER_TRI) {

#pragma omp parallel for private(iEq, i, J) schedule(static,32)
		for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
			for (i = EqSystem->I[iEq]; i < EqSystem->I[iEq+1]; ++i) {
				J = EqSystem->J[i];
				if (J!=iEq)
					Residual[J] += - (EqSystem->V[i]*EqSystem->x[iEq]);// Wrong
			}
		}

	}

	compute norm_b = 0;
	for (iEq = 0; iEq < EqSystem->nEq; ++iEq) {
		EqSystem->normResidual += Residual[iEq]*Residual[iEq];
		norm_b += EqSystem->b[iEq]*EqSystem->b[iEq];
	}
	norm_b = sqrt(norm_b);
	EqSystem->normResidual = sqrt(EqSystem->normResidual);


	// Normalize the residual
	EqSystem->normResidual /= norm_b; // Normalize the norm of the residual by the norm of the right hand side



	free(Residual);
}

