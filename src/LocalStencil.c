/*
 * LocalStencil.c
 *
 *  Created on: Jul 27, 2016
 *      Author: abauville
 */


#include "stokes.h"



void LocalStencil_Call(StencilType Stencil, int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{
	/*
	 * LocalStencil_Call switches between the different LocalStencil functions
	 *
	 *
	 */

		if (Stencil==Stencil_Stokes_Momentum_x)		{
			LocalStencil_Stokes_Momentum_x(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
		else if (Stencil==Stencil_Stokes_Momentum_y) 	{
			LocalStencil_Stokes_Momentum_y(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
		else if (Stencil==Stencil_Stokes_Continuity) 	{
			LocalStencil_Stokes_Continuity(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
#if (DARCY)
		if (Stencil==Stencil_Stokes_Darcy_Momentum_x)		{
			LocalStencil_Stokes_Darcy_Momentum_x(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
		else if (Stencil==Stencil_Stokes_Darcy_Momentum_y) 	{
			LocalStencil_Stokes_Darcy_Momentum_y(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
		else if (Stencil==Stencil_Stokes_Darcy_Continuity) 	{
			LocalStencil_Stokes_Darcy_Continuity(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
		else if (Stencil==Stencil_Stokes_Darcy_Darcy) 	{
			LocalStencil_Stokes_Darcy_Darcy(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
#endif

#if (HEAT)
		else if (Stencil==Stencil_Heat) 	{
			LocalStencil_Heat(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
		}
#endif

}


void LocalStencil_Stokes_Momentum_x(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 11;
	*IC = 2;

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

	int nxN = Grid->nxEC;
	int nxEC = Grid->nxEC;

	int nxS = Grid->nxS;

	compute EtaN, EtaS, EtaE, EtaW;
	compute ZN, ZS, ZW, ZE; // visco-elasticity factor


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);


	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXS[0]+Grid->DXS[nxS-2]);
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVx-1) {
			dxW = Grid->DXS[ix-1];
			dxE = 0.5*(Grid->DXS[0]+Grid->DXS[nxS-2]);
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
			PPeriod  = 0;//nxN;

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
			order[ 9] =  9; // PW
			order[10] = 10; // PE
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
			order[ 9] = 10; // PW
			order[10] =  9; // PE
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
	Jloc[order[ 9]]   =   ix    + (iy)*nxN  + nVxTot+nVyTot + PPeriod; // PW
	Jloc[order[10]]   =   ix+1  + (iy)*nxN  + nVxTot+nVyTot ; // PE



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

	ZN = (dt*GShearN) / (dt*GShearN + EtaN);
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









void LocalStencil_Stokes_Momentum_y(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 11;
	*IC = 6;

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

	int nxN = Grid->nxEC;
	int nxEC = Grid->nxEC;

	int nxS = Grid->nxS;

	compute EtaN, EtaS, EtaE, EtaW;
	compute ZN, ZS, ZE, ZW;
	compute sigma_yy_0_N, sigma_yy_0_S, sigma_xy_0_E, sigma_xy_0_W;
	compute GShearE, GShearW;



	compute dyS = Grid->DYS [iy-1];
	compute dyN = Grid->DYS [iy  ];
	compute dyC = 0.5*(dyS+dyN);
	compute dt = Physics->dt;

	compute dxW, dxE, dxC;
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
			dxE = Grid->DXEC[ix  ];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVx-1) {
			dxW = Grid->DXEC[ix-1];
			dxE = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
			dxC = 0.5*(dxW+dxE);

		} else {
			dxW = Grid->DXS[ix-1];
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);
		}

	} else {
		dxW = Grid->DXEC[ix-1];
		dxE = Grid->DXEC[ix  ];
		dxC = 0.5*(dxW+dxE);

	}




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
		if (ix==nxVy-2) {
			if (UPPER_TRI) {
				*shift = 5;
			}
			VxPeriod =  0;//nxVx-1		; // VxSW
			VyPeriod  = 0;//nxVy-2  	; // VyW
			PPeriod   = 0;//nxN    		; // PS
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
	Jloc[order[ 9]] =   ix      + (iy)*nxN + nVxTot+nVyTot   + PPeriod     ; // PS
	Jloc[order[10]] =   ix      + (iy+1)*nxN + nVxTot+nVyTot   + PPeriod     ; // PN



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







void LocalStencil_Stokes_Continuity(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 4;
	*IC = 0;

	// Define size variables
	// ===============================
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x

	int nVxTot = nxVx*nyVx;

	int nxN = Grid->nxEC;

	compute dx = Grid->DXS[ix-1];
	compute dy = Grid->DYS[iy-1];

		if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				dx = 0.5*(Grid->DXS[0]+Grid->DXS[Grid->nxS-2]);

			}  else {
				compute dx = Grid->DXS[ix-1];
			}

		} else {
			compute dx = Grid->DXS[ix-1];
		}


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

		if (ix==nxN-2) {
			if (UPPER_TRI) {
				*shift = 4;
			}
			order[0] = 1;
			order[1] = 0;
			order[2] = 2;
			order[3] = 3;
		}
	}


	Jloc[order[0]] = ix-1   + (iy)*nxVx              ; // VxW
	Jloc[order[1]] = ix + (iy)*nxVx              ; // VxE
	Jloc[order[2]] = ix + (iy-1)*(nxVy)      + nVxTot  ; // VyS
	Jloc[order[3]] = ix + (iy)*(nxVy)  + nVxTot  ; // VyN

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
void LocalStencil_Heat(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 5;
	*IC = 2;

	int TPeriod = 0;

	int TS, TW, TC, TE, TN;
	// Define size variables
	// ===============================
	int nxEC = Grid->nxEC;


	compute kN, kS, kW, kE;
	compute dxW, dxE, dxC;
	compute dyS = Grid->DYEC[iy-1];
	compute dyN = Grid->DYEC[iy];
	compute dyC = 0.5*(dyS+dyN);

	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[0]+Grid->DXEC[Grid->nxS-1]);
			dxE = Grid->DXEC[ix];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==Grid->nxEC-1) {
			dxW = Grid->DXEC[ix-1];
			dxE = 0.5*(Grid->DXEC[0]+Grid->DXEC[Grid->nxS-1]);
			dxC = 0.5*(dxW+dxE);

		} else {
			dxW = Grid->DXEC[ix-1];
			dxE = Grid->DXEC[ix  ];
			dxC = 0.5*(dxW+dxE);
		}

	} else {
		dxW = Grid->DXEC[ix-1];
		dxE = Grid->DXEC[ix  ];
		dxC = 0.5*(dxW+dxE);

	}

	/*
	dxW = Grid->dx;
	dxE = Grid->dx;
	dxC = 0.5*(dxW+dxE);
	*/

	//printf("ix = %i, iy = %i, dx = %.2f, dxW = %.2f, dxE = %.2f, dxC = %.2f \n", ix, iy,  Grid->dx, dxW, dxE, dxC);



	compute dt = Physics->dtT;

	if (UPPER_TRI) {
		*shift = 2;
	}
	else {
		*shift = 0;
	}


	if (SetupType==SimpleShearPeriodic) {
		if (ix==nxEC-2) {
			if (UPPER_TRI) {
				*shift = 1;
			}
			TPeriod  = 0;//(nxEC)-2  ; // VyW

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
	TE = (ix+1) +  iy   *(nxEC);
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

	// Add the contribution of the shear heating
	compute EII, sigma_xy, sigma_xx, sigmaII;
	Physics_computeStrainInvariantForOneCell(Physics, Grid, ix, iy, &EII);
	sigma_xy  = Physics->sigma_xy_0[ix-1 + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix   + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix-1 + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix   + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
	sigma_xy /= 4.0;

	sigma_xx = Physics->sigma_xx_0[ix+iy*nxEC] + Physics->Dsigma_xx_0[ix+iy*nxEC];
	sigmaII = sqrt(sigma_xx*sigma_xx + sigma_xy*sigma_xy);

	*bloc += sigmaII*EII;

}
#endif






#if (DARCY)


void LocalStencil_Stokes_Darcy_Momentum_x(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	// 1. call Stokes_Momentum_x
	LocalStencil_Stokes_Momentum_x(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);


	/*
	// 2. multiply by 2 the contribution from Pt
	Vloc[order[ 9]] += Vloc[order[ 9]]; // PW
	Vloc[order[10]] += Vloc[order[10]]; //PE
	*/


	*nLoc = 13;



	int nxN = Grid->nxC;
	int nxS = Grid->nxS;
	int PPeriod = 0;
	int nCTot = Grid->nCTot;
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y
	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;

	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXS[0]+Grid->DXS[nxS-2]);
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVx-1) {
			dxW = Grid->DXS[ix-1];
			dxE = 0.5*(Grid->DXS[0]+Grid->DXS[nxS-2]);
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





	// 3. add contributions from Pf

	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 1;
			}
			PPeriod  = nxN;


			/*
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
			*/
			order[11] = 12; // PfW
			order[12] = 11; // PfE

		}
		if (ix==Grid->nxVx-2) {
			if (UPPER_TRI) {
				*shift = 3;
			}
			/*
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
			*/

			order[11] = 11; // PfW
			order[12] = 12; // PfE

		}
	}

	Jloc[order[11]]   =   ix-1    + (iy-1)*nxN  + nVxTot+nVyTot+nCTot + PPeriod; // PfW
	Jloc[order[12]]   =   ix      + (iy-1)*nxN  + nVxTot+nVyTot+nCTot ; // PfE

	Vloc[order[11]] =  1.0/dxW;
	Vloc[order[12]] = -1.0/dxE;



}












void LocalStencil_Stokes_Darcy_Momentum_y(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	// 1. call Stokes_Momentum_x
	LocalStencil_Stokes_Momentum_y(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);



	/*
	// 2. multiply by 2 the contribution from Pt
	Vloc[order[ 9]] += Vloc[order[ 9]]; // PS
	Vloc[order[10]] += Vloc[order[10]]; // PN
	*/


	*nLoc = 13;



	int nxN = Grid->nxC;
	int nxS = Grid->nxS;
	int PPeriod = 0;
	int nCTot = Grid->nCTot;
	int nxEC = Grid->nxEC;
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y
	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYS [iy-1];
	compute dyN = Grid->DYS [iy  ];
	compute dyC = 0.5*(dyS+dyN);
	compute dt = Physics->dt;



	if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				dxW = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
				dxE = Grid->DXEC[ix  ];
				dxC = 0.5*(dxW+dxE);

			} else if (ix==nxVx-1) {
				dxW = Grid->DXEC[ix-1];
				dxE = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
				dxC = 0.5*(dxW+dxE);

			} else {
				dxW = Grid->DXS[ix-1];
				dxE = Grid->DXS[ix];
				dxC = 0.5*(dxW+dxE);
			}

		} else {
			dxW = Grid->DXEC[ix-1];
			dxE = Grid->DXEC[ix  ];
			dxC = 0.5*(dxW+dxE);

		}






	// 3. add contributions from Pf


	// Special cases for periodic BC
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			PPeriod   = nxN    		; // PS
		}
		else if (ix==nxVy-3) {

		}
	}

	Jloc[order[11]] =   ix-1    + (iy-1)*nxN + nVxTot+nVyTot+nCTot   + PPeriod     ; // PcS
	Jloc[order[12]] =   ix-1    + (iy  )*nxN + nVxTot+nVyTot+nCTot   + PPeriod     ; // PcN

	Vloc[order[11]] =  1.0/dyS; // PcS
	Vloc[order[12]] = -1.0/dyN; // PcN



}










void LocalStencil_Stokes_Darcy_Continuity(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{
	*bloc = 0; // just a security
	// 1. call Stokes_Momentum_x to build the veolocity divergence
	LocalStencil_Stokes_Continuity(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);

	*nLoc = 5;
	*IC = 5;


	if (UPPER_TRI) {
		*shift = 4;
	}
	else {
		*shift = 0;
	}

	int nxC = Grid->nxC;
	int nVxTot = Grid->nVxTot;
	int nVyTot = Grid->nVyTot;
	int nCTot = Grid->nCTot;

	int NormalC = ix+1 + (iy+1)*Grid->nECTot; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells

	compute dt = Physics->dt;
	compute Zb, ZbStar;

	Jloc[order[4]] = ix    + (iy)*nxC + nVxTot+nVyTot+nCTot; // PcC

	compute eta_b =  Physics->eta_b[NormalC];
	Zb = (dt*Physics->G     [NormalC]) / (dt*Physics->B     [NormalC] + eta_b);
	ZbStar = (1 - Physics->phi[NormalC]) * Zb;
	Vloc[order[4]] =  1.0/ (ZbStar * eta_b);
	//Vloc[order[6]] = -1.0/ (ZbStar * eta_b);

	*bloc += -     (        (1.0 - Zb)*Physics->Pc0[NormalC]       )       /      (  ZbStar*eta_b  )   ;




}












void LocalStencil_Stokes_Darcy_Darcy 	 (int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{



	*bloc = 0; // just a security
	// 1. call Stokes_Momentum_x to build the velocity divergence
	LocalStencil_Stokes_Continuity(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
	*bloc = 0; // just a security

	*nLoc = 9;
	*IC = 6;



	int nxC = Grid->nxC;
	int nxEC = Grid->nxEC;
	int nCTot = Grid->nCTot;
	int nVxTot = Grid->nVxTot;
	int nVyTot = Grid->nVyTot;


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);


	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[ix+1]+Grid->DXEC[nxEC-3]);
			dxE = Grid->DXEC[ix+1];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxC-1) {
			dxW = Grid->DXEC[ix+1-1];
			dxE = 0.5*(Grid->DXEC[ix+1]+Grid->DXS[nxEC-3]);
			dxC = 0.5*(dxW+dxE);

		} else {
			dxW = Grid->DXEC[ix+1-1];
			dxE = Grid->DXEC[ix+1];
			dxC = 0.5*(dxW+dxE);
		}

	} else {
		dxW = Grid->DXEC[ix+1-1];
		dxE = Grid->DXEC[ix+1];
		dxC = 0.5*(dxW+dxE);

	}






	if (UPPER_TRI) {
		*shift = 6;
	}
	else {
		*shift = 0;
	}


	int PPeriodL = 0;
	int PPeriodR = 0;

	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 5;
			}
			PPeriodL = nxC;
			PPeriodR = 0;
			order[4] = 4; // PfS
			order[5] = 7; // PfW
			order[6] = 5; // PfC
			order[7] = 6; // PfE
			order[8] = 8; // PfN
		}
		else if (ix == nxC-1) {
			if (UPPER_TRI) {
				*shift = 7;
			}
			PPeriodL = 0;
			PPeriodR = -nxC;
			order[4] = 4; // PfS
			order[5] = 6; // PfW
			order[6] = 7; // PfC
			order[7] = 5; // PfE
			order[8] = 8; // PfN
		}
	}



	// the first 4 slots are filled by Vx and Vy equations
	Jloc[order[4]] = ix    + (iy-1)*nxC + nVxTot+nVyTot; // PfS
	Jloc[order[5]] = ix-1  + (iy  )*nxC + nVxTot+nVyTot + PPeriodL; // PfW
	Jloc[order[6]] = ix    + (iy  )*nxC + nVxTot+nVyTot; // PfC
	Jloc[order[7]] = ix+1  + (iy  )*nxC + nVxTot+nVyTot + PPeriodR; // PfE
	Jloc[order[8]] = ix    + (iy+1)*nxC + nVxTot+nVyTot; // PfN


	int NormalS = ix+1   + (iy+1-1)*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalW = ix+1-1 + (iy+1  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalC = ix+1   + (iy+1  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalE = ix+1+1 + (iy+1  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalN = ix+1   + (iy+1+1)*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells

	compute KS 		= ((Physics->perm[NormalS] + Physics->perm[NormalC])/2.0) / (Physics->eta_f); // averaging because K has to be defined on the shear node
	compute KN 		= ((Physics->perm[NormalN] + Physics->perm[NormalC])/2.0) / (Physics->eta_f); // averaging because K has to be defined on the shear node
	compute KW	 	= ((Physics->perm[NormalW] + Physics->perm[NormalC])/2.0) / (Physics->eta_f); // averaging because K has to be defined on the shear node
	compute KE 		= ((Physics->perm[NormalE] + Physics->perm[NormalC])/2.0) / (Physics->eta_f); // averaging because K has to be defined on the shear node

	compute rhoS 	=  (Physics->rho [NormalS] + Physics->rho [NormalC])/2.0;
	compute rhoN 	=  (Physics->rho [NormalN] + Physics->rho [NormalC])/2.0;
	compute rhoW 	=  (Physics->rho [NormalW] + Physics->rho [NormalC])/2.0;
	compute rhoE 	=  (Physics->rho [NormalE] + Physics->rho [NormalC])/2.0;

	Vloc[order[4]] = KS/dyS/dyC; // PfS
	Vloc[order[5]] = KW/dxE/dxC; // PfW
	Vloc[order[6]] = -KE/dxE/dxC -KW/dxW/dxC -KS/dyS/dyC -KN/dyN/dyC; // PfC
	Vloc[order[7]] = KE/dxE/dxC; // PfE
	Vloc[order[8]] = KN/dyN/dyC; // PfN

	*bloc += KE*rhoE*Physics->g[0]/dxE/dxC - KW*rhoW*Physics->g[0]/dxW/dxC;
	*bloc += KN*rhoN*Physics->g[1]/dyN/dyC - KS*rhoS*Physics->g[1]/dyS/dyC;



}



























#endif


































