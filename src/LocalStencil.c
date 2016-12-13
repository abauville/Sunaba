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

	//compute EtaN, EtaS, EtaE, EtaW;
	//compute KhiN, KhiS, KhiE, KhiW;
	compute GN  , GS  , GE  , GW  ;
	compute ZN, ZS, ZE, ZW; // visco-elasticity factor
	//compute ZN, ZS, ZW, ZE; // visco-elasticity factor
	//compute One_ZN, One_ZS, One_ZE, One_ZW;


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);





	if (Grid->isPeriodic) {
		if (ix==0) {
			dxW = Grid->DXS[nxS-2];
			dxE = Grid->DXS[ix];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVx-1) {
			dxW = Grid->DXS[ix-1];
			dxE = Grid->DXS[0];
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



	//printf("dxW = %.2f, dyS = %.2f, Grid->dx = %.2f, ix = %i, iy = %i,Grid->nxVx = %i, Grid->nxS = %i\n",dxW, dyS,Grid->dx,ix,iy,Grid->nxVx, Grid->nxS);

	if (UPPER_TRI) {
		*shift = 2;
	}
	else {
		*shift = 0;
	}

	// Special case for periodic BC
	if (Grid->isPeriodic) {
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



	NormalE = ix+1   + (iy)*nxEC;
	NormalW = ix     + (iy)*nxEC + NormalPeriod;
	ShearN = ix      + iy*nxS;
	ShearS = ix      + (iy-1)*nxS;

	////EtaN    = shearValue(Physics->eta, ix,  iy   , nxEC); // Shear N
	////EtaS    = shearValue(Physics->eta, ix, (iy-1), nxEC); // ShearS
	//EtaN 	= Physics->etaShear[ShearN];
	//EtaS 	= Physics->etaShear[ShearS];
	//EtaE    = Physics->eta[ NormalE ]; // NormalE
	//EtaW    = Physics->eta[ NormalW ]; // NormalW

	GN = shearValue(Physics->G, ix,  iy   , nxEC);
	GS = shearValue(Physics->G, ix, (iy-1), nxEC);
	GE = Physics->G[NormalE];
	GW = Physics->G[NormalW];

	//KhiN 	= Physics->khiShear[ ShearN ]; // NormalE
	//KhiS    = Physics->khiShear[ ShearS ]; // NormalE
	//KhiE    = Physics->khi[ NormalE ]; // NormalE
	//KhiW    = Physics->khi[ NormalW ]; // NormalE

	ZN = shearValue(Physics->Z, ix,  iy   , nxEC);
	ZS = shearValue(Physics->Z, ix, (iy-1), nxEC);
	ZE = Physics->Z[NormalE];
	ZW = Physics->Z[NormalW];


	//ZN  	= 1.0/( 1.0/KhiN + 1.0/EtaN + 1.0/(GN*dt) );
	//ZS  	= 1.0/( 1.0/KhiS + 1.0/EtaS + 1.0/(GS*dt) );
	//ZE  	= 1.0/( 1.0/KhiE + 1.0/EtaE + 1.0/(GE*dt) );
	//ZW  	= 1.0/( 1.0/KhiW + 1.0/EtaW + 1.0/(GW*dt) );



#if (DARCY)
	compute phiN, phiS, phiE, phiW;

	phiN    = shearValue(Physics->phi, ix,  iy   , nxEC); // Shear N
	phiS    = shearValue(Physics->phi, ix, (iy-1), nxEC); // ShearS
	phiE    = Physics->phi[ NormalE ]; // NormalE
	phiW    = Physics->phi[ NormalW ]; // NormalW

	ZN *= (1.0-phiN);
	ZS *= (1.0-phiS);
	ZE *= (1.0-phiE);
	ZW *= (1.0-phiW);

#endif

	sigma_xx_0_E =  Physics->sigma_xx_0[NormalE];
	sigma_xx_0_W =  Physics->sigma_xx_0[NormalW];
	sigma_xy_0_N =  Physics->sigma_xy_0[ShearN ];
	sigma_xy_0_S =  Physics->sigma_xy_0[ShearS ];

	// Fill Vloc: list of coefficients
	// ================================
	Vloc[order[ 0]] =  ZS/dyS/dyC;
	Vloc[order[ 1]] =  1.0 * ZW/dxW/dxC;
	Vloc[order[ 2]] = -1.0 * ZE/dxE/dxC   -1.0 * ZW/dxW/dxC   -1.0 * ZN/dyN/dyC   -1.0 * ZS/dyS/dyC;
	Vloc[order[ 3]] =  1.0 * ZE/dxE/dxC;
	Vloc[order[ 4]] =  ZN/dyN/dyC;

	Vloc[order[ 5]] =  ZS/dxC/dyC;
	Vloc[order[ 6]] = -ZS/dxC/dyC;
	Vloc[order[ 7]] = -ZN/dxC/dyC;
	Vloc[order[ 8]] =  ZN/dxC/dyC;
	Vloc[order[ 9]] =  1.0/dxC;
	Vloc[order[10]] = -1.0/dxC;

	*bloc = - Physics->gFac[0] * 0.5 * ( Physics->rho_g[NormalE] + Physics->rho_g[NormalW] );

	// add contributions of old stresses
	*bloc += - ( sigma_xx_0_E*ZE/(GE*dt)  -   sigma_xx_0_W*ZW/(GW*dt))/dxC  -  (sigma_xy_0_N*ZN/(GN*dt)  -  sigma_xy_0_S*ZS/(GS*dt))/dyC;
}









void LocalStencil_Stokes_Momentum_y(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 11;
	*IC = 6;

	int VxPeriod = 0;
	int VyPeriod = 0;
	int PPeriod  = 0;
	//int ShearPeriod = 0;



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

	//compute EtaN, EtaS, EtaE, EtaW;
	//compute KhiN, KhiS, KhiE, KhiW;
	compute GN  , GS  , GE  , GW  ;
	compute ZN, ZS, ZE, ZW; // visco-elasticity factor
	compute sigma_yy_0_N, sigma_yy_0_S, sigma_xy_0_E, sigma_xy_0_W;




	compute dyS = Grid->DYS [iy-1];
	compute dyN = Grid->DYS [iy  ];
	compute dyC = 0.5*(dyS+dyN);
	compute dt = Physics->dt;

	compute dxW, dxE, dxC;
	/*
	if (SetupType==SimpleShearPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
			dxE = Grid->DXEC[ix  ];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVy-1) {
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
	 */

	dxW = Grid->DXEC[ix-1];
	dxE = Grid->DXEC[ix  ];
	dxC = 0.5*(dxW+dxE);




	if (UPPER_TRI) {
		*shift = 6;
	}
	else {
		*shift = 0;
	}
	// =====================================================================
	//                          Local numbering
	// =====================================================================
	NormalN = ix      + (iy+1)*nxEC ;
	NormalS = ix      + (iy)*nxEC ;
	ShearE  = ix      + iy*nxS    ;
	ShearW  = ix-1    + iy*nxS    ;

	VxPeriod = 0		; // VxSW
	VyPeriod = 0 		; // VyW
	PPeriod  = 0   		; // PS
	//ShearPeriod = 0;

	// Special cases for periodic BC
	if (Grid->isPeriodic) {
		if (ix==nxVy-2) {
			if (UPPER_TRI) {
				*shift = 5;
			}
			VxPeriod =  0;//nxVx-1		; // VxSW
			VyPeriod  = 0;//nxVy-2  	; // VyW
			PPeriod   = 0;//nxN    		; // PS
			//ShearPeriod = 0;//nxS-1;


			NormalN += 0;//nxN			;
			NormalS += 0;//nxN			;
			ShearW  += 0;//nxS-1 		;

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
	//EtaN    = Physics->eta[NormalN];
	//EtaS    = Physics->eta[NormalS];

	////EtaE    = shearValue(Physics->eta,  ix   , iy, nxEC);
	////EtaW    = shearValue(Physics->eta, (ix-1)+ShearPeriod, iy, nxEC);
	//EtaE    = Physics->etaShear[ShearE];
	//EtaW    = Physics->etaShear[ShearW];

	GN = Physics->G[NormalN];
	GS = Physics->G[NormalS];
	GE = shearValue(Physics->G,  ix   , iy, nxEC);
	GW = shearValue(Physics->G, (ix-1), iy, nxEC);

	//KhiN = Physics->khi[NormalN];
	//KhiS = Physics->khi[NormalS];
	//KhiE = Physics->khiShear[ShearE];
	//KhiW = Physics->khiShear[ShearW];

	ZN = Physics->Z[NormalN];
	ZS = Physics->Z[NormalS];
	ZE = shearValue(Physics->Z,  ix   , iy, nxEC);
	ZW = shearValue(Physics->Z, (ix-1), iy, nxEC);


	//ZN  	= 1.0/( 1.0/KhiN + 1.0/EtaN + 1.0/(GN*dt) );
	//ZS  	= 1.0/( 1.0/KhiS + 1.0/EtaS + 1.0/(GS*dt) );
	//ZE  	= 1.0/( 1.0/KhiE + 1.0/EtaE + 1.0/(GE*dt) );
	//ZW  	= 1.0/( 1.0/KhiW + 1.0/EtaW + 1.0/(GW*dt) );


#if (DARCY)
	compute phiN, phiS, phiE, phiW;

	phiN = Physics->phi[NormalN];
	phiS = Physics->phi[NormalS];
	phiE = shearValue(Physics->phi,  ix   , iy, nxEC);
	//phiW = shearValue(Physics->phi, (ix-1)+ShearPeriod, iy, nxEC);
	phiW = shearValue(Physics->phi, (ix-1), iy, nxEC);
	ZN *= (1.0-phiN);
	ZS *= (1.0-phiS);
	ZE *= (1.0-phiE);
	ZW *= (1.0-phiW);

#endif



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
	Vloc[order[ 0]] =  ZW/dxC/dyC; // VxSW
	Vloc[order[ 1]] = -ZE/dxC/dyC; // VxSE
	Vloc[order[ 2]] = -ZW/dxC/dyC; // VxNW
	Vloc[order[ 3]] =  ZE/dxC/dyC; // VxNE
	Vloc[order[ 4]] =  1.0 * ZS/dyS/dyC; // VyS
	Vloc[order[ 5]] =  ZW/dxW/dxC; 		 //VyW
	Vloc[order[ 6]] = -1.0 * ZN/dyN/dyC   -1.0 * ZS/dyS/dyC   -1.0 * ZE/dxE/dxC   -1.0 * ZW/dxW/dxC; // VyC
	Vloc[order[ 7]] =  ZE/dxE/dxC; // VyE
	Vloc[order[ 8]] =  1.0 * ZN/dyN/dyC; //VyN
	Vloc[order[ 9]] =  1.0/dyC; // PS
	Vloc[order[10]] = -1.0/dyC; // PN

	*bloc = - Physics->gFac[1] * 0.5 * ( Physics->rho_g[NormalN] + Physics->rho_g[NormalS] );

	// add contributions of old stresses
	*bloc += - (sigma_yy_0_N*ZN/(GN*dt) - sigma_yy_0_S*ZS/(GS*dt))/dyC  -  (sigma_xy_0_E*ZE/(GE*dt) - sigma_xy_0_W*ZW/(GW*dt))/dxC;

}







void LocalStencil_Stokes_Continuity(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	*nLoc = 4;
	*IC = -1;

	// Define size variables
	// ===============================
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x

	int nVxTot = nxVx*nyVx;

	int nxN = Grid->nxEC;

	compute dx = Grid->DXS[ix-1];
	compute dy = Grid->DYS[iy-1];

	/*
		if (SetupType==SimpleShearPeriodic) {
			if (ix==0) {
				dx = 0.5*(Grid->DXS[0]+Grid->DXS[Grid->nxS-2]);

			}  else {
				dx = Grid->DXS[ix-1];
			}

		} else {
			dx = Grid->DXS[ix-1];
		}

		dx = Grid->DXS[ix-1];
		*/

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
	if (Grid->isPeriodic) {

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

	compute eta, rho;
	compute sigma_xy, sigma_xx;

	/*
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
	*/

	dxW = Grid->DXEC[ix-1];
	dxE = Grid->DXEC[ix  ];
	dxC = 0.5*(dxW+dxE);

	/*
	dxW = Grid->dx;
	dxE = Grid->dx;
	dxC = 0.5*(dxW+dxE);
	*/

	//printf("ix = %i, iy = %i, dx = %.2f, dxW = %.2f, dxE = %.2f, dxC = %.2f \n", ix, iy,  Grid->dx, dxW, dxE, dxC);



	compute dt = Physics->dt;

	if (UPPER_TRI) {
		*shift = 2;
	}
	else {
		*shift = 0;
	}


	if (Grid->isPeriodic) {
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

	rho = Physics->rho_g[TC]/sqrt(Physics->g[0]*Physics->g[0]+Physics->g[1]*Physics->g[1]);

	Vloc[order[0]] =  -kS/dyS/dyC; // TS
	Vloc[order[1]] =  -kW/dxW/dxC; // TW
	Vloc[order[2]] =  -(-kW/dxW/dxC -kE/dxE/dxC -kN/dyN/dyC -kS/dyS/dyC) + rho*Physics->Cp/dt; // TC
	Vloc[order[3]] =  -kE/dxE/dxC; // TE
	Vloc[order[4]] =  -kN/dyN/dyC; // TN


	*bloc = + rho*Physics->Cp*Physics->T[TC]/dt;

	// Add the contribution of the shear heating


	sigma_xy  = Physics->sigma_xy_0[ix-1 + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy-1)*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix   + (iy-1)*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy-1)*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix-1 + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix-1 + (iy  )*Grid->nxS];
	sigma_xy += Physics->sigma_xy_0[ix   + (iy  )*Grid->nxS] + Physics->Dsigma_xy_0[ix   + (iy  )*Grid->nxS];
	sigma_xy /= 4.0;

	sigma_xx = Physics->sigma_xx_0[ix+iy*nxEC] + Physics->Dsigma_xx_0[ix+iy*nxEC];

	eta = Physics->eta[TC];


	// Shear heating = Sxx*Exx + Syy*Eyy + Sxy*Exy + Syx*Eyx
	// with Sij, deviatoric stress tensor, Eij, deviatoric strain rate tensor
	// since Exy = Eyx and Sxx = -Syy (in 2D), then we can write:
	// H = 2*Sxx*Exx + 2*Sxy*Exy, or:
	//*bloc += sigma_xx*sigma_xx/eta + sigma_xy*sigma_xy/eta;
	//printf("bloc = %.2e, Hs = %.2e, eta = %.2e\n", *bloc, sigma_xx*sigma_xx/eta + sigma_xy*sigma_xy/eta, eta);
	//printf("Vloc[0] = %.2e, Vloc[1] = %.2e, Vloc[2] = %.2e, Vloc[3] = %.2e, Vloc[4] = %.2e\n",Vloc[0], Vloc[1], Vloc[2],Vloc[3], Vloc[4]);
	// Adiabatic heat should come here
	//*bloc += sigma_xx*sigma_xx/eta + sigma_xy*sigma_xy/eta;

}
#endif






#if (DARCY)


void LocalStencil_Stokes_Darcy_Momentum_x(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	// 1. call Stokes_Momentum_x
	LocalStencil_Stokes_Momentum_x(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);





	*nLoc = 13;



	int nxN = Grid->nxEC;
	int nxS = Grid->nxS;
	int PPeriod = 0;
	int nECTot = Grid->nECTot;
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y
	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;

	compute dxW, dxE, dxC;

	//compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	//compute dyN = Grid->DYEC[iy-1];;
	//compute dyC = 0.5*(dyS+dyN);
	if (Grid->isPeriodic) {
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

	if (Grid->isPeriodic) {

		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 1;
			}
			PPeriod  = 0;//nxN;


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
			order[11] = 12; // PcW
			order[12] = 11; // PcE

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

			order[11] = 11; // PcW
			order[12] = 12; // PcE

		}
	}

	Jloc[order[11]]   =   ix    + (iy)*nxN  + nVxTot+nVyTot+nECTot + PPeriod; // PcW
	Jloc[order[12]]   =   ix+1  + (iy)*nxN  + nVxTot+nVyTot+nECTot ; // PcE

	Vloc[order[11]] =  1.0/dxC;
	Vloc[order[12]] = -1.0/dxC;



	/*
	int NormalPeriod = nxN;
	int NormalW = ix-1+1    + (iy-1+1)*Grid->nxEC + NormalPeriod;
	int NormalE = ix  +1    + (iy-1+1)*Grid->nxEC;



	// BC for water and air
	// ==============================

	// Boundary condition for air
	if (Physics->phase[NormalW]==Physics->phaseAir || Physics->phase[NormalE]==Physics->phaseAir) {
		Vloc[order[ 9]] = 0.0; // PfW
		Vloc[order[10]] = 0.0; // PfE
		if (Physics->phase[NormalW]==Physics->phaseAir && Physics->phase[NormalE]==Physics->phaseAir) {
			*bloc -= 0.0;
		} else {
			*bloc -= - Physics->PfGrad_Air_X; // Value of the lateral gradient
		}
	}

	// BC for water
	if (Physics->phase[NormalW]==Physics->phaseWater) {

		*bloc -= Vloc[order[9]]*Physics->Plitho[NormalW]; // Dirichlet value
		Vloc[order[ 9]] = 0.0;
		// *bloc -= Vloc[order[11]]* 0.0; // Dirichlet value
		//Vloc[order[11]] = 0.0;
	}

	if (Physics->phase[NormalE]==Physics->phaseWater) {
		*bloc -= Vloc[order[10]]*Physics->Plitho[NormalE]; // Dirichlet value
		Vloc[order[10]] = 0.0;
		// *bloc -= Vloc[order[12]]* 0.0; // Dirichlet value
		//Vloc[order[12]] = 0.0;
	}

*/



}












void LocalStencil_Stokes_Darcy_Momentum_y(int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{

	// 1. call Stokes_Momentum_x
	LocalStencil_Stokes_Momentum_y(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);




	*nLoc = 13;



	int nxN = Grid->nxEC;
	//int nxS = Grid->nxS;
	int PPeriod = 0;
	int nECTot = Grid->nECTot;
	//int nxEC = Grid->nxEC;
	int nxVx = Grid->nxVx; // number of Vx nodes in x
	int nyVx = Grid->nyVx; // number of Vx nodes in y
	int nxVy = Grid->nxVy; // number of Vy nodes in x
	int nyVy = Grid->nyVy; // number of Vy nodes in y
	int nVxTot = nxVx*nyVx;
	int nVyTot = nxVy*nyVy;


	//compute dxW, dxE, dxC;

	compute dyS = Grid->DYS [iy-1];
	compute dyN = Grid->DYS [iy  ];
	compute dyC = 0.5*(dyS+dyN);
	//compute dt = Physics->dt;


/*
	if (Grid->isPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[0]+Grid->DXEC[nxEC-2]);
			dxE = Grid->DXEC[ix  ];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxVy-1) {
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
*/





	// 3. add contributions from Pc


	// Special cases for periodic BC
	if (Grid->isPeriodic) {
		if (ix==0) {
			PPeriod   = 0;//nxN    		; // PS
		}
		else if (ix==nxVy-3) {

		}
	}

	Jloc[order[11]] =   ix    + (iy  )*nxN + nVxTot+nVyTot+nECTot   + PPeriod     ; // PcS
	Jloc[order[12]] =   ix    + (iy+1)*nxN + nVxTot+nVyTot+nECTot   + PPeriod     ; // PcN

	Vloc[order[11]] =  1.0/dyC; // PcS
	Vloc[order[12]] = -1.0/dyC; // PcN


	/*

	// BC for water and air
	// ==============================
	int NormalN = ix-1+1    + (iy  +1)*Grid->nxEC ;
	int NormalS = ix-1+1    + (iy-1+1)*Grid->nxEC ;
	// Boundary condition for air
	if (Physics->phase[NormalS]==Physics->phaseAir || Physics->phase[NormalN]==Physics->phaseAir) {
		Vloc[order[ 9]] = 0.0; // PfS
		Vloc[order[10]] = 0.0; // PfN
		if (Physics->phase[NormalS]==Physics->phaseAir && Physics->phase[NormalN]==Physics->phaseAir) {
			*bloc -= - (Physics->Plitho[NormalN] - Physics->Plitho[NormalS])/dyC;
		} else {
			*bloc -= - Physics->PfGrad_Air_Y; // Value of the lateral gradient
		}
	}

	// BC for water
	if (Physics->phase[NormalS]==Physics->phaseWater) {
		*bloc -= Vloc[order[ 9]]*Physics->Plitho[NormalS]; // Dirichlet value
		Vloc[order[ 9]] = 0.0;
		// *bloc -= Vloc[order[11]]* 0.0; // Dirichlet value
		//Vloc[order[11]] = 0.0;
	}

	if (Physics->phase[NormalN]==Physics->phaseWater) {
		*bloc -= Vloc[order[10]] *Physics->Plitho[NormalN]; // Dirichlet value
		Vloc[order[10]] = 0.0;
		// *bloc -= Vloc[order[12]]* 0.0; // Dirichlet value
		//Vloc[order[12]] = 0.0;
	}
	*/



}





















void LocalStencil_Stokes_Darcy_Darcy 	 (int* order, int* Jloc, compute* Vloc, compute* bloc, int ix, int iy, Grid* Grid, Physics* Physics, int SetupType, int* shift, int* nLoc, int* IC)
{



	*bloc = 0; // just a security
	// 1. call Stokes_Momentum_x to build the velocity divergence
	LocalStencil_Stokes_Continuity(order, Jloc, Vloc, bloc, ix, iy, Grid, Physics, SetupType, shift, nLoc, IC);
	*bloc = 0; // just a security

	*nLoc = 9;
	*IC = 6;



	int nxEC = Grid->nxEC;
	//int nECTot = Grid->nECTot;
	int nVxTot = Grid->nVxTot;
	int nVyTot = Grid->nVyTot;


	compute dxW, dxE, dxC;

	compute dyS = Grid->DYEC[iy-1];//Grid->dy;//
	compute dyN = Grid->DYEC[iy-1];;
	compute dyC = 0.5*(dyS+dyN);


	if (Grid->isPeriodic) {
		if (ix==0) {
			dxW = 0.5*(Grid->DXEC[ix]+Grid->DXEC[nxEC-2]);
			dxE = Grid->DXEC[ix+1];
			dxC = 0.5*(dxW+dxE);

		} else if (ix==nxEC-1) {
			dxW = Grid->DXEC[ix-1];
			dxE = 0.5*(Grid->DXEC[ix]+Grid->DXS[nxEC-2]);
			dxC = 0.5*(dxW+dxE);

		} else {
			dxW = Grid->DXEC[ix-1];
			dxE = Grid->DXEC[ix];
			dxC = 0.5*(dxW+dxE);
		}

	} else {
		dxW = Grid->DXEC[ix-1];
		dxE = Grid->DXEC[ix];
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


	if (Grid->isPeriodic) {
		printf("error in LocalStencil_Stokes_Darcy_Darcy: periodic BC not implemented properly yet");
		exit(0);
		/*
		if (ix==0) {
			if (UPPER_TRI) {
				*shift = 5;
			}
			PPeriodL = 0;//nxC;
			PPeriodR = 0;
			order[4] = 4; // PfS
			order[5] = 7; // PfW
			order[6] = 5; // PfC
			order[7] = 6; // PfE
			order[8] = 8; // PfN
		}
		else if (ix == nxEC-1) {
			if (UPPER_TRI) {
				*shift = 7;
			}
			PPeriodL = 0;
			PPeriodR = 0;//-nxC;
			order[4] = 4; // PfS
			order[5] = 6; // PfW
			order[6] = 7; // PfC
			order[7] = 5; // PfE
			order[8] = 8; // PfN
		}
		 */
	}



	// the first 4 slots are filled by Vx and Vy equations
	Jloc[order[4]] = ix    + (iy-1)*nxEC + nVxTot+nVyTot; // PfS
	Jloc[order[5]] = ix-1  + (iy  )*nxEC + nVxTot+nVyTot + PPeriodL; // PfW
	Jloc[order[6]] = ix    + (iy  )*nxEC + nVxTot+nVyTot; // PfC
	Jloc[order[7]] = ix+1  + (iy  )*nxEC + nVxTot+nVyTot + PPeriodR; // PfE
	Jloc[order[8]] = ix    + (iy+1)*nxEC + nVxTot+nVyTot; // PfN





	int NormalS = ix   + (iy-1)*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalW = ix-1 + (iy  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalC = ix   + (iy  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalE = ix+1 + (iy  )*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells
	int NormalN = ix   + (iy+1)*nxEC; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells




		compute KS 		= ((Physics->perm_eta_f[NormalS] + Physics->perm_eta_f[NormalC])/2.0); // averaging because K has to be defined on the shear node
		compute KN 		= ((Physics->perm_eta_f[NormalN] + Physics->perm_eta_f[NormalC])/2.0); // averaging because K has to be defined on the shear node
		compute KW	 	= ((Physics->perm_eta_f[NormalW] + Physics->perm_eta_f[NormalC])/2.0); // averaging because K has to be defined on the shear node
		compute KE 		= ((Physics->perm_eta_f[NormalE] + Physics->perm_eta_f[NormalC])/2.0); // averaging because K has to be defined on the shear node



		Vloc[order[4]] = -( KS/dyS/dyC ); // PfS
		Vloc[order[5]] = -( KW/dxW/dxC ) ; // PfW
		Vloc[order[6]] = -( -KE/dxE/dxC -KW/dxW/dxC -KS/dyS/dyC -KN/dyN/dyC ); // PfC
		Vloc[order[7]] = -( KE/dxE/dxC ); // PfE
		Vloc[order[8]] = -( KN/dyN/dyC ); // PfN

		/*
		if (Physics->phase[NormalC] == Physics->phaseAir) {
			*bloc = 0;
		} else {
		*/
			*bloc -= KE*Physics->rho_f_g*Physics->gFac[0]/dxC - KW*Physics->rho_f_g*Physics->gFac[0]/dxC;
			*bloc -= KN*Physics->rho_f_g*Physics->gFac[1]/dyC - KS*Physics->rho_f_g*Physics->gFac[1]/dyC;
		//}


	//	printf("C = %i, KN = %.2e, permN = %.2e, Physics->perm[NormalC] = %.2e,  KN/dyN/dyC = %.2e\n", NormalC, KN, Physics->perm[NormalN], Physics->perm[NormalC], KN/dyN/dyC);
		/*
	printf("KW = %.2e, permW = %.2e, etaf = %.2e\n", KW, Physics->perm[NormalW], Physics->eta_f);
	printf("KE = %.2e, permE = %.2e, etaf = %.2e\n", KE, Physics->perm[NormalE], Physics->eta_f);
	printf("KN = %.2e, permN = %.2e, etaf = %.2e\n", KN, Physics->perm[NormalN], Physics->eta_f);
		 */
		//printf("dyS = %.2e, dyN = %.2e, dxW = %.2e, dxE = %.2e, dxC = %.2e, dyC = %.2e\n",dyS, dyN, dxW, dxE,dxC, dyC);





		//  BC for air and water


		//printf("KS = %.2e, KS/dyS/dyC=%.2e\n", KS, KS/dyS/dyC);
		//printf("Physics->phase[NormalC] = %i, Physics->phaseWater = %i, Physics->phaseAir = %i \n",Physics->phase[NormalC] , Physics->phaseWater, Physics->phaseAir );




/*
		if (Physics->phase[NormalC]==Physics->phaseAir) {
			Vloc[order[0]] = 0.0; // VxW
			Vloc[order[1]] = 0.0; // VxE
			Vloc[order[2]] = 0.0; // VxS
			Vloc[order[3]] = 0.0; // VxN
			Vloc[order[4]] = 0.0; // PfS
			Vloc[order[5]] = 0.0; // PfW
			Vloc[order[6]] = 1.0; // PfC
			Vloc[order[7]] = 0.0; // PfE
			Vloc[order[8]] = 0.0; // PfN

			*bloc = Physics->Plitho[NormalC]; // Dummy value(?), the correct value will be gotten back from the gradient
		//	printf("C = %i, is Air\n", NormalC);
		} else if (Physics->phase[NormalC]==Physics->phaseWater) {
			Vloc[order[0]] = 0.0; // VxW
			Vloc[order[1]] = 0.0; // VxE
			Vloc[order[2]] = 0.0; // VxS
			Vloc[order[3]] = 0.0; // VxN
			Vloc[order[4]] = 0.0; // PfS
			Vloc[order[5]] = 0.0; // PfW
			Vloc[order[6]] = 1.0; // PfC
			Vloc[order[7]] = 0.0; // PfE
			Vloc[order[8]] = 0.0; // PfN

			*bloc = Physics->Plitho[NormalC];
		//	printf("C = %i, is water\n", NormalC);
			//printf("Physics->Plitho[NormalC] = %.2e\n", Physics->Plitho[NormalC]);
		} else {


			// AIR
			// ===========================
			if (Physics->phase[NormalW]==Physics->phaseAir) {
				//printf("C = %i, W is Air\n", NormalC);
				Vloc[order[5]] = 0.0; // PfW
				Vloc[order[6]] -= -(-KW/dxW/dxC); // PfC
				*bloc -= - (- Physics->PfGrad_Air_X/dxC);
			}

			if (Physics->phase[NormalE]==Physics->phaseAir) {
			//	printf("C = %i, E is Air\n", NormalC);
				Vloc[order[7]] = 0.0; // PfE
				Vloc[order[6]] -= -(-KE/dxE/dxC); // PfC
				*bloc -= - ( Physics->PfGrad_Air_X/dxC);
			}
			if (Physics->phase[NormalS]==Physics->phaseAir) {
			//	printf("C = %i, S is Air\n", NormalC);
				Vloc[order[4]] = 0.0; // PfS
				Vloc[order[6]] -= -(-KS/dyS/dyC); // PfC
				*bloc -= - (- Physics->PfGrad_Air_Y/dyC);
			}
			if (Physics->phase[NormalN]==Physics->phaseAir) {
			//	printf("C = %i, N is Air\n", NormalC);
				Vloc[order[8]] = 0.0; // PfN
				Vloc[order[6]] -= -(-KN/dyN/dyC); // PfC
				*bloc -= - ( Physics->PfGrad_Air_Y/dyC);
			}



			// WATER
			// ===========================
			if (Physics->phase[NormalW]==Physics->phaseWater) {
			//	printf("C = %i, W is Water\n", NormalC);

				*bloc -= Vloc[order[5]]*Physics->Plitho[NormalW];
				Vloc[order[5]] = 0.0; // PfW
			}

			if (Physics->phase[NormalE]==Physics->phaseWater) {
			//	printf("C = %i, E is Water\n", NormalC);
				*bloc -= Vloc[order[7]]*Physics->Plitho[NormalE];
				Vloc[order[7]] = 0.0; // PfE
			}
			if (Physics->phase[NormalS]==Physics->phaseWater) {
			//	printf("C = %i, S is Water\n", NormalC);
				*bloc -= Vloc[order[4]]*Physics->Plitho[NormalS];
				Vloc[order[4]] = 0.0; // PfE
			}
			if (Physics->phase[NormalN]==Physics->phaseWater) {
			//	printf("C = %i, N is Water\n", NormalC);
				*bloc -= Vloc[order[8]]*Physics->Plitho[NormalN];
				Vloc[order[8]] = 0.0; // PfE
			}

		}



		*/

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

	int nxN = Grid->nxEC;
	int nVxTot = Grid->nVxTot;
	int nVyTot = Grid->nVyTot;
	int nECTot = Grid->nECTot;

	int NormalC = ix + (iy)*nxN; // +1 because Physics->B etc... are stored on embedded cells while P and Pf are on cells

	compute dt = Physics->dt;
	compute Zb;//, ZbStar;

	Jloc[order[4]] = ix    + (iy)*nxN + nVxTot+nVyTot+nECTot; // PcC

	//compute Eta_b =  Physics->eta_b[NormalC];
	compute eta_b =  Physics->eta[NormalC]/Physics->phi[NormalC];


	compute B = Physics->G[NormalC]/sqrt(Physics->phi[NormalC]);
	compute Khi_b = Physics->khi_b[NormalC];
	//printf("B = %.2e, phi = %.2e\n", B, Physics->phi[NormalC]);
	Zb = 1.0/(1.0/Khi_b + 1.0/eta_b + 1.0/(B*dt));//(dt*B) / (dt*B + eta_b);
	Zb *= (1.0 - Physics->phi[NormalC]);


	//ZbStar = (1.0 - Physics->phi[NormalC]) * Zb;

	Vloc[order[4]] =  1.0/ (Zb);
	//printf("Vloc = %.2e, eta_b = %.2e, Zbstar = %.2e, phi = %.2e, Zb = %.2e, B = %.2e, dt = %.2e\n", Vloc[order[4]], eta_b, ZbStar, Physics->phi[NormalC], Zb, B, dt);
	*bloc +=  Physics->DeltaP0[NormalC]  / (B*dt) ;

	//printf("Zb = %.2e, Zb* = %.2e, eta_b = %.2e, B = %.2e, phi = %.2e, bloc = %.2e\n", Zb, ZbStar, eta_b, B,  Physics->phi[NormalC], *bloc);

	/*
	if (Physics->phase[NormalC]==Physics->phaseAir || Physics->phase[NormalC]==Physics->phaseWater) {
		//Vloc[order[0]] = 0.0;
		//Vloc[order[1]] = 0.0;
		//Vloc[order[2]] = 0.0;
		//Vloc[order[3]] = 0.0;
		Vloc[order[4]] = 1.0;
		*bloc = 0.0;
	}
	*/



}


















#endif


































