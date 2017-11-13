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




#endif /* PARTICLES_H_ */

// code for higher order stress rotation
// after advection
/*
IX = round((thisParticle->x - Grid->xmin)/Grid->dx);
IY = round((thisParticle->y - Grid->ymin)/Grid->dy);

if (thisParticle->x<Grid->xmax && thisParticle->y<Grid->ymax && thisParticle->x>Grid->xmin && thisParticle->y>Grid->ymin) {
	locX = tempx-Grid->X[IX];
	locY = tempy-Grid->Y[IY];

	if (locX<0) {
		locX = 2.0*(locX/Grid->DXS[IX-1]);
	} else {
		locX = 2.0*(locX/Grid->DXS[IX]);
	}
	if (locY<0) {
		locY = 2.0*(locY/Grid->DYS[IY-1]);
	} else {
		locY = 2.0*(locY/Grid->DYS[IY]);
	}


}

compute alpha2 = Interp_NodeVal_Node2Particle_Local(alphaArray, ix, iy, Grid->nxS, Grid->nyS, locX, locY);
alpha = 0.5*(alpha+alpha2);
sigma_xx_temp = thisParticle->sigma_xx_0*cos(alpha)*cos(alpha) - thisParticle->sigma_xx_0*sin(alpha)*sin(alpha)  -  thisParticle->sigma_xy_0*sin(2.0*alpha);
thisParticle->sigma_xy_0 = thisParticle->sigma_xy_0*cos(2.0*alpha)  +  thisParticle->sigma_xx_0*sin(2.0*alpha);
thisParticle->sigma_xx_0 = sigma_xx_temp;
*/


// An older version of the advection routine. For reference
# if (0)

void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Declarations
	// =========================
	int iNode;
	compute locX, locY, locX0, locY0;
	int Ix, Iy;
	int ix, iy;
	int i;
	int ixN, iyN;
	compute alphaArray[4];
	compute alpha;
	compute sigma_xx_temp;

	SingleParticle* thisParticle;


	int iCell;
	int iVx, iVy;


	// Index of neighbouring cells, with respect to the node ix, iy
	compute xMod[4], yMod[4];
	xMod[0] = -1; yMod[0] = -1; // ll
	xMod[1] = -1; yMod[1] =  1; // ul
	xMod[2] =  1; yMod[2] =  1; // ur
	xMod[3] =  1; yMod[3] = -1; // lr

	int IxNV[4], IyNV[4];
	IxNV[0] =  0;   IyNV[0] =  1; // ul
	IxNV[1] =  1;	IyNV[1] =  1; // ur
	IxNV[2] =  0; 	IyNV[2] =  0; // ll
	IxNV[3] =  1; 	IyNV[3] =  0; // lr

	compute xModVx[4], yModVx[4], xModVy[4], yModVy[4];
	compute weight;
	compute Z;

	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  0;	IyN[1] =  1; // upper left
	IxN[2] =  1; 	IyN[2] =  1; // upper right
	IxN[3] =  1; 	IyN[3] =  0; // lower right
	int signX, signY;
	compute Vx, Vy;



	// Loop through inner cells
	// ========================
	iNode = 0;
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX0, locY0, locX, locY, signX, signY, i, ixN, iyN, alphaArray, alpha, sigma_xx_temp, Ix, Iy, Vx, Vy) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];




			// Loop through the particles in the cell
			// ======================================
			while (thisParticle!=NULL) {


				// Advect X
				// =====================
				locX0 = thisParticle->x-Grid->X[ix];
				locY0 = thisParticle->y-Grid->Y[iy];

				if (locX0<0) {
					locX0 = 2.0*(locX0/Grid->DXS[ix-1]);
				} else {
					locX0 = 2.0*(locX0/Grid->DXS[ix]);
				}
				if (locY0<0) {
					locY0 = 2.0*(locY0/Grid->DYS[iy-1]);
				} else {
					locY0 = 2.0*(locY0/Grid->DYS[iy]);
				}

				//locX0 = (thisParticle->x-Grid->xmin)/Grid->dx - ix;
				//locY0 = (thisParticle->y-Grid->ymin)/Grid->dy - iy;

				locX = locX0*1.0; // important for using shape functions
				locY = locY0*1.0;


				if (locX<0) {
					signX = -1;
				} else {
					signX = 1;
				}
				if (locY<0) {
					signY = -1;
				} else {
					signY = 1;
				}

				// add a condition with signX signY to avoid recomputing alpha if not necessary

				locX = fabs(locX)-1;
				locY = fabs(locY)-1;

				for (i=0;i<4;i++) {
					ixN = ix+IxN[i]*signX;
					iyN = iy+IyN[i]*signY;

					if (ixN+1>Grid->nxS || ixN<0 || iyN+1>Grid->nyS || iyN<0) {
						printf("error in Particles_advect: trying to access a non existing node\n");
						printf("IX = %i, IY = %i, locX = %.3f, locY = %.3f, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX);
						printf("thisParticle->x = %.3f , y = %.3f \n", thisParticle->x, thisParticle->y);
						exit(0);
					}

					alphaArray[i]  = - 0.5*Physics->dtAdv*((Physics->Vy[ixN+1+iyN*Grid->nxVy]   - Physics->Vy[ixN+(iyN)*Grid->nxVy])/Grid->DXEC[ix]
																																				- (Physics->Vx[ixN+(iyN+1)*Grid->nxVx] - Physics->Vx[ixN+(iyN)*Grid->nxVx])/Grid->DYEC[iy]);
					//printf("ix = %i, ixC = %i, iy = %i, iyC = %i, alphaArray[i] = %.3e\n", ix, ixC, iy, iyC, alphaArray[i]);
				}

				alpha =   ( .25*(1.0-locX)*(1.0-locY)*alphaArray[0]
																 + .25*(1.0-locX)*(1.0+locY)*alphaArray[1]
																										+ .25*(1.0+locX)*(1.0+locY)*alphaArray[2]
																																			   + .25*(1.0+locX)*(1.0-locY)*alphaArray[3] );

				// Jaumann co-rotation formulas (small angle approximation)
				//sigma_xx_corr = - thisParticle->sigma_xy_0 * 2 * alpha;
				//sigma_xy_corr = + thisParticle->sigma_xx_0 * 2 * alpha;
				//thisParticle->sigma_xx_0 += sigma_xx_corr;
				//thisParticle->sigma_xy_0 += sigma_xy_corr;

				// Correction without assuming a small angle
				sigma_xx_temp = thisParticle->sigma_xx_0*cos(1*alpha)*cos(1*alpha)  -  thisParticle->sigma_xy_0*sin(2*alpha);
				thisParticle->sigma_xy_0 = thisParticle->sigma_xy_0*cos(2*alpha)  +  thisParticle->sigma_xx_0*sin(2*alpha);
				thisParticle->sigma_xx_0 = sigma_xx_temp;

				//printf("alpha = %.3e, alphaArray[0] = %.3e, alphaArray[1] = %.3e, alphaArray[2] = %.3e, alphaArray[3] = %.3e\n", alpha, alphaArray[0], alphaArray[1], alphaArray[2], alphaArray[3]);






				locX = locX0*1.0; // important for using shape functions
				locY = locY0*1.0;

				if (locX>0.0) {
					locX = locX-1.0;
					Ix = ix;
					Iy = iy;
				}
				else {
					locX = locX+1.0;
					Ix = ix-1;
					Iy = iy;
				}



				Vx = ( .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
															 + .25*(1.0-locX)*(1.0+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
																									 + .25*(1.0+locX)*(1.0+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
																																			 + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) ;


				// Advect Y
				// =====================


				locX = locX0*1.0; // important for using shape functions
				locY = locY0*1.0;


				if (locY>0.0) {
					locY = locY-1.0;
					Ix = ix;
					Iy = iy;
				}
				else {
					locY = locY+1.0;
					Ix = ix;
					Iy = iy-1;
				}
				//printf("iP=%i, Ix=%i, Iy=%i, locX=%.2f, locY=%.2f w0=%.3f, w1=%.3f, w2=%.3f, w3=%.3f \n",iP, Ix, Iy, locX, locY, .25*(1.0-locX)*(1.0-locY), .25*(1.0-locX)*(1.0+locY), .25*(1.0+locX)*(1.0+locY), .25*(1.0+locX)*(1.0-locY));

				Vy  = (.25*(1.0-locX)*(1.0-locY)*Physics->Vy[Ix  +(Iy  )*Grid->nxVy]
															 + .25*(1.0-locX)*(1.0+locY)*Physics->Vy[Ix  +(Iy+1)*Grid->nxVy]
																									 + .25*(1.0+locX)*(1.0+locY)*Physics->Vy[Ix+1+(Iy+1)*Grid->nxVy]
																																			 + .25*(1.0+locX)*(1.0-locY)*Physics->Vy[Ix+1+(Iy  )*Grid->nxVy] ) ;






#if (ADVECT_VEL_AND_VISCOSITY)
				// get etaVisc on this particle

				locX = locX0*1.0; // important for using shape functions
				locY = locY0*1.0;

				Z  = ( .25*(1.0-locX)*(1.0-locY)*Physics->Z[ix  +(iy  )*Grid->nxEC]
															+ .25*(1.0-locX)*(1.0+locY)*Physics->Z[ix  +(iy+1)*Grid->nxEC]
																								   + .25*(1.0+locX)*(1.0+locY)*Physics->Z[ix+1+(iy+1)*Grid->nxEC]
																																		  + .25*(1.0+locX)*(1.0-locY)*Physics->Z[ix+1+(iy  )*Grid->nxEC] );

#if (DARCY)
				Zb  = ( .25*(1.0-locX)*(1.0-locY)*Physics->Zb[ix  +(iy  )*Grid->nxEC]
															  + .25*(1.0-locX)*(1.0+locY)*Physics->Zb[ix  +(iy+1)*Grid->nxEC]
																									  + .25*(1.0+locX)*(1.0+locY)*Physics->Zb[ix+1+(iy+1)*Grid->nxEC]
																																			  + .25*(1.0+locX)*(1.0-locY)*Physics->Zb[ix+1+(iy  )*Grid->nxEC] );
#endif

#endif










				// Advect particles
				thisParticle->x += Vx* Physics->dtAdv;
				thisParticle->y += Vy* Physics->dtAdv;







#if (ADVECT_VEL_AND_VISCOSITY)

				//printf("ix = %i, iy = %i\n", ix, iy);

				// Interpolate velocities from particles back to nodes
				locX = locX0; // important for using shape functions
				locY = locY0;

				compute xModVx[4], yModVx[4], xModVy[4], yModVy[4];
				xModVx[0] = -1.0; // ul
				xModVx[1] =  1.0; // ur
				xModVx[2] = -1.0; // ll
				xModVx[3] =  1.0; // lr

				yModVx[0] =  1.0; // ul
				yModVx[1] =  1.0; // ur
				yModVx[2] = -1.0; // ll
				yModVx[3] = -1.0; //lr




				xModVy[0] = -1.0; // ul
				xModVy[1] =  1.0; // ur
				xModVy[2] = -1.0; // ll
				xModVy[3] =  1.0; // lr

				yModVy[0] = -1.0; // ul
				yModVy[1] = -1.0; // ur
				yModVy[2] =  1.0; // ll
				yModVy[3] =  1.0; //lr

				int ixMod;
				int iyMod;
				if (locX>=0) {
					signX = 1.0;
					ixMod = 0;
				} else {
					signX =-1.0;
					ixMod = -1;
				}
				if (locY>=0) {
					signY = 1.0;
					iyMod = 0;
				} else {
					signY =-1.0;
					iyMod = -1;
				}


				for (i=0; i<4; i++) {
					iVx = (ix+IxNV[i]+ixMod + (iy+IyNV[i]) * Grid->nxVx);
					//weight = fabs((locX + xModVx[i]*0.5)   *   (locY + yModVx[i]*0.5));

					//compute A = (0.5+signX*xModVx[i]*(fabs(locX)-0.5) );
					//compute B = fabs(locY + yModVx[i]*0.5);
					weight = (0.5+signX*xModVx[i]*(fabs(locX)-0.5) )*fabs(locY + yModVx[i]*0.5);
					//printf("locX = %.2f, locY = %.2f, ix = %i, ixN = %i, iy = %i, iyN = %i, weight = %.2f, xContrib = %.2f, yContrib = %.2f\n", locX, locY, ix, ix+IxNV[i]+ixMod, iy, (iy+IyNV[i]) , weight, A, B);
					VxGrid[iVx*4+i] += Vx * weight;
					sumOfWeights_Vx[iVx*4+i] += weight;

				}



				//	printf("Vy\n");

				for (i=0; i<4; i++) {
					iVy = (ix+IxNV[i] + (iy+IyNV[i]+iyMod) * Grid->nxVy);

					//compute A = fabs(locX + xModVy[i]*0.5);
					//compute B = (0.5+signY*yModVy[i]*(fabs(locY)-0.5));
					weight =    fabs(locX + xModVy[i]*0.5) * (0.5+signY*yModVy[i]*(fabs(locY)-0.5)) ;

					//printf("locX = %.2f, locY = %.2f, ix = %i, ixN = %i, iy = %i, iyN = %i, weight = %.2f, xContrib = %.2f, yContrib = %.2f\n", locX, locY, ix, ix+IxNV[i], iy, (iy+IyNV[i]+iyMod) , weight, A, B);


					VyGrid[iVy*4+i] += Vy * weight;
					sumOfWeights_Vy[iVy*4+i] += weight;

				}



				//exit(0);





				// Interpolate etaVisc back to Nodes
				locX = locX0; // important for using shape functions
				locY = locY0;


				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);

					weight = fabs((locX + xMod[i]*0.5)   *   (locY + yMod[i]*0.5));

					sumOfWeights_EC[4*iCell+i] += weight;
					ZGrid[4*iCell+i] += Z*weight;
#if (DARCY)
					ZbGrid[4*iCell+i] += Zb*weight;
#endif
				}



#endif












				thisParticle = thisParticle->next;
			}
		}
	}

	//printf("Z\n");
	compute sum = 0;


#if (ADVECT_VEL_AND_VISCOSITY)
	for (iVx = 0; iVx < Grid->nVxTot; ++iVx) {
		sum = sumOfWeights_Vx[4*iVx+0] + sumOfWeights_Vx[4*iVx+1] + sumOfWeights_Vx[4*iVx+2] + sumOfWeights_Vx[4*iVx+3];
		Physics->Vx[iVx] = ( VxGrid[4*iVx+0] + VxGrid[4*iVx+1] + VxGrid[4*iVx+2] + VxGrid[4*iVx+3] ) /sum;
	}

	for (iVy = 0; iVy < Grid->nVyTot; ++iVy) {
		sum = sumOfWeights_Vy[4*iVy+0] + sumOfWeights_Vy[4*iVy+1] + sumOfWeights_Vy[4*iVy+2] + sumOfWeights_Vy[4*iVy+3];
		Physics->Vy[iVy] = ( VyGrid[4*iVy+0] + VyGrid[4*iVy+1] + VyGrid[4*iVy+2] + VyGrid[4*iVy+3] ) /sum;
	}


	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		sum = sumOfWeights_EC[4*iCell+0] + sumOfWeights_EC[4*iCell+1] + sumOfWeights_EC[4*iCell+2] + sumOfWeights_EC[4*iCell+3];
		Physics->Z[iCell] = ( ZGrid[4*iCell+0] + ZGrid[4*iCell+1] + ZGrid[4*iCell+2] + ZGrid[4*iCell+3] ) / sum;
#if (DARCY)
		Physics->Zb[iCell] = ( ZbGrid[4*iCell+0] + ZbGrid[4*iCell+1] + ZbGrid[4*iCell+2] + ZbGrid[4*iCell+3] ) / sum;
#endif
	}





	free(VxGrid);
	free(VyGrid);

	free(sumOfWeights_Vx);
	free(sumOfWeights_Vy);

	free(sumOfWeights_EC);
	free(ZGrid);
#if (DARCY)
	free(ZbGrid);
#endif

	for (iy = 0; iy<Grid->nyS; iy++) {
		for (ix = 0; ix<Grid->nxS; ix++) {
			Physics->ZShear[ix + iy*Grid->nxS] = shearValue(Physics->Z,  ix   , iy, Grid->nxEC);
		}
	}
#endif
	//printf("out of Advect\n");

}

#endif // ADVECT_METHOD == 0
