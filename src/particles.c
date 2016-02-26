/*
 * particles.c
 *
 *  Created on: Feb 18, 2016
 *      Author: abauville
 */

#include "stokes.h"


//============================================================================//
//============================================================================//
//                                                                            //
//                               INIT COORD                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
void Particles_initCoord(Grid* Grid, Particles* Particles)
{
	// Declare variables
	// ==================
	int i, j, ix, iy, I, iPx, iPy, C;
	coord x, y;
	coord dxP = Grid->dx/Particles->nPCX;
	coord dyP = Grid->dy/Particles->nPCY;
	I = 0;

	// Init random number generator
	// ==================
	srand(time(NULL));







	// Loop through cells
	// ==================
	for(iy=0;iy<Grid->nyC;iy++) {
		for(ix=0;ix<Grid->nxC;ix++) {
			// Get the coordinates of the lower left corner of the cell
			x = Grid->xmin + ix*Grid->dx;
			y = Grid->ymin + iy*Grid->dy;



			// Loop through Particles in the cell
			// ==================================
			for (iPy=0;iPy<Particles->nPCY;iPy++) {
				for (iPx=0;iPx<Particles->nPCX;iPx++) {

					// Assign coordinate
					// =================
					Particles->xy[I] 	= x + 0.5*dxP + iPx*dxP;// + dxP*(0.5 - (rand() % 1000)/1000.0);
					Particles->xy[I+1] 	= y + 0.5*dyP + iPy*dyP;// + dyP*(0.5 - (rand() % 1000)/1000.0);
					I += 2;

				} // iPx
			} // iPy



		} // ix
	} // iy




	// Init CellId and link list
	// =========================
	C = 0;

	// Loop through cells
	for (i = 0; i < Grid->nCTot; ++i) {
		Particles->linkHead[i] = C;

		// Loop through particles in the cell
		for (j = 0; j < Particles->nPC; ++j) {
			Particles->cellId[C] = i;

			// Fill linkNext
			if (j<Particles->nPC-1) {
				Particles->linkNext[C] = C+1;
			}
			else { // if last particle of the cell, thebn linkNext = -1 (i.e. NULL)
				Particles->linkNext[C] = -1;
			}

			C++;
		}

	}

	if (DEBUG) {
		printf("Linked List\n");
		int iCell, iP;
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("Cell #%i:  ", iCell);
			iP = Particles->linkHead[iCell];
			while (iP!=-1) {
				printf("%i  ",iP);
				iP = Particles->linkNext[iP];
			}
			printf("\n");
		}
	}
}








//============================================================================//
//============================================================================//
//                                                                            //
//                                INIT PHASE                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
void Particles_initPhase(Grid* Grid, Particles* Particles)
{
	// Simple inclusion
	int i;
	coord sqrDistance;
	coord sqrRadius = 0.3*0.3;
	coord cX = 0;
	coord cY = 0;

	for (i = 0; i < Particles->n; ++i) {
		sqrDistance = (Particles->xy[2*i  ]-cX)*(Particles->xy[2*i  ]-cX)
				    		+ (Particles->xy[2*i+1]-cY)*(Particles->xy[2*i+1]-cY); // d^2 = x^2 + y^2
		//printf("i = %i, x = %.2f, y = %.2f, sqrDistance = %.2f\n",i,sqrDistance, Particles->xy[2*i  ], Particles->xy[2*i+1]);
		if (sqrDistance < sqrRadius) {
			Particles->phase[i] = 1;
		}
		else {
			Particles->phase[i] = 0;
		}
	}
}





//============================================================================//
//============================================================================//
//                                                                            //
//                      UPDATE THE PARTICLE LINKED LIST                  	  //
//                                                                            //
//============================================================================//
//============================================================================//
struct IdChanged {
	int id;
	struct IdChanged*  next;
};


void Particles_updateLinkedList(Grid* Grid, Particles* Particles)
{
	printf("Begin Update Linked List\n");

	// Dummy change in the newCellId, for testing only
	// =========================
	//Particles->newCellId[4] = 2;
	//Particles->newCellId[7] = 2;
	//Particles->newCellId[9] = 4;
	//Particles->newCellId[15] = 0;
	// Declarations
	// =========================
	int iCell;
	int iP, iPPrevious; // particle index

	// Declare a linked list that contains the id of particles that have change cell
	LinkedNode* headIdChanged = (LinkedNode*) malloc(sizeof(LinkedNode));
	headIdChanged->data = 0;
	headIdChanged->next = NULL;

	int oldCellId;
	coord x, y;
	int ix, iy;
	// Update the link list
	// =========================
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		iP = Particles->linkHead[iCell];

		//printf("iCell = %i\n", iCell);


		// Follow the links through the cell (i.e. until next==-1)
		while (iP != -1) {
			oldCellId = Particles->cellId[iP];

			x = Particles->xy[2*iP];
			y = Particles->xy[2*iP+1];
			ix = (int) floor(x/Grid->dx);
			iy = (int) floor(y/Grid->dx);
			printf("iP=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,x, y, ix,iy);
			Particles->cellId[iP] = ix + iy*Grid->nxC;

			// If this particle has changed cell
			if (oldCellId != Particles->cellId[iP]) {

				// 1. Update info for the oldCell
				// ===========================
				if (iP != Particles->linkHead[iCell]) {
					Particles->linkNext[ iPPrevious ] = Particles->linkNext[iP];
				}
				else {
					Particles->linkHead[iCell] = Particles->linkNext[iP];
				}


				addToLinkedList(&headIdChanged, iP);


			}
			iPPrevious = iP;
			iP = Particles->linkNext[iP];
		}
	}

	//printf("End loop\n");
	// 2. Update info of the new cell, i.e. Add this particle to the head of the link list of the new cell
	// ===============================
	LinkedNode* IdChanged = NULL;
	IdChanged = headIdChanged;
	while (IdChanged->next!=NULL) {
		iP 			= IdChanged->data;
		IdChanged 	= IdChanged->next;

		Particles->linkNext[iP] = Particles->linkHead[Particles->cellId[iP]] ;
		Particles->linkHead[Particles->cellId[iP]] = iP;
	}

	freeLinkedList(headIdChanged);



	if (DEBUG) {
		// Check implementation
		// ====================
		for (iP = 0; iP < Particles->n; ++iP) {
			printf("cellId = %i, iP = %i, Next = %i\n",Particles->cellId[iP], iP,  Particles->linkNext[iP]);
		}
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("%i  ", Particles->linkHead[iCell]);
		}
		printf("\n\n\n");


		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("Cell #%i:  ", iCell);
			iP = Particles->linkHead[iCell];
			while (iP!=-1) {
				printf("%i  ",iP);
				iP = Particles->linkNext[iP];
			}
			printf("\n");
		}
	}
}


void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Declarations
	// =========================
	int iCell, iP;
	compute locX, locY;
	int Ix, Iy;
	int ix, iy;

	compute Vtemp;
	// Loop through inner cells
	// ========================
	iCell = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			iCell = ix  + (iy  )*Grid->nxC;
			iP = Particles->linkHead[iCell];

			// Loop through the particles in the cell
			// ======================================
			while (iP!=-1) {
				// Advect X
				// =====================
				locX = (Particles->xy[2*iP  ]-Grid->xmin)/Grid->dx - ix;
				locY = (Particles->xy[2*iP+1]-Grid->ymin)/Grid->dy - iy;

				locX = locX*2-1.0; // important for using shape functions
				locY = locY*2-1.0;


				if (locY>0.0) {
					locY = locY-1.0;
					Ix = ix;
					Iy = iy+1;
				}
				else {
					Ix = ix;
					Iy = iy;
				}
				Vtemp = ( .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
				 	    + .25*(1.0-locX)*(1.0+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
					    + .25*(1.0+locX)*(1.0+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
				 	    + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) ;
				printf("iP=%i, Vtemp=%.2f, Vx0=%.2f, Vx1=%.2f, Vx2=%.2f, Vx3=%.2f, Coeff=%.2f, Coeff2=%.2f\n", iP, Vtemp, Physics->Vx[Ix  +(Iy  )*Grid->nxVx],
						Physics->Vx[Ix  +(Iy+1)*Grid->nxVx],Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx],Physics->Vx[Ix+1+(Iy  )*Grid->nxVx], .25*(1.0-locX)*(1.0-locY),  .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]);
				Particles->xy[iP*2] += ( 1/4*(1-locX)*(1-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
									   + 1/4*(1-locX)*(1+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
									   + 1/4*(1+locX)*(1+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
					                   + 1/4*(1+locX)*(1-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) * Physics->dt;


				// Advect Y
				// =====================
				locX = (Particles->xy[2*iP  ]-Grid->xmin)/Grid->dx - ix;
				locY = (Particles->xy[2*iP+1]-Grid->ymin)/Grid->dy - iy;

				locX = locX*2-1.0; // important for using shape functions
				locY = locY*2-1.0;


				if (locX>0.0) {
					locX = locX-1.0;
					Ix = ix+1;
					Iy = iy;
				}
				else {
					Ix = ix;
					Iy = iy;
				}
				Particles->xy[iP*2+1] += ( 1/4*(1-locX)*(1-locY)*Physics->Vy[Ix  +(Iy  )*Grid->nxVy]
										 + 1/4*(1-locX)*(1+locY)*Physics->Vy[Ix  +(Iy+1)*Grid->nxVy]
									     + 1/4*(1+locX)*(1+locY)*Physics->Vy[Ix+1+(Iy+1)*Grid->nxVy]
					                     + 1/4*(1+locX)*(1-locY)*Physics->Vy[Ix+1+(Iy  )*Grid->nxVy] ) * Physics->dt;


				iP = Particles->linkNext[iP];
			}
		}
	}

}




















