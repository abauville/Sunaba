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
			Particles->oldCellId[C] = i;
			Particles->newCellId[C] = i;

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

	// Update the link list
	// =========================
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		iP = Particles->linkHead[iCell];

		//printf("iCell = %i\n", iCell);


		// Follow the links through the cell (i.e. until next==-1)
		while (iP != -1) {

			// If this particle has changed cell
			if (Particles->oldCellId[iP] != Particles->newCellId[iP]) {

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

		Particles->linkNext[iP] = Particles->linkHead[Particles->newCellId[iP]] ;
		Particles->linkHead[Particles->newCellId[iP]] = iP;
	}

	freeLinkedList(headIdChanged);



	if (DEBUG) {
		// Check implementation
		// ====================
		for (iP = 0; iP < Particles->n; ++iP) {
			printf("oCellId = %i, nCellId = %i, iP = %i, Next = %i\n",Particles->oldCellId[iP], Particles->newCellId[iP], iP,  Particles->linkNext[iP]);
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























