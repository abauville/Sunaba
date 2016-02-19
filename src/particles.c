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
			x = ix*Grid->dx;
			y = iy*Grid->dy;



			// Loop through Particles in the cell
			// ==================================
			for (iPy=0;iPy<Particles->nPCY;iPy++) {
				for (iPx=0;iPx<Particles->nPCX;iPx++) {

					// Assign coordinate
					// =================
					Particles->xy[I] 	= x + 0.5*dxP + iPx*dxP + dxP*(0.5 - (rand() % 1000)/1000.0);
					Particles->xy[I+1] 	= y + 0.5*dyP + iPy*dyP + dyP*(0.5 - (rand() % 1000)/1000.0);
					I += 2;

				} // iPx
			} // iPy



		} // ix
	} // iy




	// Init CellId and link list
	// =========================
	C = 0;

	// Loop through cells
	for (i = 0; i < Grid.nCTot; ++i) {
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

	for (i = 0; i < Particles->n; ++i) {
		sqrDistance = Particles->xy[2*i]*Particles->xy[2*i]  +  Particles->xy[2*i+1]*Particles->xy[2*i+1]; // d^2 = x^2 + y^2
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
//                                INIT PHASE                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
void Particles_updateLinkedList(Grid* Grid, Particles* Particles)
{

	// Dummy change in the newCellId, for testing only
	// =========================
	Particles->newCellId[4] = 2;

	// Declarations
	// =========================
	int iCell;
	int iP; // particule index
	int it;


	// Update the link list
	// =========================
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		iP = Particles->linkHead[iCell];



		// Follow the links through the cell (i.e. until next==-1)
		while (Particles->linkNext[iP] != -1) {

			// If this particle has changed cell
			if (Particles->oldCellId[iP] != Particles->newCellId[iP]) {

				if (Particles->linkNext[iP]==-1) {
					Particles->linkNext[iP] = Particles->linkNext[ Particles->linkNext[iP] ];
				}

				//if (iP == Particles->linkHead[iCell]) {

				//}


			}
			else {

			}

			iP = Particles->linkNext[iP];
		}


	}







}


void getPhysicsFromParticles2Grid(Grid* Grid, Particles* Particles, Physics* Physics, MatProps* MatProps)
{

	// Update the link list

}



















