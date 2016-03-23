/*
 * particles.c
 *
 *  Created on: Feb 18, 2016
 *      Author: abauville
 */

#include "stokes.h"

// Example of sweeping through the Particles:
/*
SingleParticle* thisParticle = NULL;
int iCell;
for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
	thisParticle = Particles->linkHead[iCell];
	while (thisParticle != NULL) {


		thisParticle = thisParticle->next;
	}
}
 */



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
	int ix, iy, iPx, iPy;

	coord x, y;
	coord dxP = Grid->dx/Particles->nPCX;
	coord dyP = Grid->dy/Particles->nPCY;


	// Init random number generator
	// ==================
	srand(time(NULL));




	coord xP, yP;


	// Loop through cells
	// ==================
	int iCell = 0;
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
					//printf("Rand1 = %.4f, Rand2 = %.4f\n",(0.5 - (rand() % 1000)/1000.0),(0.5 - (rand() % 1000)/1000.0));
					xP 	= x + 0.5*dxP + iPx*dxP + Particles->noiseFactor*dxP*(0.5 - (rand() % 1000)/1000.0);
					yP 	= y + 0.5*dyP + iPy*dyP + Particles->noiseFactor*dyP*(0.5 - (rand() % 1000)/1000.0);

					// Create a particle
					addSingleParticle(&Particles->linkHead[iCell], xP, yP, 0, iCell);


					//	C++;



				} // iPx
			} // iPy


			// Fill the Head array




			iCell++;
		} // ix
	} // iy



	/*
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
	 */

	/*
	if (DEBUG) {
		printf("Linked List\n");
		int iCell, iP;
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("Cell #%i:  ", iCell);
			iP = Particles->linkHead[iCell];
			while (thisParticle!=NULL) {
				printf("%i  ",iP);
				iP = Particles->linkNext[iP];
			}
			printf("\n");
		}
	}
	 */
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
	int Setup = 2;
	srand(time(NULL));

	if (Setup==0) {

		// Simple inclusion
		int object = 1; // 0 = circle, 1 = square

		coord sqrDistance;
		coord radius = (0.3*(Grid->ymax-Grid->ymin)/2);
		coord sqrRadius =  radius * radius;
		//coord sqrRadius = 0.3*0.3;
		coord cX = 0;
		coord cY = Grid->ymin + (Grid->ymax-Grid->ymin)*0.8;//Grid->ymin + 0.0*(Grid->ymax-Grid->ymin)/2.0;



		FOR_PARTICLES


		if (object == 0) {
			sqrDistance = (thisParticle->x-cX)*(thisParticle->x-cX) + (thisParticle->y-cY)*(thisParticle->y-cY);
			if (sqrDistance < sqrRadius) {
				thisParticle->phase = 1;
			}
			else {
				thisParticle->phase = 0;
			}
		}
		else if (object == 1) {
			if (abs(thisParticle->x-cX) < radius && abs(thisParticle->y-cY) <radius) {
				thisParticle->phase = 1;
			}
			else {
				thisParticle->phase = 0;
			}
		}

		END_PARTICLES




	}


	else if (Setup==1) {
		// Sinusoidal basement
		compute WaveNumber = 3; // Wavelength
		compute phase = 0.5*PI;
		compute Amplitude = 0.1*(Grid->ymax-Grid->ymin);
		compute Thickness = 0.3*(Grid->ymax-Grid->ymin);
		compute x,y;



		FOR_PARTICLES
		x = thisParticle->x;
		y = thisParticle->y;

		x = (x-Grid->xmin)/(Grid->xmax-Grid->xmin); // x is now between 0 and 1
		x = x*WaveNumber*2.0*PI;
		if (y<(Amplitude*sin(x+phase) + Grid->ymin+Thickness)) {
			thisParticle->phase = 1;
		}
		else {
			thisParticle->phase = 0;
		}
		END_PARTICLES

	}
	else if (Setup==2) {
		// MultiLayer
		int iL;
		int nLayers = 8; // Wavelength
		compute spacing = 0.015*(Grid->ymax-Grid->ymin);

		compute spaceBelow = 0.4*(Grid->ymax-Grid->ymin);

		compute Thickness = 0.02*(Grid->ymax-Grid->ymin);
		compute thicknessNoiseFactor = 0.5;
		compute spacingNoiseFactor = 0.5; // between 0 and 1
		compute layerNoiseFactor = 0.0;
		compute noise;

		// sinusoidal perturbation
		compute WaveNumber = 1; // Wavelength
		compute phase = 0.5*PI;
		compute Amplitude = 0.0*(Grid->ymax-Grid->ymin);


		compute x,y;
		compute yLayTop, yLayBot;

		compute yBot = Grid->ymin + spaceBelow;
		compute yTop;


		FOR_PARTICLES
		thisParticle->phase = 0;
		END_PARTICLES


		for (iL = 0; iL < nLayers; ++iL) {
			yTop = yBot+Thickness+   Thickness*thicknessNoiseFactor*(0.5 - (rand() % 1000)/1000.0);
			FOR_PARTICLES
			x = thisParticle->x;
			y = thisParticle->y;

			x = (x-Grid->xmin)/(Grid->xmax-Grid->xmin); // x is now between 0 and 1
			x = x*WaveNumber*2.0*PI;

			yLayTop = (Amplitude*sin(x+phase) + yTop);
			yLayBot = (Amplitude*sin(x+phase) + yBot);
			noise = Thickness*layerNoiseFactor*(0.5 - (rand() % 1000)/1000.0);
			if (y>yLayBot+noise &&  y<yLayTop+noise) {
				thisParticle->phase = 1;
			}

			END_PARTICLES
			yBot = yTop + spacing +   spacing*spacingNoiseFactor*(0.5 - (rand() % 1000)/1000.0);
		}



	}
	else {
		printf("Unknwon Setup %i in Particles_initPhase\n", Setup);
		exit(0);



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

	// Dummy change in the newCellId, for testing only
	// =========================
	//Particles->newCellId[4] = 2;
	//Particles->newCellId[7] = 2;
	//Particles->newCellId[9] = 4;
	//Particles->newCellId[15] = 0;
	// Declarations
	// =========================
	//int iCell;

	// Declare a linked list that contains the id of particles that have change cell
	ParticlePointerList* headIdChanged = (ParticlePointerList*) malloc(sizeof(ParticlePointerList*));
	headIdChanged->pointer = NULL;
	headIdChanged->next = NULL;

	int oldCellId;
	coord x, y;
	int ix, iy;
	SingleParticle* previousParticle;
	// Update the link list
	// =========================


	int ParticleCounter = 0;
	SingleParticle* thisParticle = NULL;
	int iCell = 0;
	int phase;
	int TotNumParticles = 0;
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		thisParticle = Particles->linkHead[iCell];
		ParticleCounter = 0;
		while (thisParticle != NULL) {

			ParticleCounter++;


			oldCellId = thisParticle->cellId;

			x = thisParticle->x;
			y = thisParticle->y;
			ix = (int) floor((x-Grid->xmin)/Grid->dx);
			iy = (int) floor((y-Grid->ymin)/Grid->dy);

			//printf("x = %.1f , Grid->xmin = %.1f", Particles->xy[2*iP],Grid->xmin );

			thisParticle->cellId = ix + iy*Grid->nxC;
			//printf("iP=%i, oid=%i, nid=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,oldCellId, Particles->cellId[iP],x, y, ix,iy);
			// If this particle has changed cell
			if (oldCellId != thisParticle->cellId) {
				//printf("iP=%i, oid=%i, nid=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,oldCellId, Particles->cellId[iP],x, y, ix,iy);
				// 1. Update info for the oldCell
				// ===========================
				if (thisParticle != Particles->linkHead[iCell]) {
					previousParticle->next = thisParticle->next;
				}
				else {
					Particles->linkHead[iCell] = thisParticle->next;
				}


				addToParticlePointerList(&headIdChanged, thisParticle);


			}
			else {
				previousParticle = thisParticle;
			}


			thisParticle = thisParticle->next;
		}


	}

	/*
	printf("\n\n\n\n ========== Second part\n");
	for (iP=0;iP<Particles->n;iP++) {
		printf("cellId[%i] = %i\n",iP,Particles->cellId[iP]);
	}
	 */

	//printf("End loop\n");
	// 2. Update info of the new cell, i.e. Add this particle to the head of the link list of the new cell
	// ===============================
	ParticlePointerList* IdChanged = NULL;
	IdChanged = headIdChanged;
	while (IdChanged->next!=NULL) {
		thisParticle 	= IdChanged->pointer;
		IdChanged 		= IdChanged->next;
		thisParticle->next = Particles->linkHead[thisParticle->cellId] ;
		Particles->linkHead[thisParticle->cellId] = thisParticle;
	}

	freeParticlePointerList(headIdChanged);



	// Extra sweep to inject or delete particle
	// note: Not optimal this could be done during another sweep, for example during interpolation
	TotNumParticles = 0;
	for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
		thisParticle = Particles->linkHead[iCell];
		ParticleCounter=0;
		while (thisParticle != NULL) {
			ParticleCounter++;
			TotNumParticles++;

			thisParticle = thisParticle->next;
		}

		if (ParticleCounter==0) {
			printf("Warning: cell #%i is empty\n", iCell);
		}

		else if (ParticleCounter<Particles->nPC/2) { // integer division, should be rounded properly
			thisParticle = Particles->linkHead[iCell];


			ix = iCell%Grid->nxC;
			iy = (iCell-ix)/Grid->nxC;
			x = Grid->xmin + ix*Grid->dx + 0.5*Grid->dx ;
			y = Grid->ymin + iy*Grid->dy + 0.5*Grid->dy ;
			phase = thisParticle->phase; // the phase given to the particles is the phase of the head particle. Easy and fast but not optimal

			addSingleParticle(&Particles->linkHead[iCell], x, y, phase, iCell);
			Particles->n+=1;

		}

		else if (ParticleCounter>Particles->nPC*3) { // integer division, should be rounded properly
			thisParticle = Particles->linkHead[iCell];
			Particles->linkHead[iCell] = Particles->linkHead[iCell]->next;
			free(thisParticle);
			Particles->n-=1;

		}





	}
	//printf("\n");
	printf("TotNumParticles = %i\n", TotNumParticles);


	if (DEBUG) {
		// Check implementation
		// ====================
		/*
		for (iP = 0; iP < Particles->n; ++iP) {
			printf("cellId = %i, iP = %i, Next = %i\n",Particles->cellId[iP], iP,  Particles->linkNext[iP]);
		}
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("%i  ", Particles->linkHead[iCell]);
		}
		printf("\n\n\n");
		 */

		/*
		for (iCell = 0; iCell < Grid->nCTot; ++iCell) {
			printf("Cell #%i:  ", iCell);
			thisParticle = Particles->linkHead[iCell];
			while (thisParticle!=NULL) {
				printf("%i  ",iP);
				thisParticle = thisParticle->next;
			}
			printf("\n");
		}
		 */
	}
}


void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Declarations
	// =========================
	int iCell;
	compute locX, locY, locX0, locY0;
	int Ix, Iy;
	int ix, iy;

	SingleParticle* thisParticle;

	// Loop through inner cells
	// ========================
	iCell = 0;
	for (iy = 0; iy < Grid->nyC; ++iy) {
		for (ix = 0; ix < Grid->nxC; ++ix) {
			iCell = ix  + (iy  )*Grid->nxC;
			thisParticle = Particles->linkHead[iCell];




			// Loop through the particles in the cell
			// ======================================
			while (thisParticle!=NULL) {
				// Advect X
				// =====================
				locX0 = (thisParticle->x-Grid->xmin)/Grid->dx - ix;
				locY0 = (thisParticle->y-Grid->ymin)/Grid->dy - iy;

				locX = locX0*2-1.0; // important for using shape functions
				locY = locY0*2-1.0;


				if (locY>0.0) {
					locY = locY-1.0;
					Ix = ix;
					Iy = iy+1;
				}
				else {
					locY = locY+1.0;
					Ix = ix;
					Iy = iy;
				}



				//if (iP == 3) {
				//printf("iP=%i, Ix=%i, Iy=%i, locX=%.2f, locY=%.2f w0=%.3f, w1=%.3f, w2=%.3f, w3=%.3f \n",iP, Ix, Iy, locX, locY, .25*(1.0-locX)*(1.0-locY), .25*(1.0-locX)*(1.0+locY), .25*(1.0+locX)*(1.0+locY), .25*(1.0+locX)*(1.0-locY));
				//}

				thisParticle->x += ( .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
																		   + .25*(1.0-locX)*(1.0+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
																												   + .25*(1.0+locX)*(1.0+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
																																						   + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) * Physics->dt;
				/*
				if (iP == 3) {
					printf("x1=%.2f, Vtemp=%.2f, id=%i\n", Particles->xy[iP*2], ( .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
																			   + .25*(1.0-locX)*(1.0+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
																			   + .25*(1.0+locX)*(1.0+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
																			   + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ), Particles->cellId[iP]);

				}
				 */


				// Advect Y
				// =====================
				//locX = (Particles->xy[2*iP  ]-Grid->xmin)/Grid->dx - ix;
				//locY = (Particles->xy[2*iP+1]-Grid->ymin)/Grid->dy - iy;

				locX = locX0*2-1.0; // important for using shape functions
				locY = locY0*2-1.0;


				if (locX>0.0) {
					locX = locX-1.0;
					Ix = ix+1;
					Iy = iy;
				}
				else {
					locX = locX+1.0;
					Ix = ix;
					Iy = iy;
				}
				//printf("iP=%i, Ix=%i, Iy=%i, locX=%.2f, locY=%.2f w0=%.3f, w1=%.3f, w2=%.3f, w3=%.3f \n",iP, Ix, Iy, locX, locY, .25*(1.0-locX)*(1.0-locY), .25*(1.0-locX)*(1.0+locY), .25*(1.0+locX)*(1.0+locY), .25*(1.0+locX)*(1.0-locY));
				thisParticle->y  += (.25*(1.0-locX)*(1.0-locY)*Physics->Vy[Ix  +(Iy  )*Grid->nxVy]
																		   + .25*(1.0-locX)*(1.0+locY)*Physics->Vy[Ix  +(Iy+1)*Grid->nxVy]
																												   + .25*(1.0+locX)*(1.0+locY)*Physics->Vy[Ix+1+(Iy+1)*Grid->nxVy]
																																						   + .25*(1.0+locX)*(1.0-locY)*Physics->Vy[Ix+1+(Iy  )*Grid->nxVy] ) * Physics->dt;

				/*
				if (iP == 3) {
					printf("y1=%.2f, Vtemp=%.2f, id=%i\n", Particles->xy[iP*2+1], Vtemp, Particles->cellId[iP]);

				}
				 */
				thisParticle = thisParticle->next;
			}
		}
	}

}




void Particles_Periodicize(Grid* Grid, Particles* Particles, BC* BC)
{
	// Make particles do the loop

	//if (BC->VxT < BC->VxB) {
	// sinistral simple shear:
	// particles go out through the left boundary and renter through the right one
	FOR_PARTICLES
	if ((thisParticle->x-Grid->xmin)<0 ) {
		//printf("#### Loop the loop ####\n");
		thisParticle->x += Grid->xmax-Grid->xmin;
	}
	else if ((thisParticle->x-Grid->xmax)>0 ) {
		thisParticle->x -= Grid->xmax-Grid->xmin;
	}
	END_PARTICLES


}


void addSingleParticle(SingleParticle** pointerToHead, coord x, coord y, int phase, int cellId)
{
	// Adds a Particle at the beginning of a linked list
	SingleParticle* thisParticle = (SingleParticle*) malloc(sizeof(SingleParticle));
	thisParticle->x = x;
	thisParticle->y = y;
	thisParticle->phase = phase;
	thisParticle->cellId = cellId;

	thisParticle->T = 0;

	thisParticle->next = NULL;
	if (*pointerToHead != NULL) {
		thisParticle->next = *pointerToHead;
	}
	*pointerToHead = thisParticle;

}



void Particles_freeAllSingleParticles(Particles* Particles, Grid* Grid)
{
	int iCell;
	SingleParticle* temp;
	for (iCell=0;iCell<Grid->nCTot;iCell++) {
		while (Particles->linkHead[iCell] != NULL)
		{
			temp = Particles->linkHead[iCell];
			Particles->linkHead[iCell] = Particles->linkHead[iCell]->next;
			free(temp);
		}
	}


}


void addToParticlePointerList(ParticlePointerList** pointerToHead, SingleParticle* thisParticle)
{
	// Adds a node at the beginning of a linked list
	ParticlePointerList* temp = (ParticlePointerList*) malloc(sizeof(ParticlePointerList));
	temp->pointer = thisParticle;
	temp->next = NULL;
	if (*pointerToHead != NULL) {
		temp->next = *pointerToHead;
	}
	*pointerToHead = temp;

}

void freeParticlePointerList(ParticlePointerList* head)
{
	ParticlePointerList* temp;

	while (head != NULL)
	{
		temp = head;
		head = head->next;
		free(temp);
	}

}







