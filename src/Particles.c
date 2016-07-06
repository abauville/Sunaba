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

void Particles_allocateMemory(Particles* Particles, Grid* Grid) {
	Particles->linkHead 	= (SingleParticle**) malloc( Grid->nSTot 		* sizeof(  SingleParticle*  ) ); // array of pointers to particles

	int i;
	//SingleParticle* A=NULL;
	for (i=0;i<Grid->nSTot;i++) {
		Particles->linkHead[i] = NULL;
	}
}

void Particles_freeMemory(Particles* Particles, Grid* Grid) {
	printf("Free Particles..\n");
	Particles_freeAllSingleParticles(Particles, Grid);
	free( Particles->linkHead );
}




//============================================================================//
//============================================================================//
//                                                                            //
//                               INIT COORD                             	  //
//                                                                            //
//============================================================================//
//============================================================================//
void Particles_initCoord(Particles* Particles, Grid* Grid)
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

	SingleParticle modelParticle;


	modelParticle.x = 0;
	modelParticle.y = 0;
	modelParticle.nodeId = 0;
	modelParticle.T = 0;
	modelParticle.sigma_xx_0 = 0;
	modelParticle.sigma_xy_0 = 0;
	modelParticle.phase = 0;
	modelParticle.passive = 1;
	modelParticle.psi = 0;
	modelParticle.next = NULL;
	modelParticle.faulted = false;

	// Loop through nodes
	// ==================
	int iNode = 0;
	for(iy=0;iy<Grid->nyC;iy++) {
		for(ix=0;ix<Grid->nxC;ix++) {
			// Get the coordinates of the lower left corner of the shifted cell (i.e. cell centered on the node ix, iy)
			x = Grid->xmin + ix*Grid->dx;// - 0.5*Grid->dx;
			y = Grid->ymin + iy*Grid->dy;// - 0.5*Grid->dy;


			// Loop through Particles in the cell
			// ==================================
			for (iPy=0;iPy<Particles->nPCY;iPy++) {
				for (iPx=0;iPx<Particles->nPCX;iPx++) {



					// Assign coordinate
					// =================
					//printf("Rand1 = %.4f, Rand2 = %.4f\n",(0.5 - (rand() % 1000)/1000.0),(0.5 - (rand() % 1000)/1000.0));
					modelParticle.x 	= x + 0.5*dxP + iPx*dxP + Particles->noiseFactor*dxP*(0.5 - (rand() % 1000)/1000.0);
					modelParticle.y 	= y + 0.5*dyP + iPy*dyP + Particles->noiseFactor*dyP*(0.5 - (rand() % 1000)/1000.0);

					if (ix == 0 && modelParticle.x<Grid->xmin) {
						modelParticle.x += 0.5*Grid->dx;
					} else if (ix == Grid->nxC-1 && modelParticle.x>Grid->xmax){
						modelParticle.x -= 0.5*Grid->dx;
					}  else if (iy == 0  && modelParticle.y<Grid->ymin){
						modelParticle.y += 0.5*Grid->dy;
					} else if (iy == Grid->nyC-1 && modelParticle.y>Grid->ymax){
						modelParticle.y -= 0.5*Grid->dy;
					}



					iNode = (int) round((modelParticle.x-Grid->xmin)/Grid->dx) + round((modelParticle.y-Grid->ymin)/Grid->dy) * Grid->nxS;

					modelParticle.nodeId = iNode;
					// Create a particle
					if (modelParticle.x>Grid->xmin && modelParticle.x<Grid->xmax && modelParticle.y>Grid->ymin && modelParticle.y<Grid->ymax) {
						addSingleParticle(&Particles->linkHead[iNode], &modelParticle);
						//printf("iNode = %i, xP = %.3f, yP = %.3f\n", iNode, xP, yP);
					}




				} // iPx
			} // iPy


			//iNode++;
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

void Particles_initPassive(Particles* Particles, Grid* Grid)
{
	// Init a passive grid
	coord DX, DY;

	DY = (Grid->ymax-Grid->ymin)/16.0;
	DX = DY;//(Grid->xmax-Grid->xmin)/32.0;
	int passive;
	int dum;
	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, dum, passive) schedule(static,32)
	FOR_PARTICLES
	//if (thisParticle->phase>-1) {
	dum = (int)((thisParticle->x-Grid->xmin)/DX);

	passive = dum%2;
	//printf("x = %.2f, dum = %i, passive = %i\n", thisParticle->x-Grid->xmin, dum, passive);
	dum = (int)((thisParticle->y-Grid->ymin)/DY);
	passive += (dum)%2;
	if (passive==1) {
		//if (thisParticle->phase != 0) { // quick fix for sticky air visualization
		//thisParticle->passive = 0;
		//} else {
		thisParticle->passive = 0;
		//}

	} else {
		thisParticle->passive = 1;
	}
	//}
	END_PARTICLES
}









void Particles_initPhysics(Particles* Particles, Grid* Grid, BC* BCThermal)
{
	compute locY;
	compute H = (Grid->ymax-Grid->ymin);
	INIT_PARTICLE
	FOR_PARTICLES
	locY = (thisParticle->y-Grid->ymin)/H;
	thisParticle->T = 0*(  (1-locY)*BCThermal->TB + (locY)*BCThermal->TT  );


	END_PARTICLES
}





void Particles_teleportInsideTheDomain(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Due to advection error particles might end up outside the model boundaries. This function teleports them back inside
	bool change = false;

	SingleParticle* loopingParticle = NULL;
	int ParticleCounter = 0;
	compute locX, locY;
	int Imin, i, ix, iy;
	compute Min;

	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, change, ix, iy, loopingParticle, ParticleCounter, locX, locY, Min, Imin, i) schedule(static,32)

	FOR_PARTICLES
	change = false;
	if (thisParticle->x<Grid->xmin) {
		thisParticle->x = Grid->xmin+0.1*Grid->dx;
		change = true;
	} else if (thisParticle->x>Grid->xmax) {
		thisParticle->x = Grid->xmax-0.1*Grid->dx;
		change = true;
	}

	if (thisParticle->y<Grid->ymin) {
		thisParticle->y = Grid->ymin+0.1*Grid->dy;
		change = true;
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.1*Grid->dy;
		change = true;
	}







	// Interpolate new properties to the particle
	if (change==true) {

		ix = (int) round((thisParticle->x-Grid->xmin)/Grid->dx);
		iy = (int) round((thisParticle->y-Grid->ymin)/Grid->dy);

		//printf("x = %.2f y = %.2f, ix, = %i, iy = %i\n", thisParticle->x, thisParticle->y, ix, iy);

		// find the closest particle to the node
		loopingParticle = Particles->linkHead[ix+iy*Grid->nxS];

		ParticleCounter = 0;
		while (loopingParticle != NULL) {


			locX = loopingParticle->x - thisParticle->x;
			locY = loopingParticle->y - thisParticle->y;

			if ( (locX*locX + locY*locY) < Min) {
				Min = (locX*locX + locY*locY);
				Imin = ParticleCounter;
			}

			ParticleCounter++;
			loopingParticle = loopingParticle->next;
		}

		// sweep again up to the closest particle
		loopingParticle = Particles->linkHead[iNode];
		for (i=0; i<Imin; i++) {
			loopingParticle = loopingParticle->next;
		}



		thisParticle->phase = loopingParticle->phase; // the phase given to the particles is the phase of the head particle. Easy and fast but not optimal
		thisParticle->passive = loopingParticle->passive; // the phase given to the particles is the phase of the head particle. Easy and fast but not optimal

		// This could be ok, but right now it's probably done with temperature not advected or something, which gives bad results;
		//thisParticle->T = (Physics->T[(ix)+(iy+1)*Grid->nxEC] + Physics->T[ix+1+(iy+1)*Grid->nxEC] + Physics->T[(ix)+(iy)*Grid->nxEC] + Physics->T[ix+1    +(iy)*Grid->nxEC])/4;
		//thisParticle->sigma_xx_0 = (Physics->sigma_xx_0[(ix)+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[(ix)+(iy)*Grid->nxEC] + Physics->sigma_xx_0[ix+1    +(iy)*Grid->nxEC])/4;


		thisParticle->T = loopingParticle->T;
		thisParticle->sigma_xx_0 = loopingParticle->sigma_xx_0;


		thisParticle->sigma_xy_0 = loopingParticle->sigma_xy_0; // not ideal


		thisParticle->nodeId = ix+iy*Grid->nxS;










		//thisParticle->sigma_xx_0 = 0;
		//thisParticle->sigma_xy_0 = 0;
	}

	END_PARTICLES



}

void Particles_deleteIfOutsideTheDomain(Particles* Particles, Grid* Grid)
{
	//printf("In deletion\n");
	SingleParticle* nextParticle = NULL;
	// Due to advection error particles might end up outside the model boundaries. This function teleports them back inside
	SingleParticle* thisParticle = NULL;
	//SingleParticle* prevParticle = NULL;
	int iNode = 0;

	int justDeleted = 0;
	bool change = false;
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		thisParticle = Particles->linkHead[iNode];

		while (thisParticle != NULL) {

			/*
			if (thisParticle->x>0 && thisParticle!=Particles->linkHead[iNode]) {
				printf("============= thisP: x = %.3f, y = %.3f, C = %i, D = %i\n", thisParticle->x, thisParticle->y, C, D);
				exit(0);
			}

			 */
			nextParticle = thisParticle->next;

			if (nextParticle!=NULL) {
				/*
				if (nextParticle->x>0) {
					printf("D: x = %.3f, y = %.3f, C = %i, D = %i\n", nextParticle->x, nextParticle->y, C, D);
				}
				 */
				//printf("nextP: x = %.3f, y = %.3f, C = %i, D = %i\n", nextParticle->x, nextParticle->y, C, D);
				if (nextParticle->x<Grid->xmin || nextParticle->x>Grid->xmax
						|| nextParticle->y<Grid->ymin || nextParticle->y>Grid->ymax	) {
					thisParticle->next = nextParticle->next;
					//printf("Particle to be deleted: x = %.3f, y = %.3f\n", nextParticle->x, nextParticle->y);
					free(nextParticle);
					Particles->n-=1;

					justDeleted = 1;
				}

			}

			if (justDeleted==1) {
				justDeleted = 0;
			} else {
				thisParticle = thisParticle->next;
			}



		}

		// check if head is out
		thisParticle = Particles->linkHead[iNode];
		if (thisParticle == NULL) {
			printf("error in Particles_deleteIfOutsideTheDomain: Empty node\n");
			exit(0);
		}
		change = false;
		if (thisParticle->x<Grid->xmin) {
			thisParticle->x = Grid->xmin+0.1*Grid->dx;
			change = true;
		} else if (thisParticle->x>Grid->xmax) {
			thisParticle->x = Grid->xmax-0.1*Grid->dx;
			change = true;
		}

		if (thisParticle->y<Grid->ymin) {
			thisParticle->y = Grid->ymin+0.1*Grid->dy;
			change = true;
		} else if (thisParticle->y>Grid->ymax) {
			thisParticle->y = Grid->ymax-0.1*Grid->dy;
			change = true;
		}



		// Interpolate new properties to the particle
		if (change==true) {


			thisParticle->sigma_xx_0 = 0;
			thisParticle->sigma_xy_0 = 0;
		}


		/*
		if (thisParticle->x<Grid->xmin || thisParticle->x>Grid->xmax
				|| thisParticle->y<Grid->ymin || thisParticle->y>Grid->ymax	) {
			Particles->linkHead[iNode] = thisParticle->next;
			free(thisParticle);
			Particles->n-=1;
		}
		 */


	}

	//printf("End deletion");

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


void Particles_updateLinkedList(Particles* Particles, Grid* Grid, Physics* Physics)
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
	ParticlePointerList* headIdChanged = (ParticlePointerList*) malloc(sizeof(ParticlePointerList));
	headIdChanged->pointer = NULL;
	headIdChanged->next = NULL;

	int oldNodeId;
	coord x, y;
	int ix, iy;
	SingleParticle* previousParticle;
	// Update the link list
	// =========================


	int ParticleCounter = 0;
	SingleParticle* thisParticle = NULL;
	int iNode = 0;

	printf("First loop\n");

	//#pragma omp parallel for private(iNode, thisParticle, ParticleCounter, oldNodeId, x, y, ix, iy, previousParticle) schedule(static,32)
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		thisParticle = Particles->linkHead[iNode];
		ParticleCounter = 0;
		while (thisParticle != NULL) {

			ParticleCounter++;


			oldNodeId = thisParticle->nodeId;



			x = thisParticle->x;
			y = thisParticle->y;


			ix = (int) round((x-Grid->xmin)/Grid->dx);
			iy = (int) round((y-Grid->ymin)/Grid->dy);


			//printf("x = %.1f , Grid->xmin = %.1f", Particles->xy[2*iP],Grid->xmin );

			thisParticle->nodeId = ix + iy*Grid->nxS;




			//printf("iP=%i, oid=%i, nid=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,oldCellId, Particles->cellId[iP],x, y, ix,iy);
			// If this particle has changed cell
			if (oldNodeId != thisParticle->nodeId) {
				//printf("iP=%i, oid=%i, nid=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,oldCellId, Particles->cellId[iP],x, y, ix,iy);
				// 1. Update info for the oldCell
				// ===========================
				if (thisParticle != Particles->linkHead[iNode]) {
					previousParticle->next = thisParticle->next;
				}
				else {
					Particles->linkHead[iNode] = thisParticle->next;
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
	printf("Update info\n");
	// 2. Update info of the new cell, i.e. Add this particle to the head of the link list of the new cell
	// ==============================
	ParticlePointerList* IdChanged = NULL;
	IdChanged = headIdChanged;
	while (IdChanged->next!=NULL) {
		thisParticle 	= IdChanged->pointer;
		IdChanged 		= IdChanged->next;

		ix = (int) round((thisParticle->x-Grid->xmin)/Grid->dx);
		iy = (int) round((thisParticle->y-Grid->ymin)/Grid->dy);

		thisParticle->nodeId = ix + iy*Grid->nxS;


		thisParticle->next = Particles->linkHead[thisParticle->nodeId] ;
		Particles->linkHead[thisParticle->nodeId] = thisParticle;
	}
	freeParticlePointerList(headIdChanged);













































	if (DEBUG) {
		// Check implementation
		// ====================
		/*
		for (iP = 0; iP < Particles->n; ++iP) {
			printf("cellId = %i, iP = %i, Next = %i\n",Particles->cellId[iP], iP,  Particles->linkNext[iP]);
		}
		for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
			printf("%i  ", Particles->linkHead[iNode]);
		}
		printf("\n\n\n");
		 */

		/*
		for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
			printf("Cell #%i:  ", iNode);
			thisParticle = Particles->linkHead[iNode];
			while (thisParticle!=NULL) {
				printf("%i  ",iP);
				thisParticle = thisParticle->next;
			}
			printf("\n");
		}
		 */
	}
}






void Particles_injectOrDelete(Particles* Particles, Grid* Grid)
{
	// If the node is empty add a particle at the node with the same values as the closest one
	// only the neighbour cells up, down, left and right are checked for neighbour particles


	int IxN[5], IyN[5];
	int iNodeNeigh;
	int ix, iy, i, iNode;
	int numPart;

	compute x, y;

	SingleParticle* thisParticle 	= NULL;
	SingleParticle* closestParticle = NULL;
	SingleParticle* neighParticle 	= NULL;
//	SingleParticle* thisParticle = NULL;
	i = 0;
	while(Particles->linkHead[i]==NULL) {
		i++;
	}
	SingleParticle* modelParticle = Particles->linkHead[i];

	compute dist, minDist;
	printf("Start injection loop\n");
	int iBlock; //loop index for left, right, up, down sides + inner
	int ix0, ixMax, iy0, iyMax;
	compute xMod, yMod;
	int nNeighbours;

	compute minNumPart = Particles->nPCX*Particles->nPCY*Particles->minPartPerCellFactor;
	compute maxNumPart = Particles->nPCX*Particles->nPCY*Particles->maxPartPerCellFactor;

	int* PartAdded = (int*) malloc(Grid->nSTot*sizeof(int));
	for (i = 0; i < Grid->nSTot; ++i) {
		PartAdded[i] = 0;
	}


	for (iBlock = 0; iBlock<9;++iBlock) {
		// note:: all sides are of length of nodes-1 and the xMod and yMod are shifted so that even in the corners, the new particle is not on a side
		switch (iBlock) {
		case 0: // Inner nodes
			iy0 = 1;
			iyMax = Grid->nyS-1;
			ix0 = 1;
			ixMax = Grid->nxS-1;
			IxN[0] =  0; IyN[0] =  0;
			IxN[1] = -1; IyN[1] =  0;
			IxN[2] =  1; IyN[2] =  0;
			IxN[3] =  0; IyN[3] = -1;
			IxN[4] =  0; IyN[4] =  1;
			nNeighbours = 5;
			xMod = 0; yMod = 0;
			break;
		case 1: // inner lower nodes
			iy0 = 0;
			iyMax = 1;
			ix0 = 1;
			ixMax = Grid->nxS-1;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =  -1; IyN[1] =  0;
			IxN[2] =   1; IyN[2] =  0;
			IxN[3] =   0; IyN[3] =  1;
			nNeighbours = 4;
			xMod = 0; yMod =  0.25*Grid->dy;
			break;
		case 2: // inner upper nodes
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0 = 1;
			ixMax = Grid->nxS-1;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =  -1; IyN[1] =  0;
			IxN[2] =   1; IyN[2] =  0;
			IxN[3] =   0; IyN[3] = -1;
			nNeighbours = 4;
			xMod = 0; yMod = -0.25*Grid->dy;
			break;
		case 3: // inner left nodes
			iy0 = 1;
			iyMax = Grid->nyS-1;
			ix0 = 0;
			ixMax = 1;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =   0; IyN[1] = -1;
			IxN[2] =   0; IyN[2] =  1;
			IxN[3] =   1; IyN[3] =  0;
			nNeighbours = 4;
			xMod =  0.25*Grid->dx; yMod = 0;
			break;
		case 4: // inner right nodes
			iy0 = 1;
			iyMax = Grid->nyS-1;
			ix0 = Grid->nxS-1;
			ixMax = Grid->nxS;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =   0; IyN[1] = -1;
			IxN[2] =   0; IyN[2] =  1;
			IxN[3] =  -1; IyN[3] =  0;
			nNeighbours = 4;
			xMod = -0.25*Grid->dx; yMod =  0;
			break;
		case 5: // upper left corner
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0  = 0;
			ixMax = 1;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =   1; IyN[1] =  0;
			IxN[2] =   0; IyN[2] = -1;
			IxN[3] =   1; IyN[3] = -1;
			nNeighbours = 4;
			xMod = 0.25*Grid->dx; yMod = -0.25*Grid->dy;
			break;
		case 6: // upper right corner
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0  = Grid->nxS-1;
			ixMax = Grid->nxS;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =  -1; IyN[1] =  0;
			IxN[2] =   0; IyN[2] = -1;
			IxN[3] =  -1; IyN[3] = -1;
			nNeighbours = 4;
			xMod = -0.25*Grid->dx; yMod = -0.25*Grid->dy;
			break;
		case 7: // lower right corner
			iy0 = 0;
			iyMax = 1;
			ix0  = Grid->nxS-1;
			ixMax = Grid->nxS;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =  -1; IyN[1] =  0;
			IxN[2] =   0; IyN[2] =  1;
			IxN[3] =  -1; IyN[3] =  1;
			nNeighbours = 4;
			xMod = -0.25*Grid->dx; yMod = 0.25*Grid->dy;
			break;
		case 8: // lower left corner
			iy0 = 0;
			iyMax = 1;
			ix0  = 0;
			ixMax = 1;
			IxN[0] =   0; IyN[0] =  0;
			IxN[1] =   1; IyN[1] =  0;
			IxN[2] =   0; IyN[2] = 1;
			IxN[3] =   1; IyN[3] = 1;
			nNeighbours = 4;
			xMod =  0.25*Grid->dx; yMod = 0.25*Grid->dy;
			break;

		}
#pragma omp parallel for private(iy, ix, iNode, thisParticle, numPart, i, minDist, x, y, iNodeNeigh, neighParticle, dist, closestParticle) schedule(static,32)
		for (iy = iy0; iy < iyMax; ++iy) {
			for (ix = ix0; ix < ixMax; ++ix) {
				iNode = ix + iy*Grid->nxS;

				numPart = 0;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle != NULL && numPart<minNumPart) {
					thisParticle = thisParticle->next;
					++numPart;
				}


				if (numPart<minNumPart) {
					//printf("************* A particle is about to be injected!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ****************\n");
					minDist = (100*Grid->dx)*(100*Grid->dx);

					x = Grid->xmin + ix*Grid->dx + xMod;
					y = Grid->ymin + iy*Grid->dy + yMod;


					for (i=0;i<nNeighbours;i++) {
						iNodeNeigh = ix+IxN[i] + (iy+IyN[i])*Grid->nxS;
						neighParticle = Particles->linkHead[iNodeNeigh];
						while (neighParticle != NULL) {
							dist = (neighParticle->x - x)*(neighParticle->x - x) + (neighParticle->y - y)*(neighParticle->y - y);
							if (dist<minDist) {
								closestParticle = neighParticle;
								minDist = dist;
							}

							neighParticle = neighParticle->next;
						}
					}


					addSingleParticle(&Particles->linkHead[iNode], closestParticle);
					Particles->linkHead[iNode]->x = x;
					Particles->linkHead[iNode]->y = y;
					Particles->linkHead[iNode]->nodeId = iNode;
					PartAdded[iNode] += 1;


				}

				/*
				else if (numPart>maxNumPart) {
					//printf("Delete part\n");
					thisParticle = Particles->linkHead[iNode];
					Particles->linkHead[iNode] = Particles->linkHead[iNode]->next;
					free(thisParticle);
					PartAdded[iNode] -= 1;
				}
				*/


			}
		}
	}
	//printf("Out\n");


	for (i = 0; i < Grid->nSTot; ++i) {
		Particles->n += PartAdded[i];
	}
	free(PartAdded);
}










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

	// Index of neighbouring cells, with respect to the node ix, iy
	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  0;	IyN[1] =  1; // upper left
	IxN[2] =  1; 	IyN[2] =  1; // upper right
	IxN[3] =  1; 	IyN[3] =  0; // lower right

	int signX, signY;

	// Loop through inner cells
	// ========================
	iNode = 0;
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX0, locY0, locX, locY, signX, signY, i, ixN, iyN, alphaArray, alpha, sigma_xx_temp, Ix, Iy) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];




			// Loop through the particles in the cell
			// ======================================
			while (thisParticle!=NULL) {


				// Advect X
				// =====================
				locX0 = (thisParticle->x-Grid->xmin)/Grid->dx - ix;
				locY0 = (thisParticle->y-Grid->ymin)/Grid->dy - iy;

				locX = locX0*2.0; // important for using shape functions
				locY = locY0*2.0;


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

					alphaArray[i]  = - 0.5*Physics->dtAdv*((Physics->Vy[ixN+1+iyN*Grid->nxVy]   - Physics->Vy[ixN+(iyN)*Grid->nxVy])/Grid->dx
							- (Physics->Vx[ixN+(iyN+1)*Grid->nxVx] - Physics->Vx[ixN+(iyN)*Grid->nxVx])/Grid->dy);
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






				locX = locX0*2.0; // important for using shape functions
				locY = locY0*2.0;

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



				thisParticle->x += ( .25*(1.0-locX)*(1.0-locY)*Physics->Vx[Ix  +(Iy  )*Grid->nxVx]
								   + .25*(1.0-locX)*(1.0+locY)*Physics->Vx[Ix  +(Iy+1)*Grid->nxVx]
								   + .25*(1.0+locX)*(1.0+locY)*Physics->Vx[Ix+1+(Iy+1)*Grid->nxVx]
								   + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) * Physics->dtAdv;


				// Advect Y
				// =====================


				locX = locX0*2; // important for using shape functions
				locY = locY0*2;


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

				thisParticle->y  += (.25*(1.0-locX)*(1.0-locY)*Physics->Vy[Ix  +(Iy  )*Grid->nxVy]
							       + .25*(1.0-locX)*(1.0+locY)*Physics->Vy[Ix  +(Iy+1)*Grid->nxVy]
					               + .25*(1.0+locX)*(1.0+locY)*Physics->Vy[Ix+1+(Iy+1)*Grid->nxVy]
								   + .25*(1.0+locX)*(1.0-locY)*Physics->Vy[Ix+1+(Iy  )*Grid->nxVy] ) * Physics->dtAdv;


				thisParticle = thisParticle->next;
			}
		}
	}

}




void Particles_Periodicize(Particles* Particles, Grid* Grid)
{
	// Make particles do the loop

	//if (BC->VxT < BC->VxB) {
	// sinistral simple shear:
	// particles go out through the left boundary and renter through the right one
	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle) schedule(static,32)
	FOR_PARTICLES
	if (thisParticle->x<Grid->xmin ) {
		thisParticle->x += Grid->xmax-Grid->xmin;
	}
	else if (thisParticle->x>Grid->xmax ) {
		thisParticle->x -= Grid->xmax-Grid->xmin;
	}

	if (thisParticle->y<Grid->ymin) {
		thisParticle->y = Grid->ymin+0.05*Grid->dy;
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.05*Grid->dy;
	}
	END_PARTICLES


}


void addSingleParticle(SingleParticle** pointerToHead, SingleParticle* modelParticle)
{
	// Adds a Particle at the beginning of a linked list
	SingleParticle* thisParticle = (SingleParticle*) malloc(sizeof(SingleParticle));
	thisParticle->x = modelParticle->x;
	thisParticle->y = modelParticle->y;
	thisParticle->phase = modelParticle->phase;
	thisParticle->passive = modelParticle->passive;
	thisParticle->nodeId = modelParticle->nodeId;

	thisParticle->T = modelParticle->T;
	thisParticle->sigma_xx_0 = modelParticle->sigma_xx_0;
	thisParticle->sigma_xy_0 = modelParticle->sigma_xy_0;

	thisParticle->psi = modelParticle->psi;

	thisParticle->faulted = modelParticle->faulted;


	thisParticle->next = NULL;
	if (*pointerToHead != NULL) {
		thisParticle->next = *pointerToHead;
	}
	*pointerToHead = thisParticle;

}



void Particles_freeAllSingleParticles(Particles* Particles, Grid* Grid)
{
	int iNode;
	SingleParticle* temp;
	for (iNode=0;iNode<Grid->nSTot;iNode++) {
		while (Particles->linkHead[iNode] != NULL)
		{
			temp = Particles->linkHead[iNode];
			Particles->linkHead[iNode] = Particles->linkHead[iNode]->next;
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







