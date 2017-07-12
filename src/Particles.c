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


void findNodeForThisParticle(SingleParticle* thisParticle, Grid* Grid);


void Particles_allocateMemory(Particles* Particles, Grid* Grid) {
	Particles->linkHead 	= (SingleParticle**) malloc( Grid->nSTot 		* sizeof(  SingleParticle*  ) ); // array of pointers to particles

	Particles->dispAtBoundL = (compute*) malloc(Grid->nyS * sizeof(compute));
	Particles->dispAtBoundR = (compute*) malloc(Grid->nyS * sizeof(compute));
	Particles->currentPassiveAtBoundL = (int*) malloc(Grid->nyS * sizeof(int));
	Particles->currentPassiveAtBoundR = (int*) malloc(Grid->nyS * sizeof(int));

	int i;
	//SingleParticle* A=NULL;
	for (i=0;i<Grid->nSTot;i++) {
		Particles->linkHead[i] = NULL;
	}

	for (i=0;i<Grid->nyS;i++) {
		Particles->dispAtBoundL[i] = 0.0;
		Particles->dispAtBoundR[i] = 0.0;
		Particles->currentPassiveAtBoundL[i] = 0;
		Particles->currentPassiveAtBoundR[i] = 0;
	}
}

void Particles_freeMemory(Particles* Particles, Grid* Grid) {
	printf("Free Particles..\n");
	Particles_freeAllSingleParticles(Particles, Grid);
	free( Particles->linkHead );
	free(Particles->dispAtBoundL);
	free(Particles->dispAtBoundR);
	free(Particles->currentPassiveAtBoundL);
	free(Particles->currentPassiveAtBoundR);
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
	coord dxP;
	coord dyP;


	// Init random number generator
	// ==================
	srand(time(NULL));

	SingleParticle modelParticle;


	modelParticle.x = 0;
	modelParticle.y = 0;
	modelParticle.nodeId = 0;

	modelParticle.sigma_xx_0 = 0;
	modelParticle.sigma_xy_0 = 0;
	modelParticle.phase = 0;
	modelParticle.passive = 1;
	modelParticle.Vx = 0.0;
	modelParticle.Vy = 0.0;
#if (STRAIN_SOFTENING)
	modelParticle.strain = 0.0;
#endif


	modelParticle.next = NULL;
#if (HEAT)
	modelParticle.T = 0.0;
#endif
#if (DARCY)
	modelParticle.DeltaP0 = 0;
	modelParticle.phi = 0;
#endif
	//modelParticle.faulted = false;

	// Loop through nodes
	// ==================
	int iNode = 0;
	for(iy=0;iy<Grid->nyC;iy++) {
		for(ix=0;ix<Grid->nxC;ix++) {
			iNode = ix + iy*Grid->nxS;
			// Get the coordinates of the lower left corner of the shifted cell (i.e. cell centered on the node ix, iy)
			dxP = Grid->DXS[ix]/Particles->nPCX;
			dyP = Grid->DYS[iy]/Particles->nPCY;
			x = Grid->X[ix];// - 0.5*Grid->dx;
			y = Grid->Y[iy];// - 0.5*Grid->dy;


			// Loop through Particles in the cell
			// ==================================
			for (iPy=0;iPy<Particles->nPCY;iPy++) {
				for (iPx=0;iPx<Particles->nPCX;iPx++) {



					// Assign coordinate
					// =================
					//printf("Rand1 = %.4f, Rand2 = %.4f\n",(0.5 - (rand() % 1000)/1000.0),(0.5 - (rand() % 1000)/1000.0));
					modelParticle.x 	= x + 0.5*dxP + iPx*dxP + Particles->noiseFactor*dxP*(0.5 - (rand() % 1000)/1000.0);
					modelParticle.y 	= y + 0.5*dyP + iPy*dyP + Particles->noiseFactor*dyP*(0.5 - (rand() % 1000)/1000.0);


					//printf("x = %.2f, y = %.2f\n",modelParticle.x, modelParticle.y);
					//exit(0);

					if (ix == 0 && modelParticle.x<Grid->xmin) {
						modelParticle.x += 0.5*Grid->DXS[0];
					} else if (ix == Grid->nxC-1 && modelParticle.x>Grid->xmax){
						modelParticle.x -= 0.5*Grid->DXS[Grid->nxC-1];
					}  else if (iy == 0  && modelParticle.y<Grid->ymin){
						modelParticle.y += 0.5*Grid->DYS[0];
					} else if (iy == Grid->nyC-1 && modelParticle.y>Grid->ymax){
						modelParticle.y -= 0.5*Grid->DYS[Grid->nyC-1];
					}

#if (STORE_PARTICLE_POS_INI)
					modelParticle.xIni = modelParticle.x;
					modelParticle.yIni = modelParticle.y;
#endif

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

void Particles_initPassive(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Init a passive grid
	coord DX, DY;
	if (Particles->passiveGeom==PartPassive_Grid || Particles->passiveGeom==PartPassive_Grid_w_Layers) {
		DY = Particles->passiveDy;//(Grid->ymax-Grid->ymin)*Particles->passiveRes;
		DX = Particles->passiveDx;//DY;//(Grid->xmax-Grid->xmin)/32.0;
		int passive;
		int dum;


		int iy;
		compute x;
		compute y;
		for (iy = 0; iy < Grid->nyS; ++iy) {
			y = Grid->ymin + (iy)*Grid->dy;
			// Left boundary
			x = Grid->xmin;
			dum = (int)((x-Grid->xmin)/DX);
			passive = dum%2;
			dum = (int)((y-Grid->ymin)/DY);
			passive += (dum)%2;
			if (passive==1) {
				Particles->currentPassiveAtBoundL[iy] = 0;
			} else {
				Particles->currentPassiveAtBoundL[iy] = 1;
			}
			if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
				dum = (int)((y-Grid->ymin)/DY);
				passive = (dum)%2;
				if (passive==1) {
					Particles->currentPassiveAtBoundL[iy] += 2;
				} else {
					Particles->currentPassiveAtBoundL[iy] += 0;
				}
			}
			Particles->dispAtBoundL[iy] = DX;

			// Right boundary
			x = Grid->xmax;
			dum = (int)((x-Grid->xmin)/DX);
			Particles->dispAtBoundR[iy] = (Grid->xmax-Grid->xmin) - dum*DX;
			passive = dum%2;
			dum = (int)((y-Grid->ymin)/DY);
			passive += (dum)%2;
			if (passive==1) {
				Particles->currentPassiveAtBoundR[iy] = 0;
			} else {
				Particles->currentPassiveAtBoundR[iy] = 1;
			}
			if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
				dum = (int)((y-Grid->ymin)/DY);
				passive = (dum)%2;
				if (passive==1) {
					Particles->currentPassiveAtBoundR[iy] += 2;
				} else {
					Particles->currentPassiveAtBoundR[iy] += 0;
				}
			}

			//printf("iy = %i, dispL = %.2e, dispR = %.2e, passL = %i, passR = %i\n",iy, Particles->dispAtBoundL[iy], Particles->dispAtBoundR[iy], Particles->currentPassiveAtBoundL[iy], Particles->currentPassiveAtBoundR[iy]);

		}

		Particles->currentPassiveAtBoundR[Grid->nyS-1] = Particles->currentPassiveAtBoundR[Grid->nyS-2]; // overwrite the uppermost one to avoid a disgracious passive switch just at the corner
		Particles->currentPassiveAtBoundL[Grid->nyS-1] = Particles->currentPassiveAtBoundL[Grid->nyS-2];
		Particles->currentPassiveAtBoundR[0] = Particles->currentPassiveAtBoundR[1]; // overwrite the uppermost one to avoid a disgracious passive switch just at the corner
		Particles->currentPassiveAtBoundL[0] = Particles->currentPassiveAtBoundL[1];



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



		END_PARTICLES
	}




	if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
		int dum, passive;

		INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, dum, passive) schedule(static,32)
		FOR_PARTICLES
		//if (thisParticle->phase>-1) {
		//dum = (int)((thisParticle->x-Grid->xmin)/DX);

		//passive = dum%2;
		//printf("x = %.2f, dum = %i, passive = %i\n", thisParticle->x-Grid->xmin, dum, passive);
		dum = (int)((thisParticle->y-Grid->ymin)/DY);
		passive = (dum)%2;
		if (passive==1) {
			//if (thisParticle->phase != 0) { // quick fix for sticky air visualization
			//thisParticle->passive = 0;
			//} else {
			thisParticle->passive += 2;
			//}

		} else {
			thisParticle->passive += 0;
		}



		END_PARTICLES
	}
}













void Particles_teleportInsideTheDomain(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Due to advection error particles might end up outside the model boundaries. This function teleports them back inside
	//bool change = false;

	//SingleParticle* loopingParticle = NULL;
	//int ParticleCounter = 0;
	//compute locX, locY;
	//int Imin, i, ix, iy;
	//compute Min;

	INIT_PARTICLE
	//#pragma omp parallel for private(iNode, thisParticle, change, ix, iy, loopingParticle, ParticleCounter, locX, locY, Min, Imin, i) schedule(static,32)

	FOR_PARTICLES

	//change = false;
	if (thisParticle->x<Grid->xmin) {
		thisParticle->x = Grid->xmin+0.05*Grid->DXEC[0];
		//change = true;
	} else if (thisParticle->x>Grid->xmax) {
		thisParticle->x = Grid->xmax-0.05*Grid->DXEC[Grid->nxS-1];
		//change = true;
	}

	if (thisParticle->y<Grid->ymin) {
		thisParticle->y = Grid->ymin+0.05*Grid->DYEC[0];
		//change = true;
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.05*Grid->DYEC[Grid->nyS-1];
		//change = true;
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
			thisParticle->x = Grid->xmin+0.1*Grid->DXEC[0];
			change = true;
		} else if (thisParticle->x>Grid->xmax) {
			thisParticle->x = Grid->xmax-0.1*Grid->DXEC[Grid->nxS-1];
			change = true;
		}

		if (thisParticle->y<Grid->ymin) {
			thisParticle->y = Grid->ymin+0.1*Grid->DYEC[0];
			change = true;
		} else if (thisParticle->y>Grid->ymax) {
			thisParticle->y = Grid->ymax-0.1*Grid->DYEC[Grid->nyS-1];
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
	//coord x, y;
	//int ix, iy;
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


			//printf("x = %.1f , Grid->xmin = %.1f", Particles->xy[2*iP],Grid->xmin );


			findNodeForThisParticle(thisParticle, Grid);
			//if (iNode == 1)



			//printf("iP=%i, oid=%i, nid=%i, x=%.2f, y=%.2f, ix=%i, iy=%i\n",iP,oldCellId, Particles->cellId[iP],x, y, ix,iy);
			// If this particle has changed cell
			if (oldNodeId != thisParticle->nodeId) {
				//printf("part not ok!, oldNoedId = %i, thisParticle->nodeId = %i\n", oldNodeId, thisParticle->nodeId);
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
				//printf("part ok!, oldNoedId = %i, thisParticle->nodeId = %i\n", oldNodeId, thisParticle->nodeId);
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
	//int i;
	//printf("End loop\n");
	printf("Update info\n");
	// 2. Update info of the new cell, i.e. Add this particle to the head of the link list of the new cell
	// ==============================
	ParticlePointerList* IdChanged = NULL;
	IdChanged = headIdChanged;
	while (IdChanged->pointer!=NULL) {
		thisParticle 	= IdChanged->pointer;


		//printf("pointer = %d\n", IdChanged->pointer);
		/*
		ix = (int) round((thisParticle->x-Grid->xmin)/Grid->DXEC[0]);
		iy = (int) round((thisParticle->y-Grid->ymin)/Grid->DYEC[0]);

		thisParticle->nodeId = ix + iy*Grid->nxS;
		 */


		//printf("koko\n");
		//printf("thisParticle->nodeId = %d. x =%.2e, y = %.2e\n", thisParticle->nodeId, thisParticle->x, thisParticle->y);
		/*
		if (isnan(thisParticle->x)!=0 || isnan(thisParticle->y)!=0 ) {
			printf("nan coord: x = %.2e, y = %.2e\n", thisParticle->x, thisParticle->y);
		}
		*/
		thisParticle->next = Particles->linkHead[thisParticle->nodeId] ;
		//printf("soko\n");
		Particles->linkHead[thisParticle->nodeId] = thisParticle;
		//printf("asoko\n");
		IdChanged 		= IdChanged->next;
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
	compute numPart;

	compute x, y;

	SingleParticle* thisParticle 	= NULL;
	SingleParticle* closestParticle = NULL;
	SingleParticle* neighParticle 	= NULL;
	//	SingleParticle* thisParticle = NULL;
	i = 0;
	while(Particles->linkHead[i]==NULL) {
		i++;
	}

	compute dist, minDist;
	printf("Start injection loop\n");
	int iBlock; //loop index for left, right, up, down sides + inner
	int ix0, ixMax, iy0, iyMax;
	compute xMod, yMod;
	int nNeighbours;

	compute minNumPart = Particles->nPCX*Particles->nPCY*Particles->minPartPerCellFactor;
	//compute maxNumPart = Particles->nPCX*Particles->nPCY*Particles->maxPartPerCellFactor;

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
			xMod = 0; yMod =  0.125*Grid->DYEC[0];
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
			xMod = 0; yMod = -0.125*Grid->DYEC[Grid->nyS-1];
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
			xMod =  0.125*Grid->DXEC[0]; yMod = 0;
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
			xMod = -0.125*Grid->DXEC[Grid->nxS-1]; yMod =  0;
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
			xMod = 0.125*Grid->DXEC[0]; yMod = -0.125*Grid->DYEC[Grid->nyS-1];
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
			xMod = -0.125*Grid->DXEC[Grid->nxS-1]; yMod = -0.125*Grid->DYEC[Grid->nyS-1];
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
			xMod = -0.125*Grid->DXEC[0]; yMod = 0.125*Grid->DYEC[Grid->nyS-1];
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
			xMod =  0.125*Grid->DXEC[0]; yMod = 0.125*Grid->DYEC[0];
			break;

		}
		//#pragma omp parallel for private(iy, ix, iNode, thisParticle, numPart, i, minDist, x, y, iNodeNeigh, neighParticle, dist, closestParticle) schedule(static,32)
		for (iy = iy0; iy < iyMax; ++iy) {
			for (ix = ix0; ix < ixMax; ++ix) {
				iNode = ix + iy*Grid->nxS;

				numPart = 0.;
				thisParticle = Particles->linkHead[iNode];
				while (thisParticle != NULL && numPart<minNumPart) {
					thisParticle = thisParticle->next;
					numPart += 1.;
				}


				if (numPart<minNumPart) {
					//printf("in\n");
					//printf("************* A particle is about to be injected!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ****************\n");
					minDist = (Grid->xmax-Grid->xmin)*(Grid->xmax-Grid->xmin);

					x = Grid->X[ix] + xMod;
					y = Grid->Y[iy] + yMod;
					//printf("A\n");

					for (i=0;i<nNeighbours;i++) {
						iNodeNeigh = ix+IxN[i] + (iy+IyN[i])*Grid->nxS;
						//printf("iNode = %i, iNodeNeigh = %i\n",iNode, iNodeNeigh);
						neighParticle = Particles->linkHead[iNodeNeigh];
						while (neighParticle != NULL) {
							dist = (neighParticle->x - x)*(neighParticle->x - x) + (neighParticle->y - y)*(neighParticle->y - y);
							//printf("dist/dx = %.2e, neighParticle->phase = %i\n",dist/Grid->dx, neighParticle->phase);
							if (dist<minDist) {
								closestParticle = neighParticle;
								minDist = dist;
							}

							neighParticle = neighParticle->next;
						}
					}
					//printf("B, minDist = %.2e, closestPartId = %d, neighParticle = %d, dist = %.2e, iNodeNeigh = %i, ixN = %i, iyN = %i\n", minDist, closestParticle, neighParticle, dist, iNodeNeigh, ix+IxN[nNeighbours-1], iy+IyN[nNeighbours-1]);
					if (closestParticle!=NULL) {
					//printf("closestParticle->phase = %i\n",closestParticle->phase);
					addSingleParticle(&Particles->linkHead[iNode], closestParticle);
					Particles->linkHead[iNode]->x = x;
					Particles->linkHead[iNode]->y = y;
					Particles->linkHead[iNode]->nodeId = iNode;
					PartAdded[iNode] += 1;
					//printf("out\n");
					} else {
						printf("Error, No closest particle attached to neighbour nodes. Injection failed\n");
						exit(0);
					}

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





void Particles_injectAtTheBoundaries(Particles* Particles, Grid* Grid, Physics* Physics, MatProps* MatProps)
{
	// If the node is empty add a particle at the node with the same values as the closest one
	// only the neighbour cells up, down, left and right are checked for neighbour particles





	int ix, iy, i, iNode;
	compute numPart;

	compute x, y;

	SingleParticle* thisParticle 	= NULL;
	SingleParticle* closestParticle = NULL;
	//SingleParticle* modelParticle = NULL;
	SingleParticle* neighParticle 	= NULL;
	SingleParticle* temp 	= NULL;
	//	SingleParticle* thisParticle = NULL;
	i = 0;
	while(Particles->linkHead[i]==NULL) {
		i++;
	}

	//compute dist, minDist;
	printf("Start injection loop At the boundaries\n");
	int iBlock; //loop index for left, right, up, down sides + inner
	int ix0, ixMax, iy0, iyMax;
	compute xMod1, xMod2, yMod1, yMod2;
	//int nNeighbours;

	//compute minNumPart = Particles->nPCX*Particles->nPCY*Particles->minPartPerCellFactor;
	compute minNumPart = Particles->nPCX*Particles->nPCY/2.0;
	//compute maxNumPart = Particles->nPCX*Particles->nPCY*Particles->maxPartPerCellFactor;

	int* PartAdded = (int*) malloc(Grid->nSTot*sizeof(int));
	for (i = 0; i < Grid->nSTot; ++i) {
		PartAdded[i] = 0;
	}

	bool forcePassive;
	float passive;

	srand(time(NULL));

	int nNeighbours, iNodeNeigh, IxN, IyN;
	compute dist, minDist;

	compute Vx;
	bool inject;

	int Method = 0; // 0: copy particles from the neighbour cells; 1: inject a single particle

	for (iBlock = 0; iBlock<9;++iBlock) {
		// note:: all sides are of length of nodes-1 and the xMod and yMod are shifted so that even in the corners, the new particle is not on a side
		switch (iBlock) {
		case 0: // inner lower nodes
			iy0 = 0;
			iyMax = 1;
			ix0 = 1;
			ixMax = Grid->nxS-1;
			xMod1 =  1.0;
			xMod2 =  1.0;
			yMod1 =  0.0;
			yMod2 =  0.5;
			IxN =   0; IyN =  1;
			break;
		case 1: // inner upper nodes
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0 = 1;
			ixMax = Grid->nxS-1;
			xMod1 =  1.0;
			xMod2 =  1.0;
			yMod1 =  0.0;
			yMod2 = -0.5;
			IxN =   0; IyN =  -1;
			break;
		case 2: // inner left nodes
			iy0 = 1;
			iyMax = Grid->nyS-1;
			ix0 = 0;
			ixMax = 1;
			xMod1 =  0.0;
			xMod2 =  0.5;
			yMod1 =  1.0;
			yMod2 =  1.0;
			IxN =   1; IyN =  0;
			break;
		case 3: // inner right nodes
			iy0 = 1;
			iyMax = Grid->nyS-1;
			ix0 = Grid->nxS-1;
			ixMax = Grid->nxS;
			xMod1 =  0.0;
			xMod2 = -0.5;
			yMod1 =  1.0;
			yMod2 =  1.0;
			IxN =   -1; IyN =  0;
			break;
		case 4: // upper left corner
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0  = 0;
			ixMax = 1;
			xMod1 =  0.0;
			xMod2 =  0.5;
			yMod1 =  0.0;
			yMod2 = -0.5;
			IxN =   1; IyN =  0;
			break;
		case 5: // upper right corner
			iy0 = Grid->nyS-1;
			iyMax = Grid->nyS;
			ix0  = Grid->nxS-1;
			ixMax = Grid->nxS;
			xMod1 =  0.0;
			xMod2 = -0.5;
			yMod1 =  0.0;
			yMod2 = -0.5;
			IxN =   -1; IyN =  0;
			break;
		case 6: // lower right corner
			iy0 = 0;
			iyMax = 1;
			ix0  = Grid->nxS-1;
			ixMax = Grid->nxS;
			xMod1 =  0.0;
			xMod2 = -0.5;
			yMod1 =  0.0;
			yMod2 =  0.5;
			IxN =   -1; IyN =  0;
			break;
		case 7: // lower left corner
			iy0 = 0;
			iyMax = 1;
			ix0  = 0;
			ixMax = 1;
			xMod1 =  0.0;
			xMod2 =  0.5;
			yMod1 =  0.0;
			yMod2 =  0.5;
			IxN =   1; IyN =  0;
			break;
		}
		//#pragma omp parallel for private(iy, ix, iNode, thisParticle, numPart, i, minDist, x, y, iNodeNeigh, neighParticle, dist, closestParticle) schedule(static,32)
		for (iy = iy0; iy < iyMax; ++iy) {
			// Compute DispBound, for proper assignment of passive to the injected particles

			for (ix = ix0; ix < ixMax; ++ix) {
				iNode = ix + iy*Grid->nxS;

				if (Grid->isFixed) {
					if (Particles->passiveGeom==PartPassive_Grid || Particles->passiveGeom==PartPassive_Grid_w_Layers) {
						if (iBlock == 2 || iBlock == 4 || iBlock == 7) { // inner left nodes
							Vx = 0.5* (Physics->Vx[ix + (iy)*Grid->nxVx] + Physics->Vx[ix + (iy+1)*Grid->nxVx]);
							if (Vx>1e-8) {
								inject = true;
							} else {
								inject = false;
							}
							Particles->dispAtBoundL[iy] += Vx * Physics->dtAdv;
							if (Particles->dispAtBoundL[iy]>Particles->passiveDx) {
								Particles->dispAtBoundL[iy] -= Particles->passiveDx;
								if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
									if (Particles->currentPassiveAtBoundL[iy]<2) {
										Particles->currentPassiveAtBoundL[iy] = abs(Particles->currentPassiveAtBoundL[iy]-1); // i.e. if 1->0, if 0->1
									} else {
										Particles->currentPassiveAtBoundL[iy] = abs(Particles->currentPassiveAtBoundL[iy]-2-1)+2; // i.e. if 1->0, if 0->1
									}
								} else {
									Particles->currentPassiveAtBoundL[iy] = abs(Particles->currentPassiveAtBoundL[iy]-1); // i.e. if 1->0, if 0->1
								}
							}
							forcePassive = true;
							passive = Particles->currentPassiveAtBoundL[iy];
						} else if (iBlock == 3 || iBlock == 5 || iBlock == 6) { // inner right nodes
							if (Vx<-1e-8) {
								inject = true;
							} else {
								inject = false;
							}
							Vx = 0.5*(Physics->Vx[ix + (iy)*Grid->nxVx] + Physics->Vx[ix + (iy+1)*Grid->nxVx]);
							Particles->dispAtBoundR[iy] -= Vx * Physics->dtAdv;
							if (Particles->dispAtBoundR[iy]>Particles->passiveDx) {
								Particles->dispAtBoundR[iy] -= Particles->passiveDx;
								if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
									if (Particles->currentPassiveAtBoundR[iy]<2) {
										Particles->currentPassiveAtBoundR[iy] = abs(Particles->currentPassiveAtBoundR[iy]-1); // i.e. if 1->0, if 0->1
									} else {
										Particles->currentPassiveAtBoundR[iy] = abs(Particles->currentPassiveAtBoundR[iy]-2-1)+2; // i.e. if 1->0, if 0->1
									}
								} else {
									Particles->currentPassiveAtBoundR[iy] = abs(Particles->currentPassiveAtBoundR[iy]-1); // i.e. if 1->0, if 0->1
								}
							}
							forcePassive = true;

							passive = Particles->currentPassiveAtBoundR[iy] ;
						} else {
							forcePassive = false;
							inject = false;
						}
					}
				}


				if (inject) {
					//printf("koko\n");

					// free all Particles from this node

					numPart = 0.;
					thisParticle = Particles->linkHead[iNode];
					while (thisParticle != NULL && numPart<minNumPart) {
						thisParticle = thisParticle->next;
						numPart += 1.;
					}
					//printf("numPart = %.2e, minNumPart = %.2e\n", numPart, minNumPart);
					if (Method == 0) {
						//printf("numPart = %.2e, minNumPart = %.2e\n", numPart, minNumPart);
						//if (numPart<minNumPart) {
							while (Particles->linkHead[iNode] != NULL)
							{
								temp = Particles->linkHead[iNode];
								Particles->linkHead[iNode] = Particles->linkHead[iNode]->next;
								free(temp);
								PartAdded[iNode] -= 1;
							}
							//printf("iBlock = %i, A PartAdded[iNode] = %i\n", iBlock, PartAdded[iNode]);
							// copy the neighbour node
							iNodeNeigh = ix+IxN + (iy+IyN)*Grid->nxS;
							//printf("IxN = %i, IyN = %i\n",IxN, IyN);
							neighParticle = Particles->linkHead[iNodeNeigh] ;
							while (neighParticle != NULL) {
								compute xShiftFac = (compute)(IxN);
								compute yShiftFac = (compute)(IyN);
								x = neighParticle->x -xShiftFac*Grid->dx;
								y = neighParticle->y -yShiftFac*Grid->dy;

								if (x>Grid->xmin && x<Grid->xmax) {
									if (y>Grid->ymin && y<Grid->ymax) {
										addSingleParticle(&Particles->linkHead[iNode], neighParticle);
										//printf("x = %.2e, 2.0*(-0.5 + (rand()  1000)/1000.0) * 0.0001*Grid->dx = %.2e\n", x, 2.0*(-0.5 + (rand() % 1000)/1000.0) * 0.001*Grid->dx);
										Particles->linkHead[iNode]->x = x;// + 2.0*(-0.5 + (rand() % 1000)/1000.0) * 0.00001*Grid->dx; // +- 1% of dx
										Particles->linkHead[iNode]->y = y;
#if (STORE_PARTICLE_POS_INI)
										Particles->linkHead[iNode]->xIni = x - Vx*Physics->time;
										Particles->linkHead[iNode]->yIni = y;
#endif

										// Wipe out the stress history (not clear that it's a good idea, but for the moment, not wiping it causes instability so...)
										//Particles->linkHead[iNode]->sigma_xx_0 *= .9;
										Particles->linkHead[iNode]->sigma_xy_0 *= 0.0;
#if (DARCY)
										//Particles->linkHead[iNode]->DeltaP0 *= .9;
										//Particles->linkHead[iNode]->phi = Particles->linkHead[iNode]->phi + 0.5*(MatProps->phiIni[Particles->linkHead[iNode]->phase]-Particles->linkHead[iNode]->phi);// * ( 1.0 + 0.5*(0.5 - (rand() % 1000)/1000.0));
#endif



										PartAdded[iNode] += 1;
										Particles->linkHead[iNode]->nodeId = iNode;
										if (forcePassive) {
											Particles->linkHead[iNode]->passive = passive;
										}

									}
								}

								neighParticle = neighParticle->next;

							}
						//}
						//printf("asoko\n");


					} else if (Method == 1) {

						//if (numPart<minNumPart) {
						//printf("************* A particle is about to be injected!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ****************\n");
						int iDum;
						for (iDum = 0; iDum < 2; ++iDum) {


							x = Grid->X[ix] + 2.0*(-xMod1*0.5 + xMod2*Particles->noiseFactor*(rand() % 1000)/1000.0) * 0.25*Grid->DXEC[ix] + xMod2*0.01*Grid->DXEC[ix];
							y = Grid->Y[iy] + 2.0*(-yMod1*0.5 + yMod2*(rand() % 1000)/1000.0) * 0.25*Grid->DYEC[iy] + yMod2*0.01*Grid->DYEC[iy];


							minDist = (Grid->xmax-Grid->xmin)*(Grid->xmax-Grid->xmin);
							for (i=0;i<1;i++) {
								iNodeNeigh = ix+IxN + (iy+IyN)*Grid->nxS;
								//printf("iNode = %i, iNodeNeigh = %i\n",iNode, iNodeNeigh);
								neighParticle = Particles->linkHead[iNodeNeigh];
								while (neighParticle != NULL) {
									dist = (neighParticle->x - x)*(neighParticle->x - x) + (neighParticle->y - y)*(neighParticle->y - y);
									//printf("dist/dx = %.2e, neighParticle->phase = %i\n",dist/Grid->dx, neighParticle->phase);
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

							if (isnan(x)!=0 || isnan(y)!=0 ) {
								printf("nan coord: x = %.2e, y = %.2e\n", x, y);
							}

							/*
							// Wipe out the stress history (not clear that it's a good idea, but for the moment, not wiping it causes instability so...)
							Particles->linkHead[iNode]->sigma_xx_0 *= 0.0;
							Particles->linkHead[iNode]->sigma_xy_0 *= 0.0;
	#if (DARCY)
							Particles->linkHead[iNode]->DeltaP0 *= 0.0;
							//Particles->linkHead[iNode]->phi = MatProps->phiIni[Particles->linkHead[iNode]->phase];
	#endif
							*/

							PartAdded[iNode] += 1;
							if (forcePassive) {
								//printf("A passN = %.2e, passive = %.2e\n",Particles->linkHead[iNode]->passive,passive);
								Particles->linkHead[iNode]->passive = passive;
								//printf("B passN = %.2e, passive = %.2e\n",Particles->linkHead[iNode]->passive,passive);
							}

						}
						//printf("injP.phi = %.2e, injP.DeltaP0 = %.2e\n",Particles->linkHead[iNode]->phi, Particles->linkHead[iNode]->DeltaP0);

						//}

					}
				}
			}




		}
	}
	//printf("Out\n");


	for (i = 0; i < Grid->nSTot; ++i) {
		Particles->n += PartAdded[i];
	}

	//printf("\n\n\nParticles->n = %i\n≠\n\n", Particles->n);

	free(PartAdded);
}































# if (ADVECT_METHOD == 0)

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
#if (ADVECT_VEL_AND_VISCOSITY)
	compute* VxGrid = (compute*) malloc(4*Grid->nVxTot*sizeof(compute));
	compute* VyGrid = (compute*) malloc(4*Grid->nVyTot*sizeof(compute));

	compute* sumOfWeights_Vx = (compute*) malloc(4* Grid->nVxTot*sizeof(compute));
	compute* sumOfWeights_Vy = (compute*) malloc(4* Grid->nVyTot*sizeof(compute));

	compute* sumOfWeights_EC 	= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));

	compute* ZGrid 		= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));
#if (DARCY)
	compute* ZbGrid 		= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));
#endif

	int iCell;
	int iVx, iVy;
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nVxTot; ++i) {
		VxGrid[i] = 0;
		sumOfWeights_Vx[i] = 0;
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nVyTot; ++i) {
		VyGrid[i] = 0;
		sumOfWeights_Vy[i] = 0;
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nECTot; ++i) {
		sumOfWeights_EC[i] = 0;
		ZGrid[i] = 0;
#if (DARCY)
		ZbGrid[i] = 0;
#endif
	}
#endif

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


#if (ADVECT_METHOD==1)
void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	int ix, iy, iCell, iBound;
	int iL, iR, iU, iD;
	printf("A\n");
	compute* VxCell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* VyCell = (compute*) malloc(Grid->nECTot * sizeof(compute));

	compute* Vx0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Vy0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));

	compute* dVxCell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* dVyCell = (compute*) malloc(Grid->nECTot * sizeof(compute));

/*
#if (ADVECT_VEL_AND_VISCOSITY)
	int i;
#if (VEL_VISC_METHOD == 0)
	compute* VxGrid = (compute*) malloc(4*Grid->nVxTot*sizeof(compute));
	compute* VyGrid = (compute*) malloc(4*Grid->nVyTot*sizeof(compute));

	compute* sumOfWeights_Vx = (compute*) malloc(4* Grid->nVxTot*sizeof(compute));
	compute* sumOfWeights_Vy = (compute*) malloc(4* Grid->nVyTot*sizeof(compute));
#else
	compute* VxCell2 = (compute*) malloc(4*Grid->nECTot * sizeof(compute));
	compute* VyCell2 = (compute*) malloc(4*Grid->nECTot * sizeof(compute));
#endif
	compute* PCell   = (compute*) malloc(4*Grid->nECTot * sizeof(compute));
	compute P;
	compute* sumOfWeights_EC 	= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));

	compute* ZGrid 		= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));
#if (DARCY)
	compute* ZbGrid 		= (compute*) malloc(4 * Grid->nECTot * sizeof(compute));
#endif
	int iVx, iVy;
#if (VEL_VISC_METHOD == 0)
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nVxTot; ++i) {
		VxGrid[i] = 0;
		sumOfWeights_Vx[i] = 0;
	}
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nVyTot; ++i) {
		VyGrid[i] = 0;
		sumOfWeights_Vy[i] = 0;
	}
#endif
#pragma omp parallel for private(i) schedule(static,32)
	for (i = 0; i < 4*Grid->nECTot; ++i) {
		sumOfWeights_EC[i] = 0;
#if (VEL_VISC_METHOD == 1)
		VxCell2[i] = 0.0;
		VyCell2[i] = 0.0;
#endif
		PCell[i]   = 0.0;
		ZGrid[i] = 0;
#if (DARCY)
		ZbGrid[i] = 0;
#endif
	}
#endif
*/
	printf("B\n");
	// interp Vx on cell centers
	// =================================================

	// Loop over cells except first and last column
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix   + iy*Grid->nxEC;
			iR 		= ix   + iy*Grid->nxVx;
			iL 		= ix-1 + iy*Grid->nxVx;
			VxCell[iCell] = (Physics->Vx[iR] + Physics->Vx[iL])/2.0;
			Vx0Cell[iCell] = (Physics->Vx0[iR] + Physics->Vx0[iL])/2.0;
		}
	}
	// Loop over first and last column
	int ixCell, ixNeighCell, ixVx;
	for (iBound = 0; iBound < 2; ++iBound) {
		if (iBound == 0) {
			ixCell = 0;
			ixNeighCell = 1;
			ixVx = 0;
		} else {
			ixCell = Grid->nxEC-1;
			ixNeighCell = Grid->nxEC-2;
			ixVx = Grid->nxVx-1;
		}
		for (iy = 0; iy < Grid->nyEC; ++iy) {
			iCell = ixCell + iy*Grid->nxEC;
			VxCell[iCell] = 2.0*Physics->Vx[ixVx + iy*Grid->nxVx] - VxCell[ixNeighCell + iy*Grid->nxEC]; // i.e ix-0 at the left boundary; ix-1 at the right
			Vx0Cell[iCell] = 2.0*Physics->Vx0[ixVx + iy*Grid->nxVx] - Vx0Cell[ixNeighCell + iy*Grid->nxEC]; // i.e ix-0 at the left boundary; ix-1 at the right
		}
	}

	// interp Vy on cell centers
	// =================================================

	// Loop over cells except first and last row
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell 	= ix   +  iy   *Grid->nxEC;
			iU 		= ix   +  iy   *Grid->nxVy;
			iD 		= ix   + (iy-1)*Grid->nxVy;
			VyCell[iCell] = (Physics->Vy[iU] + Physics->Vy[iD])/2.0;
			Vy0Cell[iCell] = (Physics->Vy0[iU] + Physics->Vy0[iD])/2.0;
		}
	}
	// Loop over first and last row
	int iyCell, iyNeighCell, iyVy;
	for (iBound = 0; iBound < 2; ++iBound) {
		if (iBound == 0) {
			iyCell = 0;
			iyNeighCell = 1;
			iyVy = 0;
		} else {
			iyCell = Grid->nyEC-1;
			iyNeighCell = Grid->nyEC-2;
			iyVy = Grid->nyVy-1;
		}
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell = ix + iyCell*Grid->nxEC;
			VyCell[iCell] = 2.0*Physics->Vy[ix + iyVy*Grid->nxVy] - VyCell[ix + iyNeighCell*Grid->nxEC]; // i.e ix-0 at the left boundary; ix-1 at the right
			Vy0Cell[iCell] = 2.0*Physics->Vy0[ix + iyVy*Grid->nxVy] - Vy0Cell[ix + iyNeighCell*Grid->nxEC]; // i.e ix-0 at the left boundary; ix-1 at the right
			//VyCell[iCell] = Physics->Vy[ix + iyVy*Grid->nxVy]; // i.e ix-0 at the left boundary; ix-1 at the right
		}
	}


	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		dVxCell[iCell] = VxCell[iCell] - Vx0Cell[iCell];
		dVyCell[iCell] = VyCell[iCell] - Vy0Cell[iCell];
		//printf("dVyCell = %.2e, VyCell = %.2e, Vy0Cell = %.2e\n", dVxCell[iCell], VxCell[iCell], Vx0Cell[iCell] );
	}


	INIT_PARTICLE


	compute locX, locY;

	// Loop through nodes
//#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY) schedule(static,32)
	for (iy = 0; iy < Grid->nyS; ++iy) {
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];

			while (thisParticle!=NULL) {

				locX = thisParticle->x-Grid->X[ix];
				locY = thisParticle->y-Grid->Y[iy];

				if (locX<0) {
					locX = 2.0*(locX/Grid->DXS[ix-1]);
				} else {
					locX = 2.0*(locX/Grid->DXS[ix]);
				}
				if (locY<0) {
					locY = 2.0*(locY/Grid->DYS[iy-1]);
				} else {
					locY = 2.0*(locY/Grid->DYS[iy]);
				}

				locX = locX/2.0;
				locY = locY/2.0;
				compute locX2 = fabs(locX);
				compute locY2 = fabs(locY);
				compute Vx, Vy;
				if 		  (locX>0.0 && locY>0.0) {
					Vx = (1.0-locX2) *  Physics->Vx[ix   + (iy+1) *Grid->nxVx]  + locX2 * Physics->Vx[ix+1 + (iy+1) *Grid->nxVx]  ;
					Vy = (1.0-locY2) *  Physics->Vy[ix+1 + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy[ix+1 + (iy+1) *Grid->nxVy]  ;
				} else if (locX>0.0 && locY<=0.0) {
					Vx = (1.0-locX2) *  Physics->Vx[ix   + (iy  ) *Grid->nxVx]  + locX2 * Physics->Vx[ix+1 + (iy  ) *Grid->nxVx]  ;
					Vy = (1.0-locY2) *  Physics->Vy[ix+1 + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy[ix+1 + (iy-1) *Grid->nxVy]  ;
				} else if (locX<=0.0 && locY>0.0) {
					Vx = (1.0-locX2) *  Physics->Vx[ix   + (iy+1) *Grid->nxVx]  + locX2 * Physics->Vx[ix-1 + (iy+1) *Grid->nxVx]  ;
					Vy = (1.0-locY2) *  Physics->Vy[ix   + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy[ix   + (iy+1) *Grid->nxVy]  ;
				} else if (locX<=0.0 && locY<=0.0) {
					Vx = (1.0-locX2) *  Physics->Vx[ix   + (iy  ) *Grid->nxVx]  + locX2 * Physics->Vx[ix-1 + (iy  ) *Grid->nxVx]  ;
					Vy = (1.0-locY2) *  Physics->Vy[ix   + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy[ix   + (iy-1) *Grid->nxVy]  ;
				}
				compute Vx0, Vy0;
				if 		  (locX>0.0 && locY>0.0) {
					Vx0 = (1.0-locX2) *  Physics->Vx0[ix   + (iy+1) *Grid->nxVx]  + locX2 * Physics->Vx0[ix+1 + (iy+1) *Grid->nxVx]  ;
					Vy0 = (1.0-locY2) *  Physics->Vy0[ix+1 + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy0[ix+1 + (iy+1) *Grid->nxVy]  ;
				} else if (locX>0.0 && locY<=0.0) {
					Vx0 = (1.0-locX2) *  Physics->Vx0[ix   + (iy  ) *Grid->nxVx]  + locX2 * Physics->Vx0[ix+1 + (iy  ) *Grid->nxVx]  ;
					Vy0 = (1.0-locY2) *  Physics->Vy0[ix+1 + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy0[ix+1 + (iy-1) *Grid->nxVy]  ;
				} else if (locX<=0.0 && locY>0.0) {
					Vx0 = (1.0-locX2) *  Physics->Vx0[ix   + (iy+1) *Grid->nxVx]  + locX2 * Physics->Vx0[ix-1 + (iy+1) *Grid->nxVx]  ;
					Vy0 = (1.0-locY2) *  Physics->Vy0[ix   + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy0[ix   + (iy+1) *Grid->nxVy]  ;
				} else if (locX<=0.0 && locY<=0.0) {
					Vx0 = (1.0-locX2) *  Physics->Vx0[ix   + (iy  ) *Grid->nxVx]  + locX2 * Physics->Vx0[ix-1 + (iy  ) *Grid->nxVx]  ;
					Vy0 = (1.0-locY2) *  Physics->Vy0[ix   + (iy  ) *Grid->nxVy]  + locY2 * Physics->Vy0[ix   + (iy-1) *Grid->nxVy]  ;
				}


				/*
				compute dVx, dVy;
				dVx = ( .25*(1.0-locX)*(1.0-locY)*dVxCell[ix  +(iy  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*dVxCell[ix  +(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*dVxCell[ix+1+(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*dVxCell[ix+1+(iy  )*Grid->nxEC] )  ;

				thisParticle->Vx += dVx;



				dVy = ( .25*(1.0-locX)*(1.0-locY)*dVyCell[ix  +(iy  )*Grid->nxEC]
					  + .25*(1.0-locX)*(1.0+locY)*dVyCell[ix  +(iy+1)*Grid->nxEC]
					  + .25*(1.0+locX)*(1.0+locY)*dVyCell[ix+1+(iy+1)*Grid->nxEC]
					  + .25*(1.0+locX)*(1.0-locY)*dVyCell[ix+1+(iy  )*Grid->nxEC] )  ;


				thisParticle->Vy += dVy;
				*/



				/*
				Vx = ( .25*(1.0-locX)*(1.0-locY)*VxCell[ix  +(iy  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VxCell[ix  +(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VxCell[ix+1+(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VxCell[ix+1+(iy  )*Grid->nxEC] )  ;



				Vy = ( .25*(1.0-locX)*(1.0-locY)*VyCell[ix  +(iy  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VyCell[ix  +(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VyCell[ix+1+(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VyCell[ix+1+(iy  )*Grid->nxEC] )  ;
					 */

				//thisParticle->Vx += dVx;
				//thisParticle->Vy += dVy;

				//thisParticle->Vx += Vx-Vx0;
				//thisParticle->Vy += Vy-Vy0;

				thisParticle->Vx = Vx;
				thisParticle->Vy = Vy;
				thisParticle->x += thisParticle->Vx  * Physics->dtAdv;
				thisParticle->y += thisParticle->Vy  * Physics->dtAdv;




				/*
				int IX, IY;
				IX = round((thisParticle->x - Grid->xmin)/Grid->dx);
				IY = round((thisParticle->y - Grid->ymin)/Grid->dy);

				locX = thisParticle->x-Grid->X[IX];
				locY = thisParticle->y-Grid->Y[IY];

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



				//compute Vx, Vy;
				Vx = ( .25*(1.0-locX)*(1.0-locY)*VxCell[IX  +(IY  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VxCell[IX  +(IY+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VxCell[IX+1+(IY+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VxCell[IX+1+(IY  )*Grid->nxEC] )  ;



				Vy = ( .25*(1.0-locX)*(1.0-locY)*VyCell[IX  +(IY  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VyCell[IX  +(IY+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VyCell[IX+1+(IY+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VyCell[IX+1+(IY  )*Grid->nxEC] )  ;


				thisParticle->Vx = Vx;
				thisParticle->Vy = Vy;
				*/



				/*
				if (isnan(thisParticle->x)!=0 || isnan(thisParticle->y)!=0 ) {
					printf("adv 2 nan coord: x = %.2e, y = %.2e\n", thisParticle->x, thisParticle->y);
					printf("ix = %i, iy = %i, VyCell = %.2e, %.2e, %.2e, %.2e\n", ix, iy, VyCell[ix  +(iy  )*Grid->nxEC], VyCell[ix  +(iy+1)*Grid->nxEC], VyCell[ix+1+(iy+1)*Grid->nxEC], VyCell[ix+1+(iy  )*Grid->nxEC]);
				}
				*/
				//printf("C\n")  ;

				/*
#if (ADVECT_VEL_AND_VISCOSITY)
				compute Z;
				Z  	= ( .25*(1.0-locX)*(1.0-locY)*Physics->Z[ix  +(iy  )*Grid->nxEC]
				 	  + .25*(1.0-locX)*(1.0+locY)*Physics->Z[ix  +(iy+1)*Grid->nxEC]
				 	  + .25*(1.0+locX)*(1.0+locY)*Physics->Z[ix+1+(iy+1)*Grid->nxEC]
				 	  + .25*(1.0+locX)*(1.0-locY)*Physics->Z[ix+1+(iy  )*Grid->nxEC] )  ;
				P  	= ( .25*(1.0-locX)*(1.0-locY)*Physics->P[ix  +(iy  )*Grid->nxEC]
				 	  + .25*(1.0-locX)*(1.0+locY)*Physics->P[ix  +(iy+1)*Grid->nxEC]
				 	  + .25*(1.0+locX)*(1.0+locY)*Physics->P[ix+1+(iy+1)*Grid->nxEC]
				 	  + .25*(1.0+locX)*(1.0-locY)*Physics->P[ix+1+(iy  )*Grid->nxEC] )  ;
#if (DARCY)
				compute Zb;
				Zb  	= ( .25*(1.0-locX)*(1.0-locY)*Physics->Zb[ix  +(iy  )*Grid->nxEC]
				 	      + .25*(1.0-locX)*(1.0+locY)*Physics->Zb[ix  +(iy+1)*Grid->nxEC]
				 	      + .25*(1.0+locX)*(1.0+locY)*Physics->Zb[ix+1+(iy+1)*Grid->nxEC]
				 	      + .25*(1.0+locX)*(1.0-locY)*Physics->Zb[ix+1+(iy  )*Grid->nxEC] )  ;
#endif
	*/
				/*
				//printf("ix = %i, iy = %i\n", ix, iy);

				// Interpolate velocities from particles back to nodes
				//locX = locX0; // important for using shape functions
				//locY = locY0;
				compute signX, signY;
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

				yModVy[0] =  1.0; // ul
				yModVy[1] =  1.0; // ur
				yModVy[2] = -1.0; // ll
				yModVy[3] = -1.0; //lr

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

				int IxNV[4], IyNV[4];
				IxNV[0] =  0;   IyNV[0] =  1; // ul
				IxNV[1] =  1;	IyNV[1] =  1; // ur
				IxNV[2] =  0; 	IyNV[2] =  0; // ll
				IxNV[3] =  1; 	IyNV[3] =  0; // lr

				compute weight;

#if (VEL_VISC_METHOD == 0)
				for (i=0; i<4; i++) {
					iVx = (ix+IxNV[i]+ixMod + (iy+IyNV[i]) * Grid->nxVx);

					weight = (1.0+signX*xModVx[i]*(fabs(locX)-1.0) )*fabs(locY + yModVx[i]*1.0);

					//printf("locX = %.2f, locY = %.2f, ix = %i, ixN = %i, iy = %i, iyN = %i, weight = %.2f, xContrib = %.2f, yContrib = %.2f\n", locX, locY, ix, ix+IxNV[i]+ixMod, iy, (iy+IyNV[i]) , weight, A, B);
					VxGrid[iVx*4+i] += thisParticle->Vx * weight;
					sumOfWeights_Vx[iVx*4+i] += weight;

				}



				//	printf("Vy\n");

				for (i=0; i<4; i++) {
					iVy = (ix+IxNV[i] + (iy+IyNV[i]+iyMod) * Grid->nxVy);

					weight =    fabs(locX + xModVy[i]*1.0) * (1.0+signY*yModVy[i]*(fabs(locY)-1.0)) ;

					//printf("locX = %.2f, locY = %.2f, ix = %i, ixN = %i, iy = %i, iyN = %i, weight = %.2f, xContrib = %.2f, yContrib = %.2f\n", locX, locY, ix, ix+IxNV[i], iy, (iy+IyNV[i]+iyMod) , weight, A, B);


					VyGrid[iVy*4+i] += thisParticle->Vy * weight;
					sumOfWeights_Vy[iVy*4+i] += weight;

				}
#endif
*/


				//exit(0);

				/*



				// Interpolate etaVisc back to Nodes
				//locX = locX0; // important for using shape functions
				//locY = locY0;
				int IxN[4], IyN[4];
				IxN[0] =  0;  	IyN[0] =  0; // lower left
				IxN[1] =  0;	IyN[1] =  1; // upper left
				IxN[2] =  1; 	IyN[2] =  1; // upper right
				IxN[3] =  1; 	IyN[3] =  0; // lower right

				compute xMod[4], yMod[4];
				xMod[0] = -1; yMod[0] = -1; // ll
				xMod[1] = -1; yMod[1] =  1; // ul
				xMod[2] =  1; yMod[2] =  1; // ur
				xMod[3] =  1; yMod[3] = -1; // lr

				for (i=0; i<4; i++) {
					iCell = (ix+IxN[i] + (iy+IyN[i]) * Grid->nxEC);

					weight = fabs((locX + xMod[i])   *   (locY + yMod[i]));

					sumOfWeights_EC[4*iCell+i] += weight;
#if (VEL_VISC_METHOD == 1)
					VxCell2[4*iCell+i] += thisParticle->Vx*weight;
					VyCell2[4*iCell+i] += thisParticle->Vy*weight;
#endif
					PCell  [4*iCell+i] += P *weight;
					ZGrid[4*iCell+i] += Z*weight;
#if (DARCY)
					ZbGrid[4*iCell+i] += Zb*weight;
#endif
				}



#endif
				//printf("D\n");

				 */

				thisParticle = thisParticle->next;
			}
		}
	}
	printf("E\n");
	compute sum;

#if (ADVECT_VEL_AND_VISCOSITY)
/*
#if (VEL_VISC_METHOD == 0)
	for (iVx = 0; iVx < Grid->nVxTot; ++iVx) {
		sum = sumOfWeights_Vx[4*iVx+0] + sumOfWeights_Vx[4*iVx+1] + sumOfWeights_Vx[4*iVx+2] + sumOfWeights_Vx[4*iVx+3];
		Physics->Vx[iVx] = ( VxGrid[4*iVx+0] + VxGrid[4*iVx+1] + VxGrid[4*iVx+2] + VxGrid[4*iVx+3] ) /sum;
	}

	for (iVy = 0; iVy < Grid->nVyTot; ++iVy) {
		sum = sumOfWeights_Vy[4*iVy+0] + sumOfWeights_Vy[4*iVy+1] + sumOfWeights_Vy[4*iVy+2] + sumOfWeights_Vy[4*iVy+3];
		Physics->Vy[iVy] = ( VyGrid[4*iVy+0] + VyGrid[4*iVy+1] + VyGrid[4*iVy+2] + VyGrid[4*iVy+3] ) /sum;
	}
#endif
	for (iCell = 0; iCell < Grid->nECTot; ++iCell) {
		sum = sumOfWeights_EC[4*iCell+0] + sumOfWeights_EC[4*iCell+1] + sumOfWeights_EC[4*iCell+2] + sumOfWeights_EC[4*iCell+3];
		Physics->Z[iCell] = ( ZGrid[4*iCell+0] + ZGrid[4*iCell+1] + ZGrid[4*iCell+2] + ZGrid[4*iCell+3] ) / sum;
#if (VEL_VISC_METHOD == 1)
		VxCell[iCell] 		= ( VxCell2[4*iCell+0] + VxCell2[4*iCell+1] + VxCell2[4*iCell+2] + VxCell2[4*iCell+3] ) /sum;
		VyCell[iCell] 		= ( VyCell2[4*iCell+0] + VyCell2[4*iCell+1] + VyCell2[4*iCell+2] + VyCell2[4*iCell+3] ) /sum;
#endif
		Physics->P[iCell] 	= ( PCell  [4*iCell+0] + PCell  [4*iCell+1] + PCell  [4*iCell+2] + PCell  [4*iCell+3] ) /sum;
#if (DARCY)
		Physics->Zb[iCell] = ( ZbGrid[4*iCell+0] + ZbGrid[4*iCell+1] + ZbGrid[4*iCell+2] + ZbGrid[4*iCell+3] ) / sum;
#endif

	}
	*/

/*
#if (VEL_VISC_METHOD == 1)
	compute VxBef;
	for (iy = 0; iy < Grid->nyVx; ++iy) {
		for (ix = 0; ix < Grid->nxVx; ++ix) {
			iVx = ix + iy*Grid->nxVx;
			VxBef = Physics->Vx[iVx];
			Physics->Vx[iVx] = (VxCell[ix   + (iy  )*Grid->nxEC] + VxCell[ix+1 + (iy  )*Grid->nxEC])/2.0;

			//printf("DeltaVx = %.2e\n", (Physics->Vx[iVx]-VxBef));
		}
	}

	for (iy = 0; iy < Grid->nyVy; ++iy) {
		for (ix = 0; ix < Grid->nxVy; ++ix) {
			iVy = ix + iy*Grid->nxVy;
			Physics->Vy[iVy] = (VyCell[ix   + (iy  )*Grid->nxEC] + VyCell[ix   + (iy+1)*Grid->nxEC])/2.0;
		}
	}
#endif

*/

/*
#if (VEL_VISC_METHOD == 0)
	free(VxGrid);
	free(VyGrid);
	free(sumOfWeights_Vx);
	free(sumOfWeights_Vy);
#else
	free(VxCell2);
	free(VyCell2);
#endif
	free(PCell);

*/
	free(VxCell);
	free(VyCell);

	free(Vx0Cell);
	free(Vy0Cell);

	free(dVxCell);
	free(dVyCell);


/*
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
*/
#endif



}
#endif




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
		thisParticle->y = Grid->ymin+0.05*Grid->DYEC[0];
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.05*Grid->DYEC[Grid->nyS-1];
	}
	END_PARTICLES


}





void Particles_switchStickyAir(Particles* Particles, Grid* Grid, Physics* Physics, Numerics* Numerics, MatProps* MatProps, BC* BCStokes) {

	int iy, ix, iNode;
	SingleParticle* thisParticle = NULL;
	int phase;

	int iyTop, iyBottom;
	iyTop = floor((Numerics->stickyAirSwitchingDepth - Grid->ymin) / Grid->dy)   + 1;
	iyBottom = iyTop - ceil(Numerics->CFL_fac_Stokes)      - 1 ;

	if(iyTop>Grid->nyS-1) {
		iyTop = Grid->nyS-1;
	}

	if(iyBottom<0) {
		iyBottom = 0;
	}

	int PlusOrMinusOne;
	Numerics->stickyAirTimeSinceLastPassiveSwitch += Physics->dtAdv;
	if (Numerics->stickyAirTimeSinceLastPassiveSwitch>Numerics->stickyAirTimeSwitchPassive) {
		PlusOrMinusOne = (Numerics->stickyAirSwitchPassiveTo%2)*2 - 1;
		Numerics->stickyAirSwitchPassiveTo += PlusOrMinusOne;
	}

	//#pragma omp parallel for private(iy, ix, iNode, thisParticle) schedule(static,16)
	//printf("instickyAir loop\n");
	for (iy = iyBottom; iy < iyTop; ++iy) {
		//	printf("iy = %i\n", iy);
		for (ix = 0; ix < Grid->nxS; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			// Loop through the particles in the cell
			// ======================================
			while (thisParticle!=NULL) {
				if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
					thisParticle->phase = Numerics->stickyAirSwitchPhaseTo;
					thisParticle->passive = Numerics->stickyAirSwitchPassiveTo;
#if (DARCY)
					thisParticle->phi 	= MatProps->phiIni[thisParticle->phase];
#endif
				}
				thisParticle = thisParticle->next;
			}
		}
	}
	//printf("out of stickyAir loop\n");


	if (BCStokes->SetupType == Stokes_Sandbox && BCStokes->Sandbox_NoSlipWall == true) {
		ix = Grid->nxS-1;
		int iNodeUL;
		int iPhase;
		int contribPhase[MatProps->nPhase];
		int maxContrib;
		int majorPhase;
		for (iy = 0; iy < Grid->nyS-1; ++iy) {
			//	printf("iy = %i\n", iy);

			iNode = ix  + (iy  )*Grid->nxS;
			iNodeUL = ix-1  + (iy+1)*Grid->nxS;

			// Check upper left node (with resect to this node)
			// If the major phase is not air then switch air particles in this cell

			// Reinitialize contribs
			// ===================
			for (iPhase=0;iPhase<MatProps->nPhase;++iPhase) {
				contribPhase[iPhase] = 0;
			}


			// Count contribs
			// ===================
			thisParticle = Particles->linkHead[iNodeUL];
			while (thisParticle != NULL) {
				++contribPhase[thisParticle->phase];
				thisParticle = thisParticle->next;
			}

			// Find the most prominent phase
			// ===================
			maxContrib = 0;
			for (iPhase=0;iPhase<MatProps->nPhase;++iPhase) {
				if (contribPhase[iPhase] > maxContrib) {
					majorPhase = iPhase;
					maxContrib = contribPhase[iPhase];
				}
			}
			//printf("majorPhase = %i, ix = %i\n",majorPhase, ix);



			// If the major phase in the upper left node is not air or water, then switch air particles in this node to the major phase of the UL node
			if (majorPhase != Physics->phaseAir && majorPhase != Physics->phaseWater) {
				thisParticle = Particles->linkHead[iNode];
				// Loop through the particles in the cell
				// ======================================
				while (thisParticle!=NULL) {
					if (thisParticle->phase == Physics->phaseAir || thisParticle->phase == Physics->phaseWater) {
						thisParticle->phase = majorPhase;// Numerics->stickyAirSwitchPhaseTo;
						//thisParticle->passive = Numerics->stickyAirSwitchPassiveTo;
#if (DARCY)
						thisParticle->phi 	= MatProps->phiIni[majorPhase];
#endif
					}
					thisParticle = thisParticle->next;
				}
			}

		}
	}

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


	thisParticle->sigma_xx_0 = modelParticle->sigma_xx_0;
	thisParticle->sigma_xy_0 = modelParticle->sigma_xy_0;
#if (CRANK_NICHOLSON_VEL || INERTIA)
	thisParticle->Vx = modelParticle->Vx;
	thisParticle->Vy = modelParticle->Vy;
#endif

#if (STRAIN_SOFTENING)
	thisParticle->strain = modelParticle->strain;
#endif

#if (STORE_PARTICLE_POS_INI)
	thisParticle->xIni = modelParticle->xIni;
	thisParticle->yIni = modelParticle->yIni;
#endif

#if (HEAT)
	thisParticle->T = modelParticle->T;
#endif
#if (DARCY)
	thisParticle->DeltaP0 = modelParticle->DeltaP0;
	thisParticle->phi = modelParticle->phi;
	//	thisParticle->faulted = modelParticle->faulted;
#endif


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


void findNodeForThisParticle(SingleParticle* thisParticle, Grid* Grid)
{
	int ix = (int) round((thisParticle->x-Grid->xmin)/Grid->DXEC[0]);
	int iy = (int) round((thisParticle->y-Grid->ymin)/Grid->DYEC[0]);
	thisParticle->nodeId = ix + iy*Grid->nxS;
	//thisParticle->nodeId = (int) round((thisParticle->x-Grid->xmin)/Grid->dx) + round((thisParticle->y-Grid->ymin)/Grid->dy) * Grid->nxS;
	//printf("DXEC[0] = %.2e, DXS[0] = %.2e, Grid->dx = %.2e, DYEC[0] = %.2e, DYS[0] = %.2e, Grid->dy = %.2e\n",Grid->DXEC[0], Grid->DXS[0], Grid->dx, Grid->DYEC[0], Grid->DYS[0], Grid->dy);

	// pseudo code for the variable grid spacing
	/*
	ix = some_initial_ix_from_input;
	OK = false;
	while (!OK) {
		locX0 = xPart-xNode[ix];
		if (x<xN - dxS(ix)/2) {
			ix -= 1;
		} else if (x>xN + dxS(ix+1)/2) {
			ix+=1;
		} else {
			OK = true;
		}
	}
	 */





}




