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


void Particles_Memory_allocate(Particles* Particles, Grid* Grid) {
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

void Particles_Memory_free(Particles* Particles, Grid* Grid) {
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

	modelParticle.sigma_xx_0 = 0.0;
	modelParticle.sigma_xy_0 = 0.0;

	modelParticle.Dsigma_xx_0 = 0.0;
	modelParticle.Dsigma_xy_0 = 0.0;

	modelParticle.phase = 0;
	modelParticle.passive = 1;
#if (INERTIA || CRANK_NICHOLSON_VEL)
	modelParticle.Vx = 0.0;
	modelParticle.Vy = 0.0;
#endif
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

		}

		Particles->currentPassiveAtBoundR[Grid->nyS-1] = Particles->currentPassiveAtBoundR[Grid->nyS-2]; // overwrite the uppermost one to avoid a disgracious passive switch just at the corner
		Particles->currentPassiveAtBoundL[Grid->nyS-1] = Particles->currentPassiveAtBoundL[Grid->nyS-2];
		Particles->currentPassiveAtBoundR[0] = Particles->currentPassiveAtBoundR[1]; // overwrite the uppermost one to avoid a disgracious passive switch just at the corner
		Particles->currentPassiveAtBoundL[0] = Particles->currentPassiveAtBoundL[1];



		INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, dum, passive) OMP_SCHEDULE
		FOR_PARTICLES
		dum = (int)((thisParticle->x-Grid->xmin)/DX);

		passive = dum%2;
		dum = (int)((thisParticle->y-Grid->ymin)/DY);
		passive += (dum)%2;
		if (passive==1) {
			thisParticle->passive = 0;
		} else {
			thisParticle->passive = 1;
		}

		END_PARTICLES
	}




	if (Particles->passiveGeom==PartPassive_Grid_w_Layers) {
		int dum, passive;

		INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle, dum, passive) OMP_SCHEDULE
		FOR_PARTICLES

		dum = (int)((thisParticle->y-Grid->ymin)/DY);
		passive = (dum)%2;
		if (passive==1) {
			thisParticle->passive += 2;
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



	INIT_PARTICLE
	//#pragma omp parallel for private(iNode, thisParticle, change, ix, iy, loopingParticle, ParticleCounter, locX, locY, Min, Imin, i) OMP_SCHEDULE

	FOR_PARTICLES
	if (thisParticle->x<Grid->xmin) {
		thisParticle->x = Grid->xmin+0.05*Grid->DXEC[0];
	} else if (thisParticle->x>Grid->xmax) {
		thisParticle->x = Grid->xmax-0.05*Grid->DXEC[Grid->nxS-1];
	}

	if (thisParticle->y<Grid->ymin) {
		thisParticle->y = Grid->ymin+0.05*Grid->DYEC[0];
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.05*Grid->DYEC[Grid->nyS-1];
	}
	END_PARTICLES

}

void Particles_deleteIfOutsideTheDomain(Particles* Particles, Grid* Grid)
{
	SingleParticle* nextParticle = NULL;
	// Due to advection error particles might end up outside the model boundaries. This function teleports them back inside
	SingleParticle* thisParticle = NULL;
	int iNode = 0;

	int justDeleted = 0;
	bool change = false;
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {
		thisParticle = Particles->linkHead[iNode];

		while (thisParticle != NULL) {

			nextParticle = thisParticle->next;

			if (nextParticle!=NULL) {
				if (nextParticle->x<Grid->xmin || nextParticle->x>Grid->xmax
						|| nextParticle->y<Grid->ymin || nextParticle->y>Grid->ymax	) {
					thisParticle->next = nextParticle->next;
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


	// Declare a linked list that contains the id of particles that have change cell
	ParticlePointerList* headIdChanged = (ParticlePointerList*) malloc(sizeof(ParticlePointerList));
	headIdChanged->pointer = NULL;
	headIdChanged->next = NULL;

	int oldNodeId;
	SingleParticle* previousParticle;
	// Update the link list
	// =========================


	int ParticleCounter = 0;
	SingleParticle* thisParticle = NULL;
	int iNode = 0;

	printf("First loop\n");

	//#pragma omp parallel for private(iNode, thisParticle, ParticleCounter, oldNodeId, x, y, ix, iy, previousParticle) OMP_SCHEDULE
	for (iNode = 0; iNode < Grid->nSTot; ++iNode) {



		thisParticle = Particles->linkHead[iNode];
		ParticleCounter = 0;
		while (thisParticle != NULL) {

			ParticleCounter++;


			oldNodeId = thisParticle->nodeId;


			findNodeForThisParticle(thisParticle, Grid);

			// If this particle has changed cell
			if (oldNodeId != thisParticle->nodeId) {
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



	printf("Update info\n");
	// 2. Update info of the new cell, i.e. Add this particle to the head of the link list of the new cell
	// ==============================
	ParticlePointerList* IdChanged = NULL;
	IdChanged = headIdChanged;
	while (IdChanged->pointer!=NULL) {
		thisParticle 	= IdChanged->pointer;

		thisParticle->next = Particles->linkHead[thisParticle->nodeId] ;
		//printf("soko\n");
		Particles->linkHead[thisParticle->nodeId] = thisParticle;
		//printf("asoko\n");
		IdChanged 		= IdChanged->next;
	}
	freeParticlePointerList(headIdChanged);

	printf("... info updated\n");
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
		//#pragma omp parallel for private(iy, ix, iNode, thisParticle, numPart, i, minDist, x, y, iNodeNeigh, neighParticle, dist, closestParticle) OMP_SCHEDULE
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

					minDist = (Grid->xmax-Grid->xmin)*(Grid->xmax-Grid->xmin);

					x = Grid->X[ix] + xMod;
					y = Grid->Y[iy] + yMod;
					//printf("A\n");

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
					if (closestParticle!=NULL) {
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
			}
		}
	}

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
		//#pragma omp parallel for private(iy, ix, iNode, thisParticle, numPart, i, minDist, x, y, iNodeNeigh, neighParticle, dist, closestParticle) OMP_SCHEDULE
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
										//Particles->linkHead[iNode]->sigma_xy_0 *= 0.0;
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
					}
				}
			}
		}
	}


	for (i = 0; i < Grid->nSTot; ++i) {
		Particles->n += PartAdded[i];
	}


	free(PartAdded);
}






























void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	INIT_PARTICLE
	int ix, iy, iCell, iBound;
	int iL, iR, iU, iD;
	printf("A\n");
	compute* VxCell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* VyCell = (compute*) malloc(Grid->nECTot * sizeof(compute));

	compute* Vx0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* Vy0Cell = (compute*) malloc(Grid->nECTot * sizeof(compute));

	compute* dVxCell = (compute*) malloc(Grid->nECTot * sizeof(compute));
	compute* dVyCell = (compute*) malloc(Grid->nECTot * sizeof(compute));


	int IxN[4], IyN[4];
	IxN[0] =  0;  	IyN[0] =  0; // lower left
	IxN[1] =  0;	IyN[1] =  1; // upper left
	IxN[2] =  1; 	IyN[2] =  1; // upper right
	IxN[3] =  1; 	IyN[3] =  0; // lower right


	// interp Vx on cell centers
	// =================================================

	// Loop over cells except first and last column
#pragma omp parallel for private(iy, ix, iCell, iR, iL) OMP_SCHEDULE
	for (iy = 0; iy < Grid->nyEC; ++iy) {
		for (ix = 1; ix < Grid->nxEC-1; ++ix) {
			iCell 	= ix   + iy*Grid->nxEC;
			iR 		= ix   + iy*Grid->nxVx;
			iL 		= ix-1 + iy*Grid->nxVx;
			VxCell[iCell] = (Physics->Vx[iR] + Physics->Vx[iL])/2.0;
			//Vx0Cell[iCell] = (Physics->Vx0[iR] + Physics->Vx0[iL])/2.0;
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
		}
	}

	// interp Vy on cell centers
	// =================================================

	// Loop over cells except first and last row
#pragma omp parallel for private(iy, ix, iCell, iU, iD) OMP_SCHEDULE
	for (iy = 1; iy < Grid->nyEC-1; ++iy) {
		for (ix = 0; ix < Grid->nxEC; ++ix) {
			iCell 	= ix   +  iy   *Grid->nxEC;
			iU 		= ix   +  iy   *Grid->nxVy;
			iD 		= ix   + (iy-1)*Grid->nxVy;
			VyCell[iCell] = (Physics->Vy[iU] + Physics->Vy[iD])/2.0;
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
		}
	}


	// Compute the Alpha array
					// add a condi	ztion with signX signY to avoid recomputing alpha if not necessary
	compute* alphaArray = (compute*) malloc(Grid->nSTot * sizeof(compute));
	compute alpha;
#pragma omp parallel for private(iy, ix, iNode) OMP_SCHEDULE
	for (iy=0; iy<Grid->nyS; iy++) {
		for (ix=0; ix<Grid->nxS; ix++) {
			iNode = ix + iy*Grid->nxS;
			alphaArray[iNode]  = .5*Physics->dtAdv*((Physics->Vy[ix+1 + (iy  )*Grid->nxVy] - Physics->Vy[ix   +(iy  )*Grid->nxVy])/Grid->DXEC[ix]
											 	  - (Physics->Vx[ix   + (iy+1)*Grid->nxVx] - Physics->Vx[ix   +(iy  )*Grid->nxVx])/Grid->DYEC[iy]);
			
		}
	}
				

	




	


	compute Vx, Vy, Vx2, Vy2;
	compute locX, locY;
	int IX, IY;
	compute sigma_xx_temp;
	compute tempx, tempy;
				
	// Loop through nodes
#pragma omp parallel for private(iy, ix, iNode, thisParticle, locX, locY, alpha, Vx, Vy, Vx2, Vy2, IX, IY, sigma_xx_temp, tempx, tempy) OMP_SCHEDULE
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








				// Correction without assuming a small angle
				alpha = Interp_NodeVal_Node2Particle_Local(alphaArray, ix, iy, Grid->nxS, Grid->nyS, locX, locY);				
				sigma_xx_temp = thisParticle->sigma_xx_0*cos(1.0*alpha)*cos(1.0*alpha)  -  thisParticle->sigma_xy_0*sin(2.0*alpha);
				thisParticle->sigma_xy_0 = thisParticle->sigma_xy_0*cos(2.0*alpha)  +  thisParticle->sigma_xx_0*sin(2.0*alpha);
				thisParticle->sigma_xx_0 = sigma_xx_temp;




				// =====================================================


				
				Vx = ( .25*(1.0-locX)*(1.0-locY)*VxCell[ix  +(iy  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VxCell[ix  +(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VxCell[ix+1+(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VxCell[ix+1+(iy  )*Grid->nxEC] )  ;



				Vy = ( .25*(1.0-locX)*(1.0-locY)*VyCell[ix  +(iy  )*Grid->nxEC]
					 + .25*(1.0-locX)*(1.0+locY)*VyCell[ix  +(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0+locY)*VyCell[ix+1+(iy+1)*Grid->nxEC]
					 + .25*(1.0+locX)*(1.0-locY)*VyCell[ix+1+(iy  )*Grid->nxEC] )  ;


				
				tempx = thisParticle->x+Vx*Physics->dtAdv;
				tempy = thisParticle->y+Vy*Physics->dtAdv;


				IX = round((tempx - Grid->xmin)/Grid->dx);
				IY = round((tempy - Grid->ymin)/Grid->dy);

				if (tempx<Grid->xmax && tempy<Grid->ymax && tempx>Grid->xmin && tempy>Grid->ymin) {
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

					Vx2 = ( .25*(1.0-locX)*(1.0-locY)*VxCell[IX  +(IY  )*Grid->nxEC]
						+ .25*(1.0-locX)*(1.0+locY)*VxCell[IX  +(IY+1)*Grid->nxEC]
						+ .25*(1.0+locX)*(1.0+locY)*VxCell[IX+1+(IY+1)*Grid->nxEC]
						+ .25*(1.0+locX)*(1.0-locY)*VxCell[IX+1+(IY  )*Grid->nxEC] )  ;

					Vy2 = ( .25*(1.0-locX)*(1.0-locY)*VyCell[IX  +(IY  )*Grid->nxEC]
						+ .25*(1.0-locX)*(1.0+locY)*VyCell[IX  +(IY+1)*Grid->nxEC]
						+ .25*(1.0+locX)*(1.0+locY)*VyCell[IX+1+(IY+1)*Grid->nxEC]
						+ .25*(1.0+locX)*(1.0-locY)*VyCell[IX+1+(IY  )*Grid->nxEC] )  ;

					Vx = .5*(Vx+Vx2);
					Vy = .5*(Vy+Vy2);

				}



#if (INERTIA || CRANK_NICHOLSON_VEL)
				thisParticle->Vx = Vx;
				thisParticle->Vy = Vy;
				thisParticle->x += thisParticle->Vx  * Physics->dtAdv;
				thisParticle->y += thisParticle->Vy  * Physics->dtAdv;
#else
				thisParticle->x += Vx  * Physics->dtAdv;
				thisParticle->y += Vy  * Physics->dtAdv;
#endif

				thisParticle = thisParticle->next;
			}
		}
	}
	printf("E\n");
	compute sum;

	free(VxCell);
	free(VyCell);

	free(Vx0Cell);
	free(Vy0Cell);

	free(dVxCell);
	free(dVyCell);

	free(alphaArray);


}





void Particles_Periodicize(Particles* Particles, Grid* Grid)
{
	// Make particles do the loop

	//if (BC->VxT < BC->VxB) {
	// sinistral simple shear:
	// particles go out through the left boundary and renter through the right one
	INIT_PARTICLE
#pragma omp parallel for private(iNode, thisParticle) OMP_SCHEDULE
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

	//#pragma omp parallel for private(iy, ix, iNode, thisParticle) OMP_SCHEDULE
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

	thisParticle->Dsigma_xx_0 = modelParticle->Dsigma_xx_0;
	thisParticle->Dsigma_xy_0 = modelParticle->Dsigma_xy_0;

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




