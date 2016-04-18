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

	SingleParticle modelParticle;


	coord xP, yP;

	modelParticle.x = 0;
	modelParticle.y = 0;
	modelParticle.nodeId = 0;
	modelParticle.T = 0;
	modelParticle.sigma_xx_0 = 0;
	modelParticle.sigma_xy_0 = 0;
	modelParticle.phase = 0;
	modelParticle.passive = 0;
	modelParticle.next = NULL;

	// Loop through nodes
	// ==================
	int iNode = 0;
	for(iy=0;iy<Grid->nyC;iy++) {
		for(ix=0;ix<Grid->nxC;ix++) {
			// Get the coordinates of the lower left corner of the shifted cell (i.e. cell centered on the node ix, iy)
			x = Grid->xmin + ix*Grid->dx;// - 0.5*Grid->dx;
			y = Grid->ymin + iy*Grid->dy;// - 0.5*Grid->dy ;



			// Loop through Particles in the cell
			// ==================================
			for (iPy=0;iPy<Particles->nPCY;iPy++) {
				for (iPx=0;iPx<Particles->nPCX;iPx++) {



					// Assign coordinate
					// =================
					//printf("Rand1 = %.4f, Rand2 = %.4f\n",(0.5 - (rand() % 1000)/1000.0),(0.5 - (rand() % 1000)/1000.0));
					modelParticle.x 	= x + 0.5*dxP + iPx*dxP + Particles->noiseFactor*dxP*(0.5 - (rand() % 1000)/1000.0);
					modelParticle.y 	= y + 0.5*dyP + iPy*dyP + Particles->noiseFactor*dyP*(0.5 - (rand() % 1000)/1000.0);

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
void Particles_initPhase(Grid* Grid, Particles* Particles)
{
	int Setup = 5;
	srand(time(NULL));

	if (Setup==0) {
		INIT_PARTICLE
		FOR_PARTICLES
			thisParticle->phase = 0;
		END_PARTICLES

	}
	else if (Setup==1) {

		// Simple inclusion
		int object = 0; // 0 = circle, 1 = square
		int nObjects = 1;
		int i, A;
		coord sqrDistance;
		coord radius = (0.2*(Grid->ymax-Grid->ymin)/2);
		coord rx = (0.4*(Grid->ymax-Grid->ymin)/2);
		coord ry = (0.4*(Grid->ymax-Grid->ymin)/2);
		coord sqrRadius =  radius * radius;
		coord alpha;
		//coord sqrRadius = 0.3*0.3;
		coord cX = Grid->xmin+(Grid->xmax-Grid->xmin)/2;
		coord cY = Grid->ymin + (Grid->ymax-Grid->ymin)*0.5;//Grid->ymin + 0.0*(Grid->ymax-Grid->ymin)/2.0;

		coord x, y, Ex, Ey;

		compute thickCrust = Grid->ymax-Grid->ymin;

		int iObj;

		//INIT_TIMER
		//TIC

		INIT_PARTICLE
		FOR_PARTICLES
			thisParticle->phase = 0;
			if (thisParticle->y<Grid->ymin+thickCrust) {
				thisParticle->phase = 1;
			}
			if (thisParticle->y<0.05*(Grid->ymax-Grid->ymin)) {
				//thisParticle->phase = 2;
			}

		END_PARTICLES


		coord* CX = (coord*) malloc(nObjects*sizeof(coord));
		coord* CY = (coord*) malloc(nObjects*sizeof(coord));

		for (iObj = 0; iObj<nObjects; iObj++) {
			CX[iObj] =  Grid->xmin + (Grid->xmax-Grid->xmin)*(rand() % 1000)/1000.0;
			CY[iObj] =  Grid->ymin + (Grid->ymax-Grid->ymin)*(rand() % 1000)/1000.0;
		}

		if (nObjects == 1) {
			CX[0] = cX;
			CY[0] = cY;
		}




#pragma omp parallel for private(iNode, thisParticle, iObj, x, y ) schedule(static,32)
		FOR_PARTICLES
		for (iObj = 0; iObj<nObjects; iObj++) {
			if (object == 0) {
				x = (thisParticle->x-CX[iObj]);
				y = (thisParticle->y-CY[iObj]);
				//if (sqrDistance < sqrRadius) {
				if ( ( (x*x)/(rx*rx) + (y*y)/(ry*ry) )<= 1 ) { // note: infinity shape condition: sqrDistance<rx*cos(alpha)*ry*sin(alpha)
					if (thisParticle->y<0.5*(Grid->ymax-Grid->ymin))
					thisParticle->phase = 2;
					break;
				}

			}
			else if (object == 1) {
				if (abs(thisParticle->x-cX) < rx && abs(thisParticle->y-cY) <ry) {
					if (thisParticle->y<0.5*(Grid->ymax-Grid->ymin))
					thisParticle->phase = 2;
					break;
				}

			}
		}

		END_PARTICLES

		//TOC
		//printf("looping through the particles take: %.5f s\n", toc);

		free(CX);
		free(CY);

	}


	else if (Setup==2) {
		// Sinusoidal basement
		compute WaveNumber = 3; // Wavelength
		compute phase = 0.5*PI;
		compute Amplitude = 0.00*(Grid->ymax-Grid->ymin);
		compute Thickness = 0.05*(Grid->ymax-Grid->ymin);
		compute x,y;


		INIT_PARTICLE
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
	else if (Setup==3) {
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

		INIT_PARTICLE
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

	else if (Setup==4) {
			// Sandbox
			compute Hmodel = (Grid->ymax-Grid->ymin);
			compute Wmodel = (Grid->xmax-Grid->xmin);
			compute thickCrust = 0.4*(Grid->ymax-Grid->ymin);
			compute thickBase = 0.1*thickCrust;

			compute tanAngleCorner  = tan(10*PI/180);
			compute lengthCorner = (Grid->ymax-Grid->ymin);
			compute xCorner = Grid->xmax-lengthCorner - thickCrust/tanAngleCorner; // position of the triangle at the bottom of the box
			compute xP_from_xCorner;

			int nLayers = 0;
			compute HLayers[2] = {0.1*Hmodel, 0.18*Hmodel};
			compute TLayers[2] = {0.03*Hmodel, 0.03*Hmodel};

			compute Hseamount = 1.5*thickCrust;
			compute xSeamount = Grid->xmax-Wmodel/1.8;
			compute tanAngleSeamount = tan(30*PI/180);

			compute xPfromSeamountCenter;


			int i;
			INIT_PARTICLE
			FOR_PARTICLES
			thisParticle->phase = 0;
			END_PARTICLES


			FOR_PARTICLES

				xP_from_xCorner = (thisParticle->x-xCorner);
				if (xP_from_xCorner>0 && (thisParticle->y < xP_from_xCorner*tanAngleCorner) ) {
					thisParticle->phase = 1;
				}

				if (thisParticle->y<thickCrust) {
					thisParticle->phase = 1;
				}

				if (thisParticle->y<thickBase) {
					thisParticle->phase = 2;
				}


				for (i=0;i<nLayers;++i) {
					if (thisParticle->y>HLayers[i] && thisParticle->y<HLayers[i]+TLayers[i]) {
						thisParticle->phase = 2;
					}
				}

				xPfromSeamountCenter = fabs(thisParticle->x-xSeamount);
				if (-thisParticle->y+Hseamount>xPfromSeamountCenter*tanAngleSeamount) {
					//thisParticle->phase = 3;
				}



				END_PARTICLES
	}

	else if (Setup==5) {
		// Snadbox with sinusoidal basement
		compute WaveNumber = 10; // Wavelength
		compute phase = 1.5*PI;

		compute thickCrust = 0.2*(Grid->ymax-Grid->ymin);

		compute Thickness = 0.2*thickCrust;
		compute Amplitude = Thickness;
		compute x,y;

		compute thickDecollement = 0.05*thickCrust;

		compute tanAngleCorner  = tan(8*PI/180);
		compute lengthCorner = (Grid->ymax-Grid->ymin);
		compute xCorner = Grid->xmax-lengthCorner - 2*thickCrust/tanAngleCorner; // position of the triangle at the bottom of the box
		compute xP_from_xCorner;

		INIT_PARTICLE
		FOR_PARTICLES
		x = thisParticle->x;
		y = thisParticle->y;

		x = (x-Grid->xmin)/(Grid->xmax-Grid->xmin); // x is now between 0 and 1
		x = x*WaveNumber*2.0*PI;



		if (thisParticle->y<thickCrust) {
			//thisParticle->phase = 3;
		}
		else {
			thisParticle->phase = 0;
		}

		xP_from_xCorner = (thisParticle->x-xCorner);
		if (xP_from_xCorner>0 && (thisParticle->y < xP_from_xCorner*tanAngleCorner) ) {
			thisParticle->phase = 1;

		}
		if (y<thickDecollement) {
			thisParticle->phase = 2;
		}




		if (y<(Amplitude*sin(x+phase) + Grid->ymin+Thickness + thickDecollement) && xP_from_xCorner<0) {
			thisParticle->phase = 2;
		}


		if (y<(Amplitude*sin(x+phase) + Grid->ymin+Thickness) && xP_from_xCorner<0) {
			thisParticle->phase = 3;
		}


		END_PARTICLES


	}


	else {
		printf("Unknwon Setup %i in Particles_initPhase\n", Setup);
		exit(0);



	}
}

void Particles_initPassive(Grid* Grid, Particles* Particles)
{
	// Init a passive grid
	coord DX, DY;

	DY = (Grid->ymax-Grid->ymin)/20.0;
	DX = DY;//(Grid->xmax-Grid->xmin)/32.0;
	int passive;
	int dum;
	INIT_PARTICLE
	FOR_PARTICLES
	dum = (int)((thisParticle->x-Grid->xmin)/DX);

	passive = dum%2;
	//printf("x = %.2f, dum = %i, passive = %i\n", thisParticle->x-Grid->xmin, dum, passive);
	dum = (int)((thisParticle->y-Grid->ymin)/DY);
	passive += (dum)%2;
		if (passive==1) {
			if (thisParticle->phase != 0) { // quick fix for sticky air visualization
				thisParticle->passive = 0;
			} else {
				thisParticle->passive = 1;
			}

		} else {
			thisParticle->passive = 1;
		}
	END_PARTICLES
}









void Particles_initPhysics(Grid* Grid, Particles* Particles, BC* BCThermal)
{
	compute locY;
	compute H = (Grid->ymax-Grid->ymin);
	INIT_PARTICLE
	FOR_PARTICLES
		locY = (thisParticle->y-Grid->ymin)/H;
		thisParticle->T = 0*(  (1-locY)*BCThermal->TB + (locY)*BCThermal->TT  );


	END_PARTICLES
}





void Particles_teleportInsideTheDomain(Grid* Grid, Particles* Particles, Physics* Physics)
{
	// Due to advection error particles might end up outside the model boundaries. This function teleports them back inside
	bool change = false;

	SingleParticle* loopingParticle = NULL;
	int ParticleCounter = 0;
	compute locX, locY;
	compute dx = Grid->dx;
	compute dy = Grid->dy;
	int Imin, i, ix, iy;
	compute Min;

	INIT_PARTICLE
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
			thisParticle->T = (Physics->T[(ix)+(iy+1)*Grid->nxEC] + Physics->T[ix+1+(iy+1)*Grid->nxEC] + Physics->T[(ix)+(iy)*Grid->nxEC] + Physics->T[ix+1    +(iy)*Grid->nxEC])/4;
			thisParticle->sigma_xx_0 = (Physics->sigma_xx_0[(ix)+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[(ix)+(iy)*Grid->nxEC] + Physics->sigma_xx_0[ix+1    +(iy)*Grid->nxEC])/4;




			thisParticle->sigma_xy_0 = loopingParticle->sigma_xy_0; // not ideal


			thisParticle->nodeId = ix+iy*Grid->nxS;










			thisParticle->sigma_xx_0 = 0;
			thisParticle->sigma_xy_0 = 0;
		}

		END_PARTICLES



}

void Particles_deleteIfOutsideTheDomain(Grid* Grid, Particles* Particles)
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


void Particles_updateLinkedList(Grid* Grid, Particles* Particles, Physics* Physics)
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

	int oldNodeId;
	coord x, y;
	int ix, iy;
	SingleParticle* previousParticle;
	// Update the link list
	// =========================


	int ParticleCounter = 0;
	SingleParticle* thisParticle = NULL;
	int iNode = 0;
	int phase, passive;
	compute T;
	int TotNumParticles = 0;

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
	int dum;
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


	// Extra sweep to inject or delete particle
	// note: Not optimal this could be done during another sweep, for example during interpolation
	coord locX, locY;
	compute dx = Grid->dx;
	compute dy = Grid->dy;
	compute min = 1;
	int Imin;
	int i;
	TotNumParticles = 0;

	SingleParticle modelParticle;
	modelParticle.x = 0;
	modelParticle.y = 0;
	modelParticle.nodeId = 0;
	modelParticle.T = 0;
	modelParticle.sigma_xx_0 = 0;
	modelParticle.sigma_xy_0 = 0;
	modelParticle.phase = 0;
	modelParticle.passive = 0;
	modelParticle.next = NULL;

//#pragma omp parallel for private(iy, ix, iNode, thisParticle, ParticleCounter, locX, locY, min, Imin, i,  x, y, phase,passive, T) schedule(static,32)
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			iNode = ix  + (iy  )*Grid->nxS;
			thisParticle = Particles->linkHead[iNode];
			ParticleCounter=0;
			min = 1.0;
			Imin = 0;
			while (thisParticle != NULL) {
				ParticleCounter++;
				//TotNumParticles++;

				thisParticle = thisParticle->next;
			}

			if (ParticleCounter==0) {
				printf("Warning: node #%i is empty\n", iNode);
			}
			else if (ParticleCounter<Particles->nPC/1.5) {

				// find the closest particle to the node
				thisParticle = Particles->linkHead[iNode];

				ParticleCounter = 0;
				while (thisParticle != NULL) {


					locX = ((thisParticle->x-Grid->xmin)/dx - ix);
					locY = ((thisParticle->y-Grid->ymin)/dy - iy);

					if ( (locX*locX + locY*locY) < min) {
						min = (locX*locX + locY*locY);
						Imin = ParticleCounter;
					}

					ParticleCounter++;
					thisParticle = thisParticle->next;
				}

				// sweep again up to the closest particle
				thisParticle = Particles->linkHead[iNode];
				for (i=0; i<Imin; i++) {
					thisParticle = thisParticle->next;
				}



				//ix = iNode%Grid->nxS;
				//iy = (iNode-ix)/Grid->nxS;
				modelParticle.x = Grid->xmin + ix*Grid->dx  ;
				modelParticle.y = Grid->ymin + iy*Grid->dy  ;
				modelParticle.phase = thisParticle->phase; // the phase given to the particles is the phase of the head particle. Easy and fast but not optimal
				modelParticle.passive = thisParticle->passive; // the phase given to the particles is the phase of the head particle. Easy and fast but not optimal

				// This could be ok, but right now it's probably done with temperature not advected or something, which gives bad results;
				modelParticle.T = (Physics->T[(ix)+(iy+1)*Grid->nxEC] + Physics->T[ix+1+(iy+1)*Grid->nxEC] + Physics->T[(ix)+(iy)*Grid->nxEC] + Physics->T[ix+1    +(iy)*Grid->nxEC])/4;
				modelParticle.sigma_xx_0 = (Physics->sigma_xx_0[(ix)+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[ix+1+(iy+1)*Grid->nxEC] + Physics->sigma_xx_0[(ix)+(iy)*Grid->nxEC] + Physics->sigma_xx_0[ix+1    +(iy)*Grid->nxEC])/4;




				modelParticle.sigma_xy_0 = thisParticle->sigma_xy_0; // not ideal


				modelParticle.nodeId = iNode;

				//T = thisParticle->T;

				addSingleParticle(&Particles->linkHead[iNode], &modelParticle);
				Particles->n+=1;

			}

			else if (ParticleCounter>Particles->nPC*2.5) {
				thisParticle = Particles->linkHead[iNode];
				Particles->linkHead[iNode] = Particles->linkHead[iNode]->next;
				free(thisParticle);
				Particles->n-=1;
			}

		}




	}
	//printf("TotNumParticles = %i\n", TotNumParticles);











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


void Particles_advect(Particles* Particles, Grid* Grid, Physics* Physics)
{
	// Declarations
	// =========================
	int iNode;
	compute locX, locY, locX0, locY0;
	int Ix, Iy;
	int ix, iy;
	int iCell, i;
	int ixN, iyN;
	compute alphaArray[4];
	compute alpha;
	compute sigma_xx_corr, sigma_xy_corr;

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

				locX = fabs(locX)*2-1;
				locY = fabs(locY)*2-1;

				for (i=0;i<4;i++) {
					ixN = ix+IxN[i]*signX;
					iyN = iy+IyN[i]*signY;

					if (ixN+1>Grid->nxS || ixN<0 || iyN+1>Grid->nyS || iyN<0) {
						printf("error in Particles_advect: trying to access a non existing node\n");
						printf("IX = %i, IY = %i, locX = %.3f, locY = %.3f, iy = %i, IyN[i] = %i, signY = %i, ix = %i, IxN[i] = %i, signX = %i\n", ix+IxN[i]*signX, iy+IyN[i]*signY, locX, locY, iy, IyN[i], signY, ix, IxN[i], signX);
						printf("thisParticle->x = %.3f , y = %.3f \n", thisParticle->x, thisParticle->y);
						exit(0);
					}

					alphaArray[i]  = 0.5*Physics->dt*((Physics->Vy[ixN+1+iyN*Grid->nxVy]   - Physics->Vy[ixN+(iyN)*Grid->nxVy])/Grid->dx
							               	    - (Physics->Vx[ixN+(iyN+1)*Grid->nxVx] - Physics->Vx[ixN+(iyN)*Grid->nxVx])/Grid->dy);
					//printf("ix = %i, ixC = %i, iy = %i, iyC = %i, alphaArray[i] = %.3e\n", ix, ixC, iy, iyC, alphaArray[i]);
				}

				alpha =   ( .25*(1.0-locX)*(1.0-locY)*alphaArray[0]
						  + .25*(1.0-locX)*(1.0+locY)*alphaArray[1]
						  + .25*(1.0+locX)*(1.0+locY)*alphaArray[2]
						  + .25*(1.0+locX)*(1.0-locY)*alphaArray[3] );

				sigma_xx_corr = - thisParticle->sigma_xy_0 * 2 * alpha;
				sigma_xy_corr = + thisParticle->sigma_xx_0 * 2 * alpha;

				//printf("alpha = %.3e, alphaArray[0] = %.3e, alphaArray[1] = %.3e, alphaArray[2] = %.3e, alphaArray[3] = %.3e\n", alpha, alphaArray[0], alphaArray[1], alphaArray[2], alphaArray[3]);

				thisParticle->sigma_xx_0 += sigma_xx_corr;
				thisParticle->sigma_xy_0 += sigma_xy_corr;


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
								   + .25*(1.0+locX)*(1.0-locY)*Physics->Vx[Ix+1+(Iy  )*Grid->nxVx] ) * Physics->dt;


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
								   + .25*(1.0+locX)*(1.0-locY)*Physics->Vy[Ix+1+(Iy  )*Grid->nxVy] ) * Physics->dt;


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
	INIT_PARTICLE
	FOR_PARTICLES
	if ((thisParticle->x-Grid->xmin)<0 ) {
		//printf("#### Loop the loop ####\n");
		thisParticle->x += Grid->xmax-Grid->xmin;
	}
	else if ((thisParticle->x-Grid->xmax)>0 ) {
		thisParticle->x -= Grid->xmax-Grid->xmin;
	}

	if (thisParticle->y<Grid->ymin) {
		thisParticle->y = Grid->ymin+0.1*Grid->dy;
	} else if (thisParticle->y>Grid->ymax) {
		thisParticle->y = Grid->ymax-0.1*Grid->dy;
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







