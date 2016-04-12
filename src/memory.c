#include "stokes.h"

void Memory_allocateMain(Grid* Grid, Particles* Particles, Physics* Physics, EqSystem* EqStokes, Numbering* NumStokes, Numbering* NumThermal)
{


	Particles->linkHead 	= (SingleParticle**) malloc( Grid->nSTot 		* sizeof(  SingleParticle*  ) ); // array of pointers to particles

	int i;
	//SingleParticle* A=NULL;
	for (i=0;i<Grid->nSTot;i++) {
		Particles->linkHead[i] = NULL;
	}


	Physics->eta 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->eta0 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->n 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->rho 			= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->k 				= (compute*) 	malloc( Grid->nECTot * sizeof(compute) );
	Physics->etaShear		= (compute*) 	malloc( Grid->nxS*Grid->nyS * sizeof(compute) );

	NumStokes->map  		= (int*) 		malloc(EqStokes->nEqIni 	* sizeof(int)); // Numbering map
	NumThermal->map  		= (int*) 		malloc((Grid->nxC+2)*(Grid->nyC+2) 	* sizeof(int)); // Numbering map

	Physics->Vx 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
	Physics->P 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->T 				= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->DT 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	Physics->G 			= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->GShear		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->sigma_xx_0  	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->sigma_xy_0		= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );
	Physics->Dsigma_xx_0 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->Dsigma_xy_0 	= (compute*) 	malloc( Grid->nSTot 		* sizeof(compute) );

	Physics->cohesion 		= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );
	Physics->frictionAngle 	= (compute*) 	malloc( Grid->nECTot 		* sizeof(compute) );

	// Initialize stuff
	//int i;
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0;
	}
	for (i = 0; i < Grid->nECTot; ++i) {
		Physics->P[i]  = 0;
		Physics->T[i]  = 0;
		Physics->DT[i] = 0;
		Physics->sigma_xx_0[i] = 0;
	}
	for (i = 0; i < Grid->nSTot; ++i) {
		Physics->sigma_xy_0[i] = 0;
	}
}


void Memory_freeMain(Particles* Particles, Physics* Physics, Numbering* NumStokes, Numbering* NumThermal, BC* BC, Grid* Grid)
{



	free( Physics->eta );
	free( Physics->eta0 );
	free( Physics->n );
	free( Physics->rho );
	free( Physics->k );
	free( Physics->etaShear );

	free( NumStokes->map );
	free( NumThermal->map );


	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P );
	free(Physics->T );
	free(Physics->DT );

	free(Physics->G );
	free(Physics->GShear );


	free(Physics->sigma_xx_0 );
	free(Physics->sigma_xy_0 );
	free(Physics->Dsigma_xx_0 );
	free(Physics->Dsigma_xy_0 );

	free(Physics->cohesion);
	free(Physics->frictionAngle);


	free(BC->list);
	free(BC->value);
	free(BC->type);

	printf("Free Particles..\n");
	Particles_freeAllSingleParticles(Particles, Grid);
	free( Particles->linkHead );
}



// Linked List function
void addToLinkedList(LinkedNode** pointerToHead, int x)
{
	// Adds a node at the beginning of a linked list
	LinkedNode* temp = (LinkedNode*) malloc(sizeof(LinkedNode));
	temp->data = x;
	temp->next = NULL;
	if (*pointerToHead != NULL) {
		temp->next = *pointerToHead;
	}
	*pointerToHead = temp;

}

void freeLinkedList(LinkedNode* head)
{
   LinkedNode* temp;

   while (head != NULL)
    {
       temp = head;
       head = head->next;
       free(temp);
    }

}





