#include "stokes.h"

void Memory_allocateMain(Grid* Grid, Particles* Particles, Physics* Physics, EqSystem* EqSystem, Numbering* Numbering)
{


	Particles->linkHead 	= (SingleParticle**) malloc( Grid->nCTot 		* sizeof(  SingleParticle*  ) ); // array of pointers to particles
	/*
	int i;
	//SingleParticle* A=NULL;
	for (i=0;i<Grid->nCTot;i++) {
		Particles->linkHead[i] = NULL;
	}
	*/

	Physics->eta 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );
	Physics->rho 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );
	Physics->etaShear		= (compute*) 	malloc( Grid->nxS*Grid->nyS * sizeof(compute) );

	Numbering->map  		= (int*) 		malloc(EqSystem->nEqIni 	* sizeof(int)); // Numbering map

	Physics->Vx 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
	Physics->P 				= (compute*) 	malloc( Grid->nCTot 		* sizeof(compute) );

	// Initialize Vx, Vy, P
	int i;
	for (i = 0; i < Grid->nVxTot; ++i) {
		Physics->Vx[i] = 0;
	}
	for (i = 0; i < Grid->nVyTot; ++i) {
		Physics->Vy[i] = 0;
	}
	for (i = 0; i < Grid->nCTot; ++i) {
		Physics->P[i] = 0;
	}

}


void Memory_freeMain(Particles* Particles, Physics* Physics, Numbering* Numbering, BC* BC, Grid* Grid)
{



	free( Physics->eta );
	free( Physics->rho );
	free( Physics->etaShear );

	free( Numbering->map );


	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P);


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





