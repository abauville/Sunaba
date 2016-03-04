#include "stokes.h"

void Memory_allocateMain(Grid* Grid, Particles* Particles, Physics* Physics, EqSystem* EqSystem, Numbering* Numbering)
{

	Particles->xy 			= (coord*) 		malloc( Particles->n * 2 	* sizeof( coord ) );
	Particles->phase 		= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );

	Particles->cellId 		= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );
	Particles->linkNext 	= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );
	Particles->linkHead 	= (int*) 		malloc( Grid->nCTot 		* sizeof(  int  ) );

	Physics->eta 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );
	Physics->rho 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );
	Physics->etaShear		= (compute*) 	malloc( Grid->nxS*Grid->nyS * sizeof(compute) );

	Numbering->map  		= (int*) 		malloc(EqSystem->nEqIni 	* sizeof(int)); // Numbering map

	Physics->Vx 			= (compute*) 	malloc( Grid->nVxTot 		* sizeof(compute) );
	Physics->Vy 			= (compute*) 	malloc( Grid->nVyTot 		* sizeof(compute) );
	Physics->P 				= (compute*) 	malloc( Grid->nCTot 		* sizeof(compute) );


}


void Memory_freeMain(Particles* Particles, Physics* Physics, Numbering* Numbering)
{

	free( Particles->phase );
	free( Particles->xy    );

	free( Particles->cellId );
	free( Particles->linkNext );
	free( Particles->linkHead );

	free( Physics->eta );
	free( Physics->rho );
	free( Physics->etaShear );

	free( Numbering->map );


	free(Physics->Vx);
	free(Physics->Vy);
	free(Physics->P);





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





