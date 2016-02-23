#include "stokes.h"

void allocateMemory(Grid* Grid, Particles* Particles, Physics* Physics)
{

	Particles->xy 			= (coord*) 		malloc( Particles->n * 2 	* sizeof( coord ) );
	Particles->phase 		= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );

	Particles->oldCellId 	= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );
	Particles->newCellId 	= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );
	Particles->linkNext 	= (int*) 		malloc( Particles->n 		* sizeof(  int  ) );
	Particles->linkHead 	= (int*) 		malloc( Grid->nCTot 		* sizeof(  int  ) );

	Physics->eta 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );
	Physics->rho 			= (compute*) 	malloc( Grid->nxC*Grid->nyC * sizeof(compute) );



}


void freeMemory(Particles* Particles, Physics* Physics)
{

	free( Particles->phase );
	free( Particles->xy    );

	free( Particles->oldCellId );
	free( Particles->newCellId );
	free( Particles->linkNext );
	free( Particles->linkHead );

	free( Physics->eta );
	free( Physics->rho );

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





