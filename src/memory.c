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






