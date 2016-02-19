/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"

void allocateVisuMemory( Visu* Visu )
{
	Visu->elements      = (GLuint*)  malloc(Visu->ntrivert    * sizeof( GLuint ));
}




void freeVisuMemory( Visu* Visu )
{
	free(Visu->elements);
}



void initVisualization(Visu* Visu, Grid* Grid)
{

	// Create the element array
	// Fill elements, loop through cells
	int ix,iy,C;
	int nxS = Grid->nxS;
	int nyS = Grid->nyS;


	C = 0;
	for (iy=0;iy<nyS-1;iy++){
		for (ix=0;ix<nxS-1;ix++){
			// Triangle 1
			Visu->elements[C+0] = ix+iy*nxS;
			Visu->elements[C+1] = ix+1+iy*nxS;
			Visu->elements[C+2] = (iy+1)*nxS+ix;
			// Triangle 2
			Visu->elements[C+3] = ix+1+iy*nxS;
			Visu->elements[C+4] = (iy+1)*nxS+ix+1;
			Visu->elements[C+5] = (iy+1)*nxS+ix;
			C = C+6;
		}
	}


}



