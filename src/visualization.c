/*
 * visualization.c
 *
 *  Created on: Feb 19, 2016
 *      Author: abauville
 */

#include "stokes.h"

void Visu_allocateMemory( Visu* Visu, Grid* Grid )
{
	Visu->elements      = (GLuint*)   malloc(Visu->ntrivert    		* sizeof( GLuint ));
	Visu->U             = (GLfloat*)  malloc(Grid->nxS*Grid->nyS    * sizeof( GLfloat ));
	Visu->vertices      = (GLfloat*)  malloc(Grid->nxS*Grid->nyS*2  * sizeof( GLfloat ));
}




void Visu_freeMemory( Visu* Visu )
{
	free(Visu->elements);
	free(Visu->U);
}



void Visu_init(Visu* Visu, Grid* Grid)
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

	C =0;
	for (iy = 0; iy < nyS; ++iy) {
		for (ix = 0; ix < nxS; ++ix) {
			Visu->vertices[C  ] = (Grid->xmin + ix*Grid->dx);
			Visu->vertices[C+1] = (Grid->ymin + iy*Grid->dy);
			C += 2;
		}
	}

}


void Visu_plotCenterValue(Visu* Visu, Grid* Grid, compute* Value)
{
	// UC is a scalar value defined on the center grid
	// Declarations
	// =========================
	int ix, iy;
	int I;

	int iNW, iNE, iSW, iSE;
	// Value interpolated on the center nodes
	// ======================================
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		for (ix = 1; ix < Grid->nxS-1; ++ix) {
			I = ix + iy*Grid->nxS;
			iNW = (ix-1)+ iy   *Grid->nxC;
			iNE = ix    + iy   *Grid->nxC;
			iSW = (ix-1)+(iy-1)*Grid->nxC;
			iSE = ix    +(iy-1)*Grid->nxC;
			Visu->U[I] = (Value[iNW] + Value[iNE] + Value[iSW] + Value[iSE])/4;
		}
	}
	// Value extrapolated on the lower boundary
	// ======================================
	// o: centered value
	// x: value extrapolated (in 1D) fron the o nodes
	// X: valu interpolated between the two x
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1b   and 2b
	// |       |       |
	// | - - - | - - - | -
	// |       |       |
	// |   o   |   o   |       nodes 1a   and 1b
	// |       |       |
	// | - - - | - - - | -
	// |       |       |
	// | - x - X - x - |       nodes tempa and tempb
	//
	iy = 0;
	compute temp1, temp2;
	int i1a, i1b, i2a, i2b;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy+1)*Grid->nxC;
		i1a = (ix-1)+ iy   *Grid->nxC;
		i2b =  ix   +(iy+1)*Grid->nxC;
		i2a =  ix   + iy   *Grid->nxC;


		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}
	// Value extrapolated on the upper boundary
	// ======================================
	//   x  X  x
	//  1a    2a
	//  1b    2b
	iy = Grid->nyS-1;
	for (ix = 1; ix < Grid->nxS-1; ++ix) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-1)+(iy-2)*Grid->nxC;
		i1a = (ix-1)+(iy-1)*Grid->nxC;
		i2b =  ix   +(iy-2)*Grid->nxC;
		i2a =  ix   +(iy-1) *Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Value extrapolated on the left boundary
	// ======================================
	//  x 1a   1b
	//  X
	//  x 2a   2b
	ix = 0;
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		I = ix + iy*Grid->nxS;
		i1b =  ix   +(iy  )*Grid->nxC;
		i1a = (ix+1)+(iy  )*Grid->nxC;
		i2b =  ix   +(iy-1)*Grid->nxC;
		i2a = (ix+1)+(iy-1)*Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Value extrapolated on the left boundary
	// ======================================
	//  1b   1a x
	//          X
	//  2b   2a x
	ix = Grid->nxS-1;
	for (iy = 1; iy < Grid->nyS-1; ++iy) {
		I = ix + iy*Grid->nxS;
		i1b = (ix-2)+(iy  )*Grid->nxC;
		i1a = (ix-1)+(iy  )*Grid->nxC;
		i2b = (ix-2)+(iy-1)*Grid->nxC;
		i2a = (ix-1)+(iy-1)*Grid->nxC;
		temp1 = Value[i1a] - (Value[i1b] - Value[i1a])/2;
		temp2 = Value[i2a] - (Value[i2b] - Value[i2a])/2;
		Visu->U[I] = (temp1+temp2)/2;
	}

	// Lower left corner
	coord diag = sqrt(Grid->dx*Grid->dx+Grid->dy*Grid->dy);
	//          1b
	//      1a
	//   X
	ix = 0; iy = 0;
	I = ix + iy*Grid->nxS;
	i1b = (ix+1)+(iy+1)*Grid->nxC;
	i1a =  ix   +(iy  )*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Lower right corner
	//  1b
	//      1a
	//          X
	ix = Grid->nxS-1; iy = 0;
	I = ix + iy*Grid->nxS;
	i1b = (ix-2)+(iy+1)*Grid->nxC;
	i1a = (ix-1)+(iy  )*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Upper left corner
	//  X
	//      1a
	//          1b
	ix = 0; iy = Grid->nyS-1;
	I = ix + iy*Grid->nxS;
	i1b = (ix+1)+(iy-2)*Grid->nxC;
	i1a =  ix   +(iy-1)*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

	// Upper right corner
	//          X
	//      1a
	//  1b
	ix = Grid->nxS-1; iy = Grid->nyS-1;
	I = ix + iy*Grid->nxS;
	i1b = (ix-2)+(iy-2)*Grid->nxC;
	i1a = (ix-1)+(iy-1)*Grid->nxC;
	Visu->U[I] = Value[i1a] - (Value[i1b] - Value[i1a])/2;

}


