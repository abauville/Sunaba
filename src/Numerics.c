/*
 * Numerics.c
 *
 *  Created on: Jun 02, 2016
 *      Author: abauville
 */


#include "stokes.h"

void Numerics_init(Numerics* Numerics)
{
	int iLS;

	Numerics->glob = (compute*) malloc((Numerics->nLineSearch+1) * sizeof(compute));

	for (iLS= 0; iLS < Numerics->nLineSearch; ++iLS) {
		//compute a, the globalization parameter;
		if (iLS!=Numerics->nLineSearch)
			Numerics->glob[iLS] = Numerics->maxCorrection - Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS);
	}

}


void Numerics_freeMemory(Numerics* Numerics)
{
	free(Numerics->glob);
}
