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
			//Numerics->glob[iLS] = Numerics->maxCorrection - Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS);
			Numerics->glob[iLS] = Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS+1);
	}

}


void Numerics_freeMemory(Numerics* Numerics)
{
	free(Numerics->glob);
}


int Numerics_updateBestGlob(Numerics* Numerics, EqSystem* EqStokes, int iLS)
{
	// returns a break condition
	// if 1, the linesearch loop should break
	int Break = 0;

	if (EqStokes->normResidual<Numerics->minRes) {
		Numerics->glob[Numerics->nLineSearch] = Numerics->glob[iLS];
		Numerics->minRes = EqStokes->normResidual;
		if (iLS==Numerics->nLineSearch-1) {// if the last one is the best one then don't recompute, i->e-> next one would be the same
			Break = 1;
		}
	}
	printf("a = %.2f, |F|/|b|: %.2e\n", Numerics->glob[iLS], EqStokes->normResidual);


	if (Numerics->itNonLin==0) {
		Numerics->normRes0 = EqStokes->normResidual;
		if (Numerics->timeStep==0)
			Numerics->normResRef = EqStokes->normResidual;

		if (Numerics->maxNonLinearIter==1) {
			Break = 1;
		}
	}

	return Break;

}
