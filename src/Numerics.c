/*
 * Numerics.c
 *
 *  Created on: Jun 02, 2016
 *      Author: abauville
 */


#include "stokes.h"



void Numerics_init(Numerics* Numerics)
{

	//Numerics->glob = (compute*) malloc((Numerics->nLineSearch+1) * sizeof(compute));


	/*
	for (iLS= 0; iLS < Numerics->nLineSearch; ++iLS) {
		//compute a, the globalization parameter;
		if (iLS!=Numerics->nLineSearch)
			Numerics->glob[iLS] = Numerics->maxCorrection - Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS);
			//Numerics->glob[iLS] = Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS+1);
	}
	*/
	Numerics->lsGlobStart 			= 0.99;

	if (Numerics->maxNonLinearIter == 1) {
		Numerics->lsGlobStart 			= 1.0;
	}

	Numerics->lsResTolImprovement 	= 0.05;
	Numerics->lsGlobMin 			= 0.02;
	Numerics->lsTolDiverge 			= 0.8; // if residual goes above (1 + lsTolDiverge)*bestRes, breaks; bestRes is reset only every time step



	Numerics->lsState = -1;
	Numerics->lsCounterUp = 0;
	Numerics->lsCounter   = 0;


	Numerics->lsBestGlob = 0.0;
	Numerics->lsBestRes = 0.0;

	Numerics->lsLowerBound = -.2;
	Numerics->lsUpperBound = 1.0;



	Numerics->lsGlob = 1.00; // should be 1.0, because of the inconsistency between BC when passing a time step (due to box update)
	Numerics->lsBestRes = 1E15;
	Numerics->lsState = 0;
	Numerics->lsLowerBound = 0.0;
	Numerics->lsUpperBound = 1.0;
	Numerics->lsState = -1;
	Numerics->lsCounter = 0;


	Numerics->lsLastRes = 1E15;
	Numerics->lsLastGlob = 0.0;
	Numerics->lsbestRes_It = 1e15;
	Numerics->stallingCounter  = 0;
	Numerics->stallingCounter = false;

}


void Numerics_Memory_free(Numerics* Numerics)
{
	//free(Numerics->glob);
}

