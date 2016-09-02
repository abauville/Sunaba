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

	//Numerics->glob = (compute*) malloc((Numerics->nLineSearch+1) * sizeof(compute));


	/*
	for (iLS= 0; iLS < Numerics->nLineSearch; ++iLS) {
		//compute a, the globalization parameter;
		if (iLS!=Numerics->nLineSearch)
			Numerics->glob[iLS] = Numerics->maxCorrection - Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS);
			//Numerics->glob[iLS] = Numerics->maxCorrection/(Numerics->nLineSearch) * (iLS+1);
	}
	*/

	Numerics->lsState = -1;
	Numerics->lsCounterUp = 0;
	Numerics->lsCounter   = 0;

	Numerics->lsGlob = 0.95;// Initial value
	Numerics->lsBestGlob = 0.0;
	Numerics->lsBestRes = 0.0;

	Numerics->lsLowerBound = 0.0;
	Numerics->lsUpperBound = 1.0;
}


void Numerics_freeMemory(Numerics* Numerics, EqSystem* EqStokes)
{
	//free(Numerics->glob);
}



void Numerics_LineSearch_chooseGlob(Numerics* Numerics, EqSystem* EqStokes) {

	int state = Numerics->lsState;
	int nextState;
	compute a; // globalization factor
	compute lowerbound, upperbound;

	compute Res = EqStokes->normResidual;
	compute bestRes = Numerics->lsBestRes;

	int Break = 0;

	switch (state) {
		case -1: 					// begin
			a = 0.95
			bestRes = Res;
			nextState = 0;
			break;
		case 0: 					// run 1.0 case
			a = 1.00;

			if (Res<bestRes) {
				// Stop and reinit
				Break = 1;
				nextState = -1;
			} else {
				nextState = 1;
			}
			break;
		case 1: 					// search for the best a

			break;

		case 2:						// exit

			break;
	}

	Numerics->lsBestRes = bestRes;

	return Break;

	state 0
	lastRes = Res
	a = 1.00
	run -> Res
	if Res<LastRes stop
	else lowerbound = 0, upperbound = 1 , go to state = 1





	state 1
	lastRes = Res
	a = (lowerbound+upperbound)/2
	run -> Res
	if Res<LastRes upperbound = a
		CounterUP = 0
	else lowerbound = a
		CountUP += 1

	if CounterUP > 1
		go to state 3
	else stay at state 1

	if fabs(Res-LastRes)/LastRes<Tol stop






	state 3
	a = best a
	run
	re init state to -1
	stop


}




int Numerics_updateBestGlob(Numerics* Numerics, EqSystem* EqStokes, int* iLS)
{




	// returns a break condition
	// if 1, the linesearch loop should break
	int Break = 0;

	compute oldBest = Numerics->minRes;

	if (*iLS==Numerics->nLineSearch) {
		Break = 1;
	}

	if (EqStokes->normResidual<Numerics->minRes) {


		Numerics->glob[Numerics->nLineSearch] = Numerics->glob[*iLS];
		Numerics->minRes = EqStokes->normResidual;
		if (*iLS==Numerics->nLineSearch-1) {// if the last one is the best one then don't recompute, i->e-> next one would be the same
			Break = 1;
		}
		printf("a = %.2f, |F|/|b|: %.2e\n", Numerics->glob[*iLS], EqStokes->normResidual);
		//printf("abs((oldBest-EqStokes->normResidual)/oldBest) = %.1f\n", fabs((oldBest-EqStokes->normResidual)/oldBest*100.0));

		if (*iLS>0 && fabs((oldBest-EqStokes->normResidual)/oldBest)<0.05) {
			// If the solution is not better by at least a certain relative amount, then end it
			//Break = 1;
			*iLS = Numerics->nLineSearch-1; // Perform it one last time at the best glob then exit

		}



	} else {
		printf("a = %.2f, |F|/|b|: %.2e\n", Numerics->glob[*iLS], EqStokes->normResidual);
		*iLS = Numerics->nLineSearch-1; // Perform it one last time at the best glob then exit
	}





	if (Numerics->itNonLin==0) {
		Numerics->normRes0 = EqStokes->normResidual;
		if (Numerics->timeStep==0)
			Numerics->normResRef = EqStokes->normResidual;

		/*
		if (Numerics->maxNonLinearIter==1) {
			Break = 1;
		}
		*/
	}

	return Break;

}
