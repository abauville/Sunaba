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


}


void Numerics_freeMemory(Numerics* Numerics)
{
	//free(Numerics->glob);
}



void Numerics_LineSearch_chooseGlob(Numerics* Numerics, EqSystem* EqStokes) {

	int currentState = Numerics->lsState;
	int nextState;
	compute a = Numerics->lsGlob; // globalization factor
	compute lowerbound = Numerics->lsLowerBound;
	compute upperbound = Numerics->lsUpperBound;

	compute Res = EqStokes->normResidual;
	compute bestRes = Numerics->lsBestRes;

	compute lastRes = Numerics->lsLastRes;
	compute lastGlob = Numerics->lsLastGlob;

	compute tolImprovement = Numerics->lsResTolImprovement; // minimum tolerated improvement

	int counterUp = Numerics->lsCounterUp;
	compute tolDiverge = Numerics->lsTolDiverge;

	compute bestRes_It = Numerics->lsbestRes_It;

	printf("a = %.2f, |F|/|b|: %.2e\n", Numerics->lsGlob, EqStokes->normResidual);




	switch (currentState) {
	case -2:
	case -1:
		//bestRes = Res;
		nextState = 0;
		break;
	case  0:
		if (Res<lastRes) {
			// Stop and reinit
			nextState = -1;

			if (Res>lastRes) {
				nextState = 1;
				upperbound = Numerics->lsGlobStart;
				//lowerbound = 0.0;// Numerics->lsGlobStart;
				//bestRes = 1e15;
			}

		} else {
			upperbound = Numerics->lsGlobStart;
			nextState = 1;
		}
		break;
	case  1:
		/*
		if (Res> bestRes*(1.0+tolImprovement)) {
			if (a>lastGlob) {
				// move back down
				upperbound = a;
			} else {
				// move back up
				lowerbound = a;
			}
			nextState = 1;
			*/
		//if (Res > bestRes_It*(1.0+tolImprovement)) {
		if (fabs(Res-lastRes)/lastRes>tolImprovement) {
			printf("bestRest_It = %.2e, bestRest_It+tol = %.2e\n", bestRes_It, bestRes_It*(1.0+tolImprovement));
			nextState = 11;
		} else {
			//bestRes =
			nextState = -1;
		}

		break;
	case 11:
		if (Res<lastRes) {
			// move up
			lowerbound = a;
		} else {
			// move down
			upperbound = a;
		}
		nextState = 1;
		/*
		if (counterUp == 1) {
			nextState = -1;
		} else {
			if (Res<bestRes) {
				if (fabs(Res-bestRes)/bestRes > tolImprovement) {
					// Go down
					upperbound = a;
					nextState = 1;
					bestRes = Res;
				} else {
					if (Res>lastRes*(1.0+tolImprovement/5.0)) {
						// Go down
						upperbound = a;
						nextState = 1;
						//bestRes = Res;
					} else {
						// Stop and reinit
						nextState = -1;
						upperbound = a;
					}
				}
			} else {
				// Go up then stop
				lowerbound = a;
				counterUp = 1;
				nextState = 1;
			}
		}
		*/



		break;
	}


	if (Res<Numerics->lsbestRes_It) {
		Numerics->lsbestRes_It = Res;
	}

	if (Numerics->maxNonLinearIter==1) {
		nextState = -1;
	}



	//if (nextState == -1 && (lastRes-Res)/lastRes< tolImprovement/5.0) {
	if (nextState == -1 && (Res - bestRes)/bestRes > tolDiverge) {
		// If despite all efforts the solution could not be improved (maybe stuck in a local minimum), then go to the next time step
		nextState = -2;
	}


	if (nextState == 1 && a<Numerics->lsGlobMin) {
		nextState = -1;
	}


	if (Numerics->itNonLin == 0) {
		nextState = -1;
	}



	Numerics->lsLastGlob = a;

	if (Res<bestRes) {
		bestRes = Res;
	}


	Numerics->lsCounter++;


	printf("nextState = %i\n", nextState);

	switch (nextState) {
	case -2: // go to the next time step
		Numerics->lsBestGlob = a;
		Numerics->lsCounter = 0;
		a = 1.00;
		bestRes = 1E15;
		lowerbound = 0.0;
		upperbound = 1.0;
		counterUp = 0;
		Numerics->lsLastRes = Res;
		Numerics->lsbestRes_It = 1E15;
		break;
	case -1:
		Numerics->lsBestGlob = a;
		a = Numerics->lsGlobStart;
		Numerics->lsCounter = 0;
		//bestRes = 1E15;
		lowerbound = 0.0;
		upperbound = 1.0;
		counterUp = 0;
		Numerics->lsLastRes = Res;
		Numerics->lsbestRes_It = 1E15;
		break;
	case 0: 					// run 1.0 case
		a = 1.00;
		Numerics->lsLastRes = Res;
		break;
	case 1: 					// search for the best a
		a = (lowerbound+upperbound)/2;
		Numerics->lsLastRes = Res;
		break;
	case 11:
		a = a+0.01;
		Numerics->lsLastRes = Res;
		break;
	}





	Numerics->lsState = nextState;
	Numerics->lsGlob = a; // globalization factor
	Numerics->lsLowerBound = lowerbound;
	Numerics->lsUpperBound = upperbound;

	Numerics->lsBestRes = bestRes;



	Numerics->lsCounterUp = counterUp;







}


/*

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


		//if (Numerics->maxNonLinearIter==1) {
		//	Break = 1;
		//}

	}

	return Break;

}
*/
