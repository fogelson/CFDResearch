/*
 * FENETimingVariables.h
 *
 *  Created on: Jul 5, 2011
 *      Author: fogelson
 */

#ifndef FENETIMING_H_
#define FENETIMING_H_

#define FeneTiming

#ifdef FeneTiming

#include <time.h>
#include <ctime>
using namespace std;

namespace CFD{
	namespace Timing{
		int stressAtPoint = 0;
		double stressAtPointTime = 0;

		int callSolver = 0;
		double callSolverTime = 0;

		int callUpdatePolymers = 0;
		double callUpdatePolymersTime = 0;
	}
}

#endif

#endif /* FENETIMING_H_ */
