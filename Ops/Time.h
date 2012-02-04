/*
 * OperatorFactory.h
 *
 *  Created on: Jan 25, 2011
 *      Author: fogelson
 */

#ifndef TIME_H_
#define TIME_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <map>
#include <blitz/tinyvec-et.h>
#include "OperatorFactory.h"

using namespace blitz;
using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

class Time;

class Time{
	double deltaT, t0;
	int current;
public:
	double getDeltaT();
	double getTime();
	void step();
	void operator++ ();
};

template<class T>
class TimeDependentOperatorFactory : public OperatorFactory<T>{
	int stepsToKeep;
	OperatorFactory<T> factories[];
	Time * time;
};


}
}

#endif /* TIME_H_ */
