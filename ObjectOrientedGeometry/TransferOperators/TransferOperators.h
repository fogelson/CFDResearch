/*
 * TransferOperators.h
 *
 *  Created on: Oct 21, 2011
 *      Author: fogelson
 */

#ifndef TRANSFEROPERATORS_H_
#define TRANSFEROPERATORS_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <map>
#include <blitz/tinyvec-et.h>
#include "../Ops/Operators.h"

using namespace blitz;
using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{

class Interpolator;
class Restrictor;

class VolumeWeightedRestrictor;

class Interpolator{
public:
	virtual ~Interpolator(){}
	virtual void doInterpolate(CellDoubleArray & uF, CellDoubleArray & uC, Grid * fine, Grid * coarse) = 0;
};

class Restrictor{
public:
	virtual ~Restrictor(){}
	virtual void doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * fine, Grid * coarse) = 0;
};

class VolumeWeightedRestrictor{
public:
	void doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * fine, Grid * coarse);
};

}
}

#endif /* TRANSFEROPERATORS_H_ */
