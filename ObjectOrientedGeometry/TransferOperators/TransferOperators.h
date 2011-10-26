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
class PiecewiseConstantInterpolator;

class Interpolator{
public:
	virtual ~Interpolator(){}
	virtual void doInterpolate(CellDoubleArray & uF, CellDoubleArray & uC, Grid * coarse, Grid * fine) = 0;
};

class PiecewiseConstantInterpolator : public Interpolator{
public:
	void doInterpolate(CellDoubleArray & uC, CellDoubleArray & uF, Grid * coarse, Grid * fine);
};

class Restrictor{
public:
	virtual ~Restrictor(){}
	virtual void doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * coarse, Grid * fine) = 0;
};

class VolumeWeightedRestrictor{
public:
	void doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * coarse, Grid * fine);
};

}
}

#endif /* TRANSFEROPERATORS_H_ */
