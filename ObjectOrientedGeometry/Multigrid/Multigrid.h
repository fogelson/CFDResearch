/*
 * Smoothers
 *
 *  Created on: Oct 26, 2011
 *      Author: fogelson
 */

#ifndef MULTIGRID_H_
#define MULTIGRID_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <map>
#include <blitz/tinyvec-et.h>
#include "../Ops/Operators.h"
#include "../Ops/OperatorFactory.h"
#include "../TransferOperators/TransferOperators.h"
#include "../Smoothers/Smoothers.h"

using namespace blitz;
using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{

class MultigridSolver;

class MultigridSolver{
	StenciledSmoother * smoother;
	Interpolator * interpolator;
	Restrictor * restrictor;

public:
	MultigridSolver(StenciledSmoother * smoother, Interpolator * interpolator, Restrictor * restrictor);
	void vCycle(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray & rhs, int v1, int v2, Grid * fine, SplitOperatorFactory<CellToCellOperator> * fac);
};

}
}

#endif /* MULTIGRID_H_ */
