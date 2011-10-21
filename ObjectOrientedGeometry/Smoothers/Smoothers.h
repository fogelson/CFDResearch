/*
 * Smoothers
 *
 *  Created on: Oct 20, 2011
 *      Author: fogelson
 */

#ifndef SMOOTHERS_H_
#define SMOOTHERS_H_

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

class StenciledSmoother;

class GSLex;

class StenciledSmoother{
public:
	virtual ~StenciledSmoother(){}
	virtual void smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, CellToCellOperator & C, int its) = 0;
	//virtual CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, CellToCellOperator & C, int its) = 0;
};

class GSLex : public StenciledSmoother{
public:
	void smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, CellToCellOperator & C, int its);
};

}
}

#endif /* OPERATORS_H_ */
