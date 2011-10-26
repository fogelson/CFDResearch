/*
 * CrankNicholsonDiffusionFactory
 *
 *  Created on: Oct 25, 2011
 *      Author: fogelson
 */


#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

CrankNicholsonDiffusionFactory::CrankNicholsonDiffusionFactory(double deltaT){
	this->deltaT = deltaT;
}

void CrankNicholsonDiffusionFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}

	Divergence div(g);
	Gradient grad(g);

	CellToCellOperator I, L, DL;
	CellToCellOperator lhs, rhs;
	I = CellToCellOperator::getIdentity(g);
	L = div(grad);
	DL = (deltaT/2)*L;
	lhs = I - DL;
	rhs = I + DL;

	lhsOperators[g] = new CellToCellOperator(lhs);
	rhsOperators[g] = new CellToCellOperator(rhs);
}

}
}
