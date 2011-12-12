/*
 * DiffusionBackwardEulerFactory.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

DiffusionBackwardEulerFactory::DiffusionBackwardEulerFactory(double D, double deltaT){
	this->D = D;
	this->deltaT = deltaT;
}

void DiffusionBackwardEulerFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	Divergence div(g);
	Gradient grad(g);
	CellToCellOperator I, L, deltaTDL;
	L = div(grad);
	deltaTDL = (deltaT*D)*L;
	I = CellToCellOperator::getIdentity(g);
	CellToCellOperator l = I - deltaTDL;
	CellToCellOperator r = I;
	lhs[g] = new CellToCellOperator(l);
	rhs[g] = new CellToCellOperator(r);
}

}
}
