/*
 * AdvectionDiffusionBackwardEulerFactory.cpp
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

AdvectionDiffusionBackwardEulerFactory::AdvectionDiffusionBackwardEulerFactory(double aX, double aY, double D, double deltaT){
	this->aX = aX;
	this->aY = aY;
	this->D = D;
	this->deltaT = deltaT;
}

void AdvectionDiffusionBackwardEulerFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}

	CellToCellOperator * L = laplacian.get(g);
	CellToCellOperator I, cL;
	cL = deltaT*D*(*L);
	I = CellToCellOperator::getIdentity(g);
	CellToCellOperator l = I - cL;
	CellToCellOperator r = I;

	lhs[g] = new CellToCellOperator(l);
	rhs[g] = new CellToCellOperator(r);
}

}
}
