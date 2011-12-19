/*
 * AdvectionBackwardEulerFactory.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

void AdvectionBackwardEulerFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}




	//CellToCellOperator I, cL;
	//cL = deltaT*D*(*L);
	//I = CellToCellOperator::getIdentity(g);
	//CellToCellOperator l = I - cL;
	//CellToCellOperator r = I;

	//lhs[g] = new CellToCellOperator(l);
	//rhs[g] = new CellToCellOperator(r);
}

void AdvectionBackwardEulerFactory::setA(double aX, double aY){
	if(!lhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
	this->aX = aX;
	this->aY = aY;
}
void AdvectionBackwardEulerFactory::setDeltaT(double deltaT){
	if(!lhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
	this->deltaT = deltaT;
}

}
}
