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
	setA(aX,aY);
	setD(D);
	setDeltaT(deltaT);
}

void AdvectionDiffusionBackwardEulerFactory::setA(double aX, double aY){
	if(this->aX == aX && this->aY == aY)
		return;
	this->aX = aX;
	this->aY = aY;
	uw.setA(aX,aY);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}

void AdvectionDiffusionBackwardEulerFactory::setD(double D){
	if(this->D == D)
		return;
	this->D = D;
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}

void AdvectionDiffusionBackwardEulerFactory::setDeltaT(double deltaT){
	if(this->deltaT == deltaT)
		return;
	this->deltaT = deltaT;
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}

void AdvectionDiffusionBackwardEulerFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}

	CellToCellOperator * L = laplacian.get(g);
	CellToCellOperator * A = uw.get(g);

	CellToCellOperator I;
	I = CellToCellOperator::getIdentity(g);


	CellToCellOperator cL, cR;
	cL = deltaT*D*(*L);
	cR = deltaT*(*A);

	CellToCellOperator l = I - cL - cR;
	CellToCellOperator r = I;

	lhs[g] = new CellToCellOperator(l);
	rhs[g] = new CellToCellOperator(r);
}

}
}
