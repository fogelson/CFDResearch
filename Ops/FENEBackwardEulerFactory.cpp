/*
 * FENEBackwardEulerFactory.cpp
 *
 *  Created on: Dec 24, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

FENEBackwardEulerFactory::FENEBackwardEulerFactory(double deltaT){
	setD(1);
	setDeltaT(deltaT);
}
void FENEBackwardEulerFactory::setDeltaT(double deltaT){
	if(this->deltaT == deltaT)
		return;
	this->deltaT = deltaT;
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENEBackwardEulerFactory::setD(double D){
	if(this->D == D)
		return;
	this->D = D;
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENEBackwardEulerFactory::setQmax(double Qmax){
	if(uw.getQmax() == Qmax)
		return;
	uw.setQmax(Qmax);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENEBackwardEulerFactory::setH(double H){
	if(uw.getH() == H)
		return;
	uw.setH(H);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENEBackwardEulerFactory::setGradU(double u11, double u12, double u21, double u22){
	uw.setGradU(u11,u12,u21,u22);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}

void FENEBackwardEulerFactory::produce(Grid * g){
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
