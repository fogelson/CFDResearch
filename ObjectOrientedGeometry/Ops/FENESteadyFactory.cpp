/*
 * FENESteadyFactory.cpp
 *
 *  Created on: Dec 26, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

FENESteadyFactory::FENESteadyFactory(){
	setD(1);
}
void FENESteadyFactory::setD(double D){
	if(this->D == D)
		return;
	this->D = D;
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENESteadyFactory::setQmax(double Qmax){
	if(uw.getQmax() == Qmax)
		return;
	uw.setQmax(Qmax);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENESteadyFactory::setH(double H){
	if(uw.getH() == H)
		return;
	uw.setH(H);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}
void FENESteadyFactory::setGradU(double u11, double u12, double u21, double u22){
	uw.setGradU(u11,u12,u21,u22);
	if(!lhs.empty() || !rhs.empty()){
		clearMap(lhs);
		clearMap(rhs);
	}
}

void FENESteadyFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}

	CellToCellOperator * L = laplacian.get(g);
	CellToCellOperator * A = uw.get(g);

	CellToCellOperator I;
	I = CellToCellOperator::getIdentity(g);


	CellToCellOperator cL, cR;
	cL = (-D)*(*L);
	cR = (*A);

	CellToCellOperator l = cL - cR;
	CellToCellOperator r = cR;

	lhs[g] = new CellToCellOperator(l);
	rhs[g] = new CellToCellOperator(r);
}

}
}
