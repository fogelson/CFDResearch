/*
 * UpwindFactory.cpp
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

void UpwindFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	CellToFaceOperator F;
	F.setGrid(g);
	vector<Face*>::iterator it;
	for(it = g->faces.begin(); it != g->faces.end(); it++){
		if((*it)->isInterior()){
			double a = 0;
			if((*it)->isNS()){
				a = aY;
			}
			if((*it)->isEW()){
				a = aX;
			}
			Cell * upwind;
			if(a > 0){
				upwind = (*it)->getFrom();
			}
			else{
				upwind = (*it)->getTo();
			}
			FaceIndex faceIndex = (*it)->getIndex();
			CellIndex cellIndex(upwind->getI(), upwind->getJ());
			F.coefficients[faceIndex][cellIndex] = -a;
		}
	}
	Divergence div(g);

	operators[g] = new CellToCellOperator(div(F));
}


void UpwindFactory::setA(double aX, double aY){
	if(this->aX == aX && this->aY == aY)
		return;
	if(!operators.empty()){
		clearMap(operators);
	}
	this->aX = aX;
	this->aY = aY;
}

}
}
