/*
 * DiffusionBackwardEulerFactory.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"

//#include "../MexTools/MexTools.h"

#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

DiffusionBackwardEulerFactory::DiffusionBackwardEulerFactory(double D, double deltaT){
	this->D = D;
	this->deltaT = deltaT;
	this->t = 0;
}

void DiffusionBackwardEulerFactory::setTime(double t){
	this->t = t;
	this->clearMap(lhs);
	this->clearMap(rhs);
}

void DiffusionBackwardEulerFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	Divergence div(g);
	Gradient grad(g);
	CellToCellOperator I, L, deltaTDL;
	if(!g->hasFine){
		vector<Face*>::iterator it;
		for(it = g->faces.begin(); it != g->faces.end(); it++){
			if((*it)->isBoundary()){
				int index = (*it)->getIndex();
				TinyVector<double,2> n = (*it)->getNormal();
				Coord c = (*it)->getCentroid();
				double x = c(0), y = c(1);
				TinyVector<double,2> g;
				g(0) = 2*x*exp(t+deltaT);
				g(1) = 2*y*exp(t+deltaT);
				grad.constantTerm(index) = n(0)*g(0) + n(1)*g(1);
			}
		}
	}
/*	OOMexTools::MexPlotTool plotter;
	plotter.newFigure();
	plotter.graphFaceDoubleArray(grad.constantTerm,g);*/
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
