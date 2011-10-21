/*
 * Gradient.cpp
 *
 *  Created on: Oct 19, 2011
 *      Author: fogelson
 */

/*
#include "Operators.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

ConstantAdvectiveFlux::ConstantAdvectiveFlux(Grid * g, double aX, double aY){
	this->g = g;
	this->aX = aX;
	this->aY = aY;
	double h = g->getH();
	constantTerm.resize(g->faces.size());
	constantTerm = 0;
	for(int k = 0; k < g->faces.size(); k++){
		if(g->faces[k] != 0 && g->faces[k]->isUncovered()){
			double x, y, nX, nY;
			x = g->faces[k]->getCentroid()(0);
			y = g->faces[k]->getCentroid()(1);
			nX = g->faces[k]->getNormal()(0);
			nY = g->faces[k]->getNormal()(1);
			double a = aX*nX + aY*nY;

			Cell * from, * to;
			from = g->faces[k]->getFrom();
			to = g->faces[k]->getTo();
			if(a > 0 && from != 0 && from->isUncovered()){
				int i = from->getI(), j = from->getJ();
				CellToFaceIndex ind(i,j,k);
			}
			if(a < 0 && to != 0 && to->isUncovered()){
				int i = to->getI(), j = to->getJ();
				CellToFaceIndex ind(i,j,k);
				coefficients[ind] = a;
			}
		}
	}
}

}
}*/
