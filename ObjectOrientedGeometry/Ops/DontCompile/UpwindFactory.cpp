/*
 * UpwindFactory.cpp
 *
 *  Created on: Oct 27, 2011
 *      Author: fogelson
 */


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
	Divergence div(g);
	CellToFaceOperator F;
	F.setGrid(g);
	vector<Face*>::iterator it;
	for(it = g->faces.begin(); it != g->faces.end(); it++){
		if((*it) != 0 && (*it)->isBoundary()){
			FaceIndex faceIndex = (*it)->getIndex();
			Coord pos = (*it)->getCentroid();
			TinyVector<double,2> aSpeed = (*c)(pos);
			double a = sum(aSpeed*((*it)->getNormal()));
			if(a > 0){
				CellIndex upwind((*it)->getFrom()->getI(),(*it)->getFrom()->getJ());
				F.coefficients[faceIndex][upwind] = -a;
			}
			if(a < 0){
				F.constantTerm(faceIndex) = -a;
			}
		}
		else if((*it) != 0 && (*it)->isInterior()){
			FaceIndex faceIndex = (*it)->getIndex();
			Coord pos = (*it)->getCentroid();
			TinyVector<double,2> aSpeed = (*c)(pos);
			double a = sum(aSpeed*((*it)->getNormal()));
			if(a > 0){
				CellIndex upwind((*it)->getFrom()->getI(),(*it)->getFrom()->getJ());
				F.coefficients[faceIndex][upwind] = -a;
			}
			else if(a < 0){
				//cout << "Negative speed through face" << endl;
				CellIndex upwind((*it)->getTo()->getI(),(*it)->getTo()->getJ());
				F.coefficients[faceIndex][upwind] = -a;
			}
		}
	}
	operators[g] = new CellToCellOperator(div(F));
}

UpwindFactory::UpwindFactory(TinyVector<double,2> (* c)(Coord pos)){
	this->c = c;
}

}
}
