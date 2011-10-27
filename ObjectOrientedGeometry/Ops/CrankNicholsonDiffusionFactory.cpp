/*
 * CrankNicholsonDiffusionFactory
 *
 *  Created on: Oct 25, 2011
 *      Author: fogelson
 */


#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

CrankNicholsonDiffusionFactory::CrankNicholsonDiffusionFactory(double deltaT){
	this->deltaT = deltaT;
}

void CrankNicholsonDiffusionFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}

	Divergence div(g);
	//cout << "Created divergence" << endl;
	Gradient grad(g);
	//cout << "Created gradient" << endl;

	CellToCellOperator I, L, DL;
	CellToCellOperator lhs, rhs;
	I = CellToCellOperator::getIdentity(g);
	//cout << "Created identity" << endl;
	L = div(grad);
	DL = (deltaT/2)*L;
	lhs = I - DL;
	rhs = I + DL;
	//cout << "Created lhs and rhs" << endl;

	lhsOperators[g] = new CellToCellOperator(lhs);
	rhsOperators[g] = new CellToCellOperator(rhs);
	//cout << "Loaded lhs and rhs into the operator maps" << endl;
	//cout << "Created new Crank Nicholson operators for the grid with h = " << g->h << endl;
}

}
}
