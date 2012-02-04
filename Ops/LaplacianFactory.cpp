/*
 * LaplacianFactory.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

void LaplacianFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	Divergence div(g);
	Gradient grad(g);
	CellToCellOperator L;
	L = div(grad);
	operators[g] = new CellToCellOperator(L);
}

}
}
