/*
 * LaplacianFactory.cpp
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

void LaplacianFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	Divergence div(g);
	Gradient grad(g);
	CellToCellOperator * L = new CellToCellOperator(div(grad));
	operators[g] = L;
}

}
}
