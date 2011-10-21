/*
 * GSLex.cpp
 *
 *  Created on: Oct 20, 2011
 *      Author: fogelson
 */


#include "Smoothers.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{

void GSLex::smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, CellToCellOperator & C, int its){
	Grid * g = C.g;
	u = u0;
	CellDoubleArray rhs = g->makeCellDoubleArray();
	CellToCellCoefficients::iterator it;
	for(it = C.coefficients.begin(); it != C.coefficients.end(); it++){
		int i = (*it).first.i, j = (*it.first.j);
		CellCoefficients::iterator it2;
		for(it2 = (*it).second.begin(); it2 != (*it).second.end(); it2++){
			if((*it).first != (*it2).first){
				// put in the right place in rhs array
			}
		}
	}
}

}
}
