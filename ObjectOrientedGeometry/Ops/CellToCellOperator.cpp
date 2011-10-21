/*
 * CellToCell
 *
 *  Created on: Oct 20, 2011
 *      Author: fogelson
 */

#include "Operators.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

CellToCellOperator & CellToCellOperator::operator= (const CellToCellOperator & rhs){
	g = rhs.g;
	coefficients = rhs.coefficients;
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = rhs.constantTerm;
	return *this;
}

CellDoubleArray CellToCellOperator::apply(CellDoubleArray & u){
	CellDoubleArray out = g->makeCellDoubleArray();

	CellToCellCoefficients::iterator fromIt;
	for(fromIt = coefficients.begin(); fromIt != coefficients.end(); fromIt++){
		int iFrom = (*fromIt).first.i, jFrom = (*fromIt).first.j;
		CellCoefficients::iterator toIt;
		for(toIt = (*fromIt).second.begin(); toIt != (*fromIt).second.end(); toIt++){
			int iTo = (*toIt).first.i, jTo = (*toIt).first.j;
			double c = (*toIt).second;
			out(iTo,jTo) += c*u(iFrom,jFrom);
		}
	}
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			out(i,j) += constantTerm(i,j);
		}
	}
	return out;
}

CellDoubleArray CellToCellOperator::operator() (CellDoubleArray & u){
	return apply(u);
}


}
}
