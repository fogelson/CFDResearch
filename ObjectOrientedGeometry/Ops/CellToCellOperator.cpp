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

CellToCellOperator::CellToCellOperator(){}

CellToCellOperator::CellToCellOperator(const CellToCellOperator & copy){
	g = copy.g;
	coefficients = copy.coefficients;
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = copy.constantTerm;
}

CellToCellOperator & CellToCellOperator::operator= (const CellToCellOperator & rhs){
	g = rhs.g;
	coefficients = rhs.coefficients;
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = rhs.constantTerm;
	return *this;
}

CellDoubleArray CellToCellOperator::apply(CellDoubleArray & u){
	CellDoubleArray out = g->makeCellDoubleArray();

	CellToCellCoefficients::iterator toIt;
	for(toIt = coefficients.begin(); toIt != coefficients.end(); toIt++){
		int iTo = (*toIt).first.i, jTo = (*toIt).first.j;
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			int iFrom = (*fromIt).first.i, jFrom = (*fromIt).first.j;
			double c = (*fromIt).second;
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

CellToCellOperator CellToCellOperator::operator+ (CellToCellOperator & B){
	if(g != B.g){
		cout << "Grid mismatch in operator addition." << endl;
	}
	CellToCellOperator out;
	out = *this;

	CellToCellCoefficients::iterator toIt;
	for(toIt = B.coefficients.begin(); toIt != B.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			out.coefficients[(*toIt).first][(*fromIt).first] += (*fromIt).second;
		}
	}
	out.constantTerm += B.constantTerm;
	return out;
}
CellToCellOperator CellToCellOperator::operator- (CellToCellOperator & B){
	CellToCellOperator out, minusB;
	minusB = (-1.0)*B;
	out = (*this) + minusB;
	return out;
}

CellToCellOperator operator* (CellDoubleArray & a, CellToCellOperator & B){
	CellToCellOperator out;
	out = B;
	CellToCellCoefficients::iterator toIt;
	for(toIt = out.coefficients.begin(); toIt != out.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			(*fromIt).second *= a((*toIt).first.i,(*toIt).first.j);
		}
	}
	out.constantTerm *= a;
	return out;
}

CellToCellOperator operator* (double a, CellToCellOperator & B){
	CellToCellOperator out;
	out = B;
	CellToCellCoefficients::iterator toIt;
	for(toIt = out.coefficients.begin(); toIt != out.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			(*fromIt).second *= a;
		}
	}
	out.constantTerm *= a;
	return out;
}

CellDoubleArray CellToCellOperator::operator() (CellDoubleArray & u){
	return apply(u);
}

CellToCellOperator CellToCellOperator::getIdentity(Grid * g){
	CellToCellOperator out;
	out.g = g;
	out.constantTerm.resize(g->xRange,g->yRange);
	out.constantTerm = 0;
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isUncovered(i,j)){
				CellIndex diag(i,j);
				out.coefficients[diag][diag] = 1;
			}
		}
	}
	return out;
}

}
}
