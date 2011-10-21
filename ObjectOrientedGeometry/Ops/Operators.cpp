/*
 * Operators.cpp
 *
 *  Created on: Oct 12, 2011
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

	CellToCellCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		CellToCellIndex c2c = (*cIt).first;
		double c = (*cIt).second;
		out(c2c.iTo,c2c.jTo) += c*u(c2c.iFrom,c2c.jFrom);
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

CellDoubleArray FaceToCellOperator::apply(FaceDoubleArray & u){
	CellDoubleArray out = g->makeCellDoubleArray();

	CellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		CellToFaceIndex c2f = (*cIt).first;
		out(c2f.i,c2f.j) += ((*cIt).second)*u(c2f.faceIndex);
	}
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			out(i,j) += constantTerm(i,j);
		}
	}
	return out;
}

/*CellDoubleArray FaceToCellOperator::operator() (FaceDoubleArray & u){
	return apply(u);
}*/
CellDoubleArray FaceToCellOperator::operator() (FaceDoubleArray u){
	return apply(u);
}

/*
 * Outputs an operator C, such that C(u) = A(B(u)), where A is the
 * current FaceToCell operator and B is a CellToFace operator.
 */
/*CellToCellOperator FaceToCellOperator::operator() (CellToFaceOperator & B){
	CellToCellOperator C;
	C.g = g;
	C.constantTerm.resize(g->xRange,g->yRange);
	C.constantTerm = constantTerm + apply(B.constantTerm);
	CellToFaceCoefficients::iterator Bit;
	for(Bit = B.coefficients.begin(); Bit != B.coefficients.end(); Bit++){
		int iFrom = (*Bit).first.i;
		int jFrom = (*Bit).first.j;
		int faceIndex = (*Bit).first.faceIndex;
		CellToFaceCoefficients::iterator Ait;
		for(Ait = coefficients.begin(); Ait != coefficients.end(); Ait++){
			if((*Ait).first.faceIndex == faceIndex){
				int iTo = (*Ait).first.i;
				int jTo = (*Ait).first.j;
				CellToCellIndex ind(iFrom,jFrom,iTo,jTo);
				//if(C.coefficients.count(ind) > 0){
					C.coefficients[ind] += ((*Ait).second)*((*Bit).second);
				//}
				//else{
				//	C.coefficients[ind] = ((*Ait).second)*((*Bit).second);
				//}
			}
		}

	}
	return C;
}*/

int FaceToCellOperator::getRows(){
	return g->cells.size();
}

int FaceToCellOperator::getCols(){
	return g->faces.size();
}

int FaceToCellOperator::getRowIndex(CellToFaceIndex ind){
	return (ind.i - g->iMin) + (ind.j - g->jMin)*(g->cells.length(0));
}
int FaceToCellOperator::getColIndex(CellToFaceIndex ind){
	return ind.faceIndex;
}

Array<double,2> FaceToCellOperator::getMatrix(){
	Array<double,2> A;
	A.resize(getRows(),getCols());
	A = 0;

	CellToFaceCoefficients::iterator cIt;

	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		A(getRowIndex((*cIt).first),getColIndex((*cIt).first)) = (*cIt).second;
	}
	return A;
}

Array<double,1> FaceToCellOperator::getVector(){
	Array<double,1> b;
	b.resize(getRows());
	b = 0;

	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			CellToFaceIndex c2f(i,j,-1);
			b(getRowIndex(c2f)) = constantTerm(i,j);
		}
	}

	return b;
}

}
}
