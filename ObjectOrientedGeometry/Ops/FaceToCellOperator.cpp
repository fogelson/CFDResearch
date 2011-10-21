/*
 * FaceToCellOperator.cpp
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

CellDoubleArray FaceToCellOperator::apply(FaceDoubleArray & u){
	CellDoubleArray out = g->makeCellDoubleArray();

	/*CellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		CellToFaceIndex c2f = (*cIt).first;
		out(c2f.i,c2f.j) += ((*cIt).second)*u(c2f.faceIndex);
	}*/

	FaceToCellCoefficients::iterator fIt;
	for(fIt = coefficients.begin(); fIt != coefficients.end(); fIt++){
		int faceIndex = (*fIt).first;
		CellCoefficients::iterator cIt;
		for(cIt = (*fIt).second.begin(); cIt != (*fIt).second.end(); cIt++){
			int i = (*cIt).first.i, j = (*cIt).first.j;
			double c = (*cIt).second;
			out(i,j) += c*u(faceIndex);
		}
	}

	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			out(i,j) += constantTerm(i,j);
		}
	}
	return out;
}

CellDoubleArray FaceToCellOperator::operator() (FaceDoubleArray u){
	return apply(u);
}

/* Returns a cell to cell operator of *this composed with B
 *
 */
CellToCellOperator FaceToCellOperator::operator() (CellToFaceOperator & B){
	CellToCellOperator C;
	C.g = g;
	C.constantTerm.resize(g->xRange,g->yRange);
	C.constantTerm = constantTerm + apply(B.constantTerm);
	CellToFaceCoefficients::iterator fromIt;
	for(fromIt = B.coefficients.begin(); fromIt != B.coefficients.end(); fromIt++){
		int iFrom = (*fromIt).first.i, jFrom = (*fromIt).first.j;
		CellIndex cellFrom(iFrom,jFrom);
		FaceCoefficients::iterator faceIt;
		for(faceIt = (*fromIt).second.begin(); faceIt != (*fromIt).second.end(); faceIt++){
			int faceIndex = (*faceIt).first;
			double cFrom = (*faceIt).second;
			CellCoefficients::iterator toIt;
			for(toIt = coefficients[faceIndex].begin(); toIt != coefficients[faceIndex].end(); toIt++){
				int iTo = (*toIt).first.i, jTo = (*toIt).first.j;
				double cTo = (*toIt).second;
				CellIndex cellTo(iTo,jTo);
				C.coefficients[cellFrom][cellTo] += cFrom*cTo;
			}
		}
	}
	return C;
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

/*int FaceToCellOperator::getRows(){
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
}*/

}
}
