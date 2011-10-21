/*
 * CellToFaceOperator.cpp
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


FaceDoubleArray CellToFaceOperator::apply(CellDoubleArray & u){
	FaceDoubleArray out = g->makeFaceDoubleArray();

	CellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		int i = (*cIt).first.i, j = (*cIt).first.j;
		FaceCoefficients::iterator fIt;
		for(fIt = (*cIt).second.begin(); fIt != (*cIt).second.end(); fIt++){
			int faceIndex = (*fIt).first;
			double c = (*fIt).second;
			out(faceIndex) += c*u(i,j);
		}
	}

	for(int k = 0; k < constantTerm.size(); k++){
		out(k) += constantTerm(k);
	}

	return out;
}

FaceDoubleArray CellToFaceOperator::operator() (CellDoubleArray u){
	return apply(u);
}

/*int CellToFaceOperator::getCols(){
	return g->cells.size();
}

int CellToFaceOperator::getRows(){
	return g->faces.size();
}

int CellToFaceOperator::getColIndex(CellToFaceIndex ind){
	return (ind.i - g->iMin) + (ind.j - g->jMin)*(g->cells.length(0));
}
int CellToFaceOperator::getRowIndex(CellToFaceIndex ind){
	return ind.faceIndex;
}

Array<double,2> CellToFaceOperator::getMatrix(){
	Array<double,2> A;
	A.resize(getRows(),getCols());
	A = 0;

	cout << getRows() << endl;
	cout << getCols() << endl;

	CellToFaceCoefficients::iterator cIt;

	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		A(getRowIndex((*cIt).first),getColIndex((*cIt).first)) = (*cIt).second;
	}
	return A;
}

Array<double,1> CellToFaceOperator::getVector(){
	Array<double,1> b;
	b.resize(getRows());
	b = 0;

	b = constantTerm;

	return b;
}*/

}
}
