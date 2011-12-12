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

void CellToFaceOperator::setGrid(Grid * g){
	this->g = g;
	constantTerm.resize(g->faces.size());
	constantTerm = 0;
}

FaceDoubleArray CellToFaceOperator::apply(CellDoubleArray & u){
	FaceDoubleArray out = g->makeFaceDoubleArray();

	CellToFaceCoefficients::iterator fIt;
	for(fIt = coefficients.begin(); fIt != coefficients.end(); fIt++){
		int faceIndex = (*fIt).first;
		CellCoefficients::iterator cIt;
		for(cIt = (*fIt).second.begin(); cIt != (*fIt).second.end(); cIt++){
			int i = (*cIt).first.i, j = (*cIt).first.j;
			double c = (*cIt).second;
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

CellToFaceOperator & CellToFaceOperator::operator= (const CellToFaceOperator & rhs){
	g = rhs.g;
	coefficients = rhs.coefficients;
	constantTerm.resize(rhs.constantTerm.size());
	constantTerm = rhs.constantTerm;
	return *this;
}
CellToFaceOperator CellToFaceOperator::operator+ (CellToFaceOperator & B){
	if(g != B.g){
		cout << "Grid mismatch in operator addition." << endl;
	}
	CellToFaceOperator out;
	out = *this;

	CellToFaceCoefficients::iterator toIt;
	for(toIt = B.coefficients.begin(); toIt != B.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			out.coefficients[(*toIt).first][(*fromIt).first] += (*fromIt).second;
		}
	}
	out.constantTerm += B.constantTerm;
	return out;

}
CellToFaceOperator CellToFaceOperator::operator- (CellToFaceOperator & B){
	CellToFaceOperator out, minusB;
	minusB = (-1.0)*B;
	out = (*this) + minusB;
	return out;
}

CellToFaceOperator operator* (FaceDoubleArray & a, CellToFaceOperator & B){
	CellToFaceOperator out;
	out = B;
	CellToFaceCoefficients::iterator toIt;
	for(toIt = out.coefficients.begin(); toIt != out.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			(*fromIt).second *= a((*toIt).first);
		}
	}
	out.constantTerm *= a;
	return out;
}
CellToFaceOperator operator* (double a, CellToFaceOperator & B){
	CellToFaceOperator out;
	out = B;
	CellToFaceCoefficients::iterator toIt;
	for(toIt = out.coefficients.begin(); toIt != out.coefficients.end(); toIt++){
		CellCoefficients::iterator fromIt;
		for(fromIt = (*toIt).second.begin(); fromIt != (*toIt).second.end(); fromIt++){
			(*fromIt).second *= a;
		}
	}
	out.constantTerm *= a;
	return out;
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
