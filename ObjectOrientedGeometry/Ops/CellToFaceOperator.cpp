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
	FaceDoubleArray outNew = g->makeFaceDoubleArray();

	CellToFaceCoefficients::iterator it;
	for(it = coefficients.begin(); it != coefficients.end(); it++){
		int i = (*it).first.i, j = (*it).first.j, faceIndex = (*it).first.faceIndex;
		double c = (*it).second;
		out(faceIndex) += c*u(i,j);
	}

/*	NewCellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		int i = (*cIt).first.i, j = (*cIt).first.j;
		FaceCoefficients::iterator fIt;
		for(fIt = (*cIt).second.begin(); fIt != (*cIt).second.end(); fIt++){
			int faceIndex = (*fIt).first;
			double c = (*fIt).second;
			//if(g->faces[faceIndex]->isRegular())
				out(faceIndex) += c*u(i,j);
		}
	}*/

	for(int k = 0; k < constantTerm.size(); k++){
		out(k) += constantTerm(k);
	}

	return out;
}

/*FaceDoubleArray CellToFaceOperator::operator() (CellDoubleArray & u){
	return apply(u);
}*/
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

CellToFaceIndex::CellToFaceIndex(int i, int j, int faceIndex){
	this->i = i;
	this->j = j;
	this->faceIndex = faceIndex;
}

bool CellToFaceIndexCompare::operator() (const CellToFaceIndex & lhs, const CellToFaceIndex & rhs) const{
	if(lhs.i != rhs.i){
		return lhs.i < rhs.i;
	}
	else if(lhs.j != rhs.j){
		return lhs.j < rhs.j;
	}
	else{
		return lhs.faceIndex < rhs.faceIndex;
	}
}


}
}
