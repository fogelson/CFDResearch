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

FaceDoubleArray CellToFaceOperator::apply(CellDoubleArray & u){
	FaceDoubleArray out = g->makeFaceDoubleArray();

	CellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		CellToFaceIndex c2f = (*cIt).first;
		out(c2f.faceIndex) += ((*cIt).second)*u(c2f.i,c2f.j);
	}

	return out;
}

FaceDoubleArray CellToFaceOperator::operator() (CellDoubleArray & u){
	return apply(u);
}

CellDoubleArray FaceToCellOperator::apply(FaceDoubleArray & u){
	CellDoubleArray out = g->makeCellDoubleArray();

	CellToFaceCoefficients::iterator cIt;
	for(cIt = coefficients.begin(); cIt != coefficients.end(); cIt++){
		CellToFaceIndex c2f = (*cIt).first;
		out(c2f.i,c2f.j) += ((*cIt).second)*u(c2f.faceIndex);
	}
	return out;
}

CellDoubleArray FaceToCellOperator::operator() (FaceDoubleArray & u){
	return apply(u);
}

}
}
