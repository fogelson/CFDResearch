/*
 * CellToFaceIndex.cpp
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
