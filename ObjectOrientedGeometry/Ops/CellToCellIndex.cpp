/*
 * CellToCellIndex.cpp
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

CellToCellIndex::CellToCellIndex(int i1, int j1, int i2, int j2){
	this->i1 = i1;
	this->j1 = j1;
	this->i2 = i2;
	this->j2 = j2;
}

bool CellToCellIndexCompare::operator() (const CellToCellIndex & lhs, const CellToCellIndex & rhs) const{
	if(lhs.i1 != rhs.i1){
		return lhs.i1 < rhs.i1;
	}
	else if(lhs.j1 != rhs.j1){
		return lhs.j1 < rhs.j1;
	}
	if(lhs.i2 != rhs.i2){
		return lhs.i2 < rhs.i2;
	}
	else if(lhs.j2 != rhs.j2){
		return lhs.j2 < rhs.j2;
	}
}
}
}
