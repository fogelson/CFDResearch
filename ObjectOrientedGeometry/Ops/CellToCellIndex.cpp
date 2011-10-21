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

CellToCellIndex::CellToCellIndex(int iFrom, int jFrom, int iTo, int jTo){
	this->iFrom = iFrom;
	this->jFrom = jFrom;
	this->iTo = iTo;
	this->jTo = jTo;
}

bool CellToCellIndexCompare::operator() (const CellToCellIndex & lhs, const CellToCellIndex & rhs) const{
	if(lhs.iFrom != rhs.iFrom){
		return lhs.iFrom < rhs.iFrom;
	}
	else if(lhs.jFrom != rhs.jFrom){
		return lhs.jFrom < rhs.jFrom;
	}
	else if(lhs.iTo != rhs.iTo){
		return lhs.iTo < rhs.iTo;
	}
	else{
		return lhs.jTo < rhs.jTo;
	}
}
}
}
