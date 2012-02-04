/*
 * CellIndex.cpp
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

CellIndex::CellIndex(int i, int j){
	this->i = i;
	this->j = j;
}

/*bool CellIndexCompare::operator() (const CellIndex & lhs, const CellIndex & rhs) const{
	if(lhs.i != rhs.i){
		return lhs.i < rhs.i;
	}
	else{
		return lhs.j < rhs.j;
	}
}*/

bool CellIndex::operator== (const CellIndex & rhs) const{
	return (i == rhs.i && j == rhs.j);
}

bool CellIndex::operator!= (const CellIndex & rhs) const{
	return !(*this == rhs);
}

bool CellIndex::operator< (const CellIndex & rhs) const{
	if(i != rhs.i){
		return i < rhs.i;
	}
	else{
		return j < rhs.j;
	}
}

bool CellIndex::operator> (const CellIndex & rhs) const{
	return !(*this == rhs) && !(*this < rhs);
}

bool CellIndex::operator<= (const CellIndex & rhs) const{
	return (*this < rhs) || (*this == rhs);
}

bool CellIndex::operator>= (const CellIndex & rhs) const{
	return (*this > rhs) || (*this == rhs);
}

}
}
