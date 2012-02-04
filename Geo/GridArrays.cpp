/*
 * GridArrays.cpp
 *
 *  Created on: Oct 4, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

void GridArray::setGrid(Grid * g){
	this->g = g;
}
Grid * GridArray::getGrid(){
	return g;
}



}
}
