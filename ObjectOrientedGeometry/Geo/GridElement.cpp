/*
 * GridElement.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

void GridElement::setIndex(int index){
	this->index = index;
}
int GridElement::getIndex(){
	return index;
}
bool GridElement::isRegular(){
	return (type == REGULAR);
}
bool GridElement::isIrregular(){
	return (type == IRREGULAR);
}
bool GridElement::isCovered(){
	return (type == COVERED);
}
bool GridElement::isUncovered(){
	return (type != COVERED);
}
void GridElement::setType(Type type){
	this->type = type;
}
Type GridElement::getType(){
	return type;
}

}
}
