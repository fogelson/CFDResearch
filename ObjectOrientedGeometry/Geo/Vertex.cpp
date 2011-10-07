/*
 * Vertex.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

double & Vertex::operator() (int i){
	return c(i);
}
Coord Vertex::getCoord(){
	return c;
}
void Vertex::setCoord(Coord c){
	this->c = c;
}

}
}
