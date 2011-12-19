/*
 * Face.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

/* Each face contains two vertices. For the case of
 * vertical cells, vertexA should be south of vertexB.
 * For the case of horizontal cells, vertexA should be
 * west of vertexB.
 */

Face::Face(){
	normal = 0;
	vertexA = 0;
	vertexB = 0;
	upToDate = false;
	area = 0;
	fromCell = 0;
	toCell = 0;
	setType(COVERED);
	isBoundaryFace = false;
	isNSFace = false;
	isEWFace = false;
}

void Face::setFrom(Cell * fromCell){
	this->fromCell = fromCell;
}
void Face::setTo(Cell * toCell){
	this->toCell = toCell;
}
Cell * Face::getFrom(){
	return fromCell;
}
Cell * Face::getTo(){
	return toCell;
}

Vertex * Face::getA(){
	return getVertexA();
}

Vertex * Face::getB(){
	return getVertexB();
}

void Face::setA(Vertex * A){
	setVertexA(A);
}

void Face::setB(Vertex * B){
	setVertexB(B);
}

Vertex * Face::getVertexA(){
	return vertexA;
}
Vertex * Face::getVertexB(){
	return vertexB;
}
void Face::setVertexA(Vertex * vertexA){
	if(hasA()){
		//cout << "Face currently has vertex A as " << this->vertexA;
		//cout << " but is now setting it to " << vertexA << endl;
	}
	this->vertexA = vertexA;
	upToDate = false;
}
void Face::setVertexB(Vertex * vertexB){
	if(hasB()){
		//cout << "Face currently has vertex B as " << this->vertexB;
		//cout << " but is now setting it to " << vertexB << endl;
	}
	this->vertexB = vertexB;
	upToDate = false;
}

void Face::update(){
	if(upToDate){
		return;
	}

	double x0, y0, x1, y1;
	x0 = (*vertexA)(0);
	y0 = (*vertexA)(1);
	x1 = (*vertexB)(0);
	y1 = (*vertexB)(1);
	area = sqrt(pow2(x1 - x0) + pow2(y1 - y0));

	centroid(0) = (x0 + x1)/2;
	centroid(1) = (y0 + y1)/2;

	//upToDate = true;
}

Coord Face::getCentroid(){
	update();
	return centroid;
}

double Face::getArea(){
	update();
	return area;
}
bool Face::hasA(){
	return vertexA != 0;
}
bool Face::hasB(){
	return vertexB != 0;
}

void Face::setNormal(TinyVector<double,2> normal){
	this->normal = normal;
}
TinyVector<double,2> Face::getNormal(){
	return normal;
}

bool Face::isBoundary(){
	return isUncovered() && isBoundaryFace;
}

bool Face::isInterior(){
	return isUncovered() && !isBoundaryFace;
}

void Face::setIsBoundary(bool isBoundaryFace){
	this->isBoundaryFace = isBoundaryFace;
}

bool Face::isEW(){
	return isEWFace;
}

bool Face::isNS(){
	return isNSFace;
}

void Face::setIsEW(bool isEWFace){
	this->isEWFace = isEWFace;
}

void Face::setIsNS(bool isNSFace){
	this->isNSFace = isNSFace;
}




}
}
