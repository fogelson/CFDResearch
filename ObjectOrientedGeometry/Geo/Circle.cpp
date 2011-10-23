/*
 * Circle.cpp
 *
 *  Created on: Oct 5, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

bool Circle::contains(Cell * c){
	return (contains(c->getFace(N)) && contains(c->getFace(E))
			&& contains(c->getFace(S)) && contains(c->getFace(W)));
}
bool Circle::contains(Face * f){
	return contains(f->getA()) && contains(f->getB());
}
bool Circle::contains(Vertex * v){
	return contains(v->getCoord());
}
bool Circle::contains(Coord c){
	return contains(c(0),c(1));
}
bool Circle::contains(double x, double y){
	return (pow2(x) + pow2(y) < pow2(r));
}

void Circle::grid(){
	double xMin, xMax, yMin, yMax;
	xMin = -r - offset;
	xMax = r;
	yMin = -r - offset;
	yMax = r;

	tile(xMin,yMin,xMax,yMax);

	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			Cell * c = cells(i,j);
			if(contains(c)){
				c->setType(REGULAR);
				c->getFace(N)->setType(REGULAR);
				c->getFace(E)->setType(REGULAR);
				c->getFace(S)->setType(REGULAR);
				c->getFace(W)->setType(REGULAR);
			}
			else{
				c->setType(COVERED);
				// North face
				if(c->getFace(N)->isUncovered()){
					c->setType(IRREGULAR);
				}
				else if(contains(c->getFace(N))){
					c->setType(IRREGULAR);
					c->getFace(N)->setType(REGULAR);
				}
				else if(contains(c->getFace(N)->getA())){
					Face * fN = c->getFace(N);

					c->setType(IRREGULAR);
					fN->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fN->getB()->getCoord();
					newCoord = oldCoord;
					newCoord(0) = sqrt(pow2(r) - pow2(newCoord(1)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fN->setB(newVertex);
				}
				else if(contains(c->getFace(N)->getB())){
					Face * fN = c->getFace(N);

					c->setType(IRREGULAR);
					fN->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fN->getA()->getCoord();
					newCoord = oldCoord;
					newCoord(0) = -sqrt(pow2(r) - pow2(newCoord(1)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fN->setA(newVertex);
				}
				else{
					c->getFace(N)->setType(COVERED);
				}

				// South face
				if(c->getFace(S)->isUncovered()){
					c->setType(IRREGULAR);
				}
				else if(contains(c->getFace(S))){
					Face * fS = c->getFace(S);

					c->setType(IRREGULAR);
					c->getFace(S)->setType(REGULAR);
				}
				else if(contains(c->getFace(S)->getA())){
					Face * fS = c->getFace(S);

					c->setType(IRREGULAR);
					fS->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fS->getB()->getCoord();
					newCoord = oldCoord;
					newCoord(0) = sqrt(pow2(r) - pow2(newCoord(1)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fS->setB(newVertex);
				}
				else if(contains(c->getFace(S)->getB())){
					Face * fS = c->getFace(S);

					c->setType(IRREGULAR);
					fS->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fS->getA()->getCoord();
					newCoord = oldCoord;
					newCoord(0) = -sqrt(pow2(r) - pow2(newCoord(1)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fS->setA(newVertex);
				}
				else{
					c->getFace(S)->setType(COVERED);
				}

				// East face
				if(c->getFace(E)->isUncovered()){
					c->setType(IRREGULAR);
				}
				else if(contains(c->getFace(E))){
					c->setType(IRREGULAR);
					c->getFace(E)->setType(REGULAR);
				}
				else if(contains(c->getFace(E)->getA())){
					Face * fE = c->getFace(E);

					c->setType(IRREGULAR);
					fE->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fE->getB()->getCoord();
					newCoord = oldCoord;
					newCoord(1) = sqrt(pow2(r) - pow2(newCoord(0)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fE->setB(newVertex);
				}
				else if(contains(c->getFace(E)->getB())){
					Face * fE = c->getFace(E);

					c->setType(IRREGULAR);
					fE->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fE->getA()->getCoord();
					newCoord = oldCoord;
					newCoord(1) = -sqrt(pow2(r) - pow2(newCoord(0)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fE->setA(newVertex);
				}
				else{
					c->getFace(E)->setType(COVERED);
				}

				// West face
				if(c->getFace(W)->isUncovered()){
					c->setType(IRREGULAR);
				}
				else if(contains(c->getFace(W))){
					c->setType(IRREGULAR);
					c->getFace(W)->setType(REGULAR);
				}
				else if(contains(c->getFace(W)->getA())){
					Face * fW = c->getFace(W);

					c->setType(IRREGULAR);
					fW->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fW->getB()->getCoord();
					newCoord = oldCoord;
					newCoord(1) = sqrt(pow2(r) - pow2(newCoord(0)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fW->setB(newVertex);
				}
				else if(contains(c->getFace(W)->getB())){
					Face * fW = c->getFace(W);

					c->setType(IRREGULAR);
					fW->setType(IRREGULAR);

					Coord oldCoord, newCoord;
					oldCoord = fW->getA()->getCoord();
					newCoord = oldCoord;
					newCoord(1) = -sqrt(pow2(r) - pow2(newCoord(0)));

					Vertex * newVertex = createVertex();
					newVertex->setCoord(newCoord);

					fW->setA(newVertex);
				}
				else{
					c->getFace(W)->setType(COVERED);
				}
			}
			Face * fB = c->createBoundary();
			addFace(fB);
			if(fB != 0){
				fB->setType(IRREGULAR);
				double xC, yC, x1, x2, y1, y2;
				xC = fB->getCentroid()(0);
				yC = fB->getCentroid()(1);
				x1 = fB->getA()->getCoord()(0);
				y1 = fB->getA()->getCoord()(1);
				x2 = fB->getB()->getCoord()(0);
				y2 = fB->getB()->getCoord()(1);
				TinyVector<double,2> normal;

				double theta = atan2(yC,xC);
				normal(0) = cos(theta);
				normal(1) = sin(theta);

				/*if(xC > 0){
					if(y1 >= y2){
						normal(0) = y1 - y2;
						normal(1) = x2 - x1;
					}
					else{
						normal(0) = y2 - y1;
						normal(1) = x1 - x2;
					}
				}
				else if(xC == 0){
					if(yC > 0){
						normal(0) = 0;
						normal(1) = 1;
					}
					else if(yC < 0){
						normal(0) = 0;
						normal(1) = -1;
					}
					else{
						cout << "Somehow the Circle.cpp grid code thinks there is"
								<< " a boundary cell at the center of the circle."
								<< " Fix that." << endl;
					}
				}
				else{
					if(y1 >= y2){
						normal(0) = y2 - y1;
						normal(1) = x1 - x2;
					}
					else{
						normal(0) = y1 - y2;
						normal(1) = x2 - x1;
					}
				}
				double r = sqrt(pow2(normal(0)) + pow2(normal(1)));
				normal = normal/r;*/
				fB->setNormal(normal);
			}
		}
	}
}

Circle::Circle(double h, double r, double offset){
	setH(h);
	setR(r);
	setOffset(offset);
	grid();
}

void Circle::setR(double r){
	this->r = r;
}

double Circle::getR(){
	return r;
}

void Circle::setOffset(double offset){
	this->offset = offset;
}

double Circle::getOffset(){
	return offset;
}

Grid * Circle::getCoarse(){
	if(hasCoarse){
		return coarseGrid;
	}
	double hC = 2*h;
	coarseGrid = new Circle(hC, r, offset);

	/*for(int iC = coarseGrid->iMin; iC <= coarseGrid->iMax; iC++){
		for(int jC = coarseGrid->jMin; jC <= coarseGrid->jMax; jC++){
			double v = 0;
			v += isUncovered(2*iC,2*jC) ? cells(2*iC,2*jC)->getVolume() : 0;
			v += isUncovered(2*iC-1,2*jC) ? cells(2*iC-1,2*jC)->getVolume() : 0;
			v += isUncovered(2*iC,2*jC-1) ? cells(2*iC,2*jC-1)->getVolume() : 0;
			v += isUncovered(2*iC-1,2*jC-1) ? cells(2*iC-1,2*jC-1)->getVolume() : 0;
			coarseGrid->cells(iC,jC)->setVolume(v);
		}
	}*/
	hasCoarse = true;
	return coarseGrid;
}

}
}
