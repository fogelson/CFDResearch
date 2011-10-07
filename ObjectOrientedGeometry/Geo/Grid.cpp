/*
 * Grid.cpp
 *
 *  Created on: Oct 4, 2011
 *      Author: fogelson
 */


#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

Grid::~Grid(){
	//cout << "Called destructor." << endl;
	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			delete cells(i,j);
		}
	}
	//cout << "Delete cells." << endl;
	vector<Face*>::iterator faceIt;
	for(faceIt = faces.begin(); faceIt != faces.end(); faceIt++){
		delete (*faceIt);
	}
	//cout << "Deleted faces." << endl;
	vector<Vertex*>::iterator vertexIt;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); vertexIt++){
		delete (*vertexIt);
	}
	//cout << "Deleted vertices." << endl;
	if(hasCoarse){
		delete coarseGrid;
	}
	if(hasFine){
		fineGrid->setHasCoarse(false);
	}
	//cout << "End of destructor" << endl;
}

Grid::Grid(){
	hasCoarse = false;
	hasFine = false;
}
void Grid::setHasCoarse(bool hasCoarse){
	this->hasCoarse = hasCoarse;
}

void Grid::tile(double xMin, double yMin, double xMax, double yMax){
	/*iMin = 1;
	jMin = 1;
	iMax = 4;
	jMax = 5;
	xRange.setRange(iMin,iMax,1);
	yRange.setRange(jMin,jMax,1);
	return;*/

	double xWidth = xMax - xMin;
	double yWidth = yMax - yMin;
	int xCells = ceil(xWidth/h);
	int yCells = ceil(yWidth/h);
	xCells = xCells + (xCells % 2);
	yCells = yCells + (yCells % 2);

	iMin = 1;
	jMin = 1;
	iMax = xCells;
	jMax = yCells;

	xRange.setRange(iMin,iMax,1);
	yRange.setRange(jMin,jMax,1);

	cells.resize(xRange,yRange);
	cells = 0;

	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			Coord c;
			c(0) = xMin + i*h - h/2;
			c(1) = yMin + j*h - h/2;
			Cell * current = createCell(i,j);
			current->setCenter(c);
		}
	}

	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			Cell * current = cells(i,j);
			if(!current->hasFace(N)){
				Face * newFace = createFace();
				current->setFace(newFace,N);
				if(j+1 <= jMax){
					cells(i,j+1)->setFace(newFace,S);
				}
			}
			if(!current->hasFace(S)){
				Face * newFace = createFace();
				current->setFace(newFace,S);
				if(j-1 >= jMin){
					cells(i,j-1)->setFace(newFace,N);
				}
			}
			if(!current->hasFace(E)){
				Face * newFace = createFace();
				current->setFace(newFace,E);
				if(i+1 <= iMax){
					cells(i+1,j)->setFace(newFace,W);
				}
			}
			if(!current->hasFace(W)){
				Face * newFace = createFace();
				current->setFace(newFace,W);
				if(i-1 >= iMin){
					cells(i-1,j)->setFace(newFace,E);
				}
			}
		}
	}

	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			if(cells(i,j)->getFace(N)->hasA()){
				cells(i,j)->getFace(W)->setB(cells(i,j)->getFace(N)->getA());
			}
			if(cells(i,j)->getFace(N)->hasB()){
				cells(i,j)->getFace(E)->setB(cells(i,j)->getFace(N)->getB());
			}
			if(cells(i,j)->getFace(E)->hasA()){
				cells(i,j)->getFace(S)->setB(cells(i,j)->getFace(E)->getA());
			}
			if(cells(i,j)->getFace(E)->hasB()){
				cells(i,j)->getFace(N)->setB(cells(i,j)->getFace(E)->getB());
			}
			if(cells(i,j)->getFace(S)->hasA()){
				cells(i,j)->getFace(W)->setA(cells(i,j)->getFace(S)->getA());
			}
			if(cells(i,j)->getFace(S)->hasB()){
				cells(i,j)->getFace(E)->setA(cells(i,j)->getFace(S)->getB());
			}
			if(cells(i,j)->getFace(W)->hasA()){
				cells(i,j)->getFace(S)->setA(cells(i,j)->getFace(W)->getA());
			}
			if(cells(i,j)->getFace(W)->hasB()){
				cells(i,j)->getFace(N)->setA(cells(i,j)->getFace(W)->getB());
			}

			Coord NE, NW, SE, SW;
			NE(0) = xMin + i*h;
			NE(1) = yMin + j*h;
			NW(0) = xMin + (i-1)*h;
			NW(1) = yMin + j*h;
			SE(0) = xMin + i*h;
			SE(1) = yMin + (j-1)*h;
			SW(0) = xMin + (i-1)*h;
			SW(1) = yMin + (j-1)*h;

			Cell * current = cells(i,j);
			//cout << current->getFace(N)->hasA() << endl;


			if(!current->getFace(N)->hasA()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(NW);
				current->getFace(N)->setVertexA(newVertex);
				current->getFace(W)->setVertexB(newVertex);
			}
			if(!current->getFace(N)->hasB()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(NE);
				current->getFace(N)->setVertexB(newVertex);
				current->getFace(E)->setVertexB(newVertex);
			}
			if(!current->getFace(E)->hasA()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(SE);
				current->getFace(S)->setVertexB(newVertex);
				current->getFace(E)->setVertexA(newVertex);
			}
			if(!current->getFace(E)->hasB()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(NE);
				current->getFace(E)->setVertexB(newVertex);
				current->getFace(N)->setVertexB(newVertex);
			}
			if(!current->getFace(S)->hasA()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(SW);
				current->getFace(W)->setVertexA(newVertex);
				current->getFace(S)->setVertexA(newVertex);
			}
			if(!current->getFace(S)->hasB()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(SE);
				current->getFace(S)->setVertexB(newVertex);
				current->getFace(E)->setVertexA(newVertex);
			}
			if(!current->getFace(W)->hasA()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(SW);
				current->getFace(S)->setVertexA(newVertex);
				current->getFace(W)->setVertexA(newVertex);
			}
			if(!current->getFace(W)->hasB()){
				Vertex * newVertex = createVertex();
				newVertex->setCoord(NW);
				current->getFace(N)->setVertexA(newVertex);
				current->getFace(W)->setVertexB(newVertex);
			}
		}
	}
	//cout << "Created vertices" << endl;

}

double Grid::getVolume(){
	double v = 0;
	for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			if(cells(i,j)->isUncovered()){
				v += cells(i,j)->getVolume();
			}
		}
	}
}

double Grid::getH(){
	return h;
}
void Grid::setH(double h){
	this->h = h;
}

Cell * Grid::createCell(int i, int j){
	Cell * out = new Cell();
	cells(i,j) = out;
	return out;
}

Face * Grid::createFace(){
	Face * out = new Face();
	out->setIndex(faces.size());
	faces.push_back(out);
	return out;
}

Vertex * Grid::createVertex(){
	Vertex * out = new Vertex();
	out->setIndex(vertices.size());
	vertices.push_back(out);
	return out;
}
void Grid::addCell(Cell * c, int i, int j){
	if(c != 0){
		if(cells(i,j) != 0){
			delete cells(i,j);
		}
		cells(i,j) = c;
	}
}
void Grid::addFace(Face * f){
	if(f != 0){
		f->setIndex(faces.size());
		faces.push_back(f);
	}
}
void Grid::addVertex(Vertex * v){
	if(v != 0){
		v->setIndex(vertices.size());
		vertices.push_back(v);
	}
}

CellDoubleArray Grid::makeCellDoubleArray(){
	CellDoubleArray out;
	return out;/*
	CellDoubleArray out;
	out.resize(xRange, yRange);
	out = 0;
	/*for(int i = iMin; i <= iMax; i++){
		for(int j = jMin; j <= jMax; j++){
			out(i,j) = i + j;
		}
	}*/
	return out;
}

FaceDoubleArray Grid::makeFaceDoubleArray(){
	FaceDoubleArray out;
	return out;
	/*faceRange.setRange(0,faces.size()-1);
	FaceDoubleArray * out;
	out = new FaceDoubleArray();
	out->resize(4);
	(*out) = 0;
	return *out;*/
}

VertexDoubleArray Grid::makeVertexDoubleArray(){
	VertexDoubleArray out;
	return out;/*
	vertexRange.setRange(0,vertices.size()-1);
	VertexDoubleArray out;
	out.resize(vertexRange);
	out = 0;
	return out;*/
}
void Grid::resizeCellDoubleArray(CellDoubleArray & c){
	c.resize(xRange,yRange);
	c = 0;
}
void Grid::resizeFaceDoubleArray(FaceDoubleArray & f){
	f.resize(faces.size());
	//f.resize(400000000);
	f = 0;
}
void Grid::resizeVertexDoubleArray(VertexDoubleArray & v){
	v.resize(vertices.size());
	v = 0;
}

}
}
