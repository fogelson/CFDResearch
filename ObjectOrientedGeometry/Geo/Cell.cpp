/*
 * Cell.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

Cell::Cell(){
	faces = 0;
}

void Cell::update(){
	if(upToDate){
		return;
	}

	numberOfFaces = 0;

	vertices.clear();

	if(faces(N) != 0 && faces(N)->isUncovered()){
		numberOfFaces++;
		Vertex * A, * B;
		A = (faces(N))->getVertexA(); // NW
		B = (faces(N))->getVertexB(); // NE
		if(vertices.empty()){
			vertices.push_back(A);
			vertices.push_back(B);
		}
		else{
			if(A != vertices.back()){
				vertices.push_back(A);
			}
			if(B != vertices.front()){
				vertices.push_back(B);
			}
		}
	}
	if(faces(E) != 0 && faces(E)->isUncovered()){
		numberOfFaces++;
		Vertex * A, * B;
		A = (faces(E))->getVertexA(); // SE
		B = (faces(E))->getVertexB(); // NE
		if(vertices.empty()){
			vertices.push_back(B);
			vertices.push_back(A);
		}
		else{
			if(B != vertices.back()){
				vertices.push_back(B);
			}
			if(A != vertices.front()){
				vertices.push_back(A);
			}
		}
	}
	if(faces(S) != 0 && faces(S)->isUncovered()){
		numberOfFaces++;
		Vertex * A, * B;
		A = (faces(S))->getVertexA();
		B = (faces(S))->getVertexB();
		if(vertices.empty()){
			vertices.push_back(B);
			vertices.push_back(A);
		}
		else{
			if(B != vertices.back()){
				vertices.push_back(B);
			}
			if(A != vertices.front()){
				vertices.push_back(A);
			}
		}
	}
	if(faces(W) != 0 && faces(W)->isUncovered()){
		numberOfFaces++;
		Vertex * A, * B;
		A = (faces(W))->getVertexA();
		B = (faces(W))->getVertexB();
		if(vertices.empty()){
			vertices.push_back(A);
			vertices.push_back(B);
		}
		else{
			if(A != vertices.back()){
				vertices.push_back(A);
			}
			if(B != vertices.front()){
				vertices.push_back(B);
			}
		}
	}
	if(faces(B) != 0){
		numberOfFaces++;
	}

	/*vector<Vertex*>::iterator it;
	cout << "Here are the vertices: " << endl;
	for(it = vertices.begin(); it != vertices.end(); it++){
		cout << (*it)->getCoord() << ". " << (*it) << endl;
	}*/
/*
	for(int k = 0; k < faces.length()-1; k++){
		if(faces(k) != 0){
			numberOfFaces++;
			Vertex * A, * B;
			A = (faces(k))->getVertexA();
			B = (faces(k))->getVertexB();
			if(vertices.empty()){
				vertices.push_back(A);
				vertices.push_back(B);
			}
			else{
				if(A != vertices.back()){
					vertices.push_back(A);
				}
				if(B != vertices.front()){
					vertices.push_back(B);
				}
			}
		}
	}*/
	numberOfVertices = vertices.size();

	double signedVolume = 0;
	for(int k = 0; k < vertices.size(); k++){
		int kP = (k+1) % vertices.size();
		Coord c, cP;
		c = vertices[k]->getCoord();
		cP = vertices[kP]->getCoord();
		signedVolume += c(0)*cP(1) - cP(0)*c(1);
	}
	signedVolume /= 2;
	volume = abs(signedVolume);

	centroid = 0;
	for(int k = 0; k < vertices.size(); k++){
		int kP = (k+1) % vertices.size();
		Coord c, cP;
		c = vertices[k]->getCoord();
		cP = vertices[kP]->getCoord();
		double temp = c(0)*cP(1) - cP(0)*c(1);
		centroid(0) += (c(0)+cP(0))*temp;
		centroid(1) += (c(1)+cP(1))*temp;
	}
	centroid /= 6*signedVolume;

	//upToDate = true;
}

Coord Cell::getCentroid(){
	update();
	return centroid;
}

double Cell::getVolume(){
	update();
	return volume;
}

double Cell::getVolumeFraction(){
	update();
	double h = g->getH();
	return getVolume()/pow2(h);
}
int Cell::getNumberOfFaces(){
	update();
	return numberOfFaces;
}
int Cell::getNumberOfVertices(){
	update();
	return numberOfVertices;
}

vector<Vertex*> Cell::getVertices(){
	update();
	return vertices;
}

Face * Cell::operator() (Direction d){
	return faces(d);
}
void Cell::setFace(Face * f, Direction d){
	faces(d) = f;
	upToDate = false;
}
bool Cell::hasFace(Direction d){
	return faces(d) != 0;
}
Face * Cell::getFace(Direction d){
	return faces(d);
}
void Cell::setCenter(Coord c){
	center = c;
}
Coord Cell::getCenter(){
	return center;
}

/* If the cell doesn't have a boundary pointer and
 * needs one, this method creates a new one and
 * returns it. Deletion is up to whoever calls it.
 *
 * If the cell already has a boundary pointer or
 * doesn't need one, this method returns a 0
 * pointer.
 */
Face * Cell::createBoundary(){
	upToDate = false;
	update();
	if(faces(B) != 0){
		return 0;
	}
	Face * fB = 0;

	vector<Face*> uncoveredFaces;
	if(faces(N) != 0 && faces(N)->isUncovered()){
		uncoveredFaces.push_back(faces(N));
	}
	if(faces(E) != 0 && faces(E)->isUncovered()){
		uncoveredFaces.push_back(faces(E));
	}
	if(faces(S) != 0 && faces(S)->isUncovered()){
		uncoveredFaces.push_back(faces(S));
	}
	if(faces(W) != 0 && faces(W)->isUncovered()){
		uncoveredFaces.push_back(faces(W));
	}
	if(uncoveredFaces.size() == 0){
		return 0;
	}

	vector<Vertex*>::iterator v;
	vector<Face*>::iterator f;
	Vertex * vA, * vB;
	vA = 0; vB = 0;
	for(v = vertices.begin(); v != vertices.end(); v++){
		int count = 0;
		for(f = uncoveredFaces.begin(); f != uncoveredFaces.end(); f++){
			if(*v == (*f)->getA() || *v == (*f)->getB()){
				count++;
			}
		}
		if(count == 1 && vA == 0){
			vA = *v;
			//cout << "Set vA to " << vA->getCoord();
		}
		else if(count == 1 && vB == 0){
			vB = *v;
			//cout << "and set vB to " << vB->getCoord() << endl;
		}
	}
	if(vA != 0 && vB != 0){
		fB = new Face();
		fB->setA(vA);
		fB->setB(vB);
		fB->setType(IRREGULAR);
		faces(B) = fB;
		//cout << "Created boundary face connecting " << vA->getCoord() << " and " << vB->getCoord() << endl;
	}
	return fB;

/*	for(int k = 0; k < uncoveredFaces.size(); k++){
		int kP = (k + 1) % uncoveredFaces.size();
		Face * f1, * f2;
		f1 = uncoveredFaces[k];
		f2 = uncoveredFaces[kP];
		if(		   f1->getA() != f2->getA()
				&& f1->getA() != f2->getB()
				&& f1->getB() != f2->getA()
				&& f1->getB() != f2->getB()){

			Vertex * vA, * vB;
			vA = 0;
			vB = 0;

			vector<Vertex*> v;
			for(v = vertices.begin(); v != vertices.end(); v++){
				if(vA == 0){
					if(f1->getA() == *v){
						vA = *v;
					}
					else if(f1->getB() == *v){
						vA = *v;
					}
					else if(f2->getA() == *v){
						vA = *v;
					}
					else if(f2->getB() == *v){
						vA = *v;
					}
				}
			}
			for(v = vertices.begin(); v != vertices.end(); v++){
				if(vB == 0){
					if(f1->getA() == *v){
						vB = *v;
					}
					else if(f1->getB() == *v){
						vB = *v;
					}
					else if(f2->getA() == *v){
						vB = *v;
					}
					else if(f2->getB() == *v){
						vB = *v;
					}
				}
			}
		}
	}*/
	return fB;
}

}
}

