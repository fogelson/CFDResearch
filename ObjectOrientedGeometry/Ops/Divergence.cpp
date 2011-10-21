/*
 * Divergence
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

Divergence::Divergence(Grid * g){
	this->g = g;
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = 0;
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			/*if(false && g->isRegular(i,j)){
				double h = g->getH();
				double V = g->cells(i,j)->getVolume();
				CellToFaceIndex nc2f(i,j,g->cells(i,j)->getFace(N)->getIndex());
				CellToFaceIndex sc2f(i,j,g->cells(i,j)->getFace(S)->getIndex());
				CellToFaceIndex ec2f(i,j,g->cells(i,j)->getFace(E)->getIndex());
				CellToFaceIndex wc2f(i,j,g->cells(i,j)->getFace(W)->getIndex());
				coefficients[nc2f] = (g->cells(i,j)->getFace(N)->getArea())/V;
				coefficients[sc2f] = -(g->cells(i,j)->getFace(N)->getArea())/V;
				coefficients[ec2f] = (g->cells(i,j)->getFace(N)->getArea())/V;
				coefficients[wc2f] = -(g->cells(i,j)->getFace(N)->getArea())/V;
			}*/
			if(g->isUncovered(i,j)){
				double h = g->getH();
				TinyVector<Face*,5> faces = g->cells(i,j)->faces;
				double V = g->cells(i,j)->getVolume();
				if(g->isFaceUncovered(i,j,N)){
					FaceIndex faceIndex = faces(N)->getIndex();
					CellIndex cellIndex(i,j);
					if(g->isFaceRegular(i,j,N)){
						coefficients[faceIndex][cellIndex] = (1/V)*h;
					}
					else{
						coefficients[faceIndex][cellIndex] = (1/V)*faces(N)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,S)){
					FaceIndex faceIndex = faces(S)->getIndex();
					CellIndex cellIndex(i,j);
					if(g->isFaceRegular(i,j,S)){
						coefficients[faceIndex][cellIndex] = -(1/V)*h;
					}
					else{
						coefficients[faceIndex][cellIndex] = -(1/V)*faces(S)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,E)){
					FaceIndex faceIndex = faces(E)->getIndex();
					CellIndex cellIndex(i,j);
					if(g->isFaceRegular(i,j,E)){
						coefficients[faceIndex][cellIndex] = (1/V)*h;
					}
					else{
						coefficients[faceIndex][cellIndex] = (1/V)*faces(E)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,W)){
					FaceIndex faceIndex = faces(W)->getIndex();
					CellIndex cellIndex(i,j);
					if(g->isFaceRegular(i,j,W)){
						coefficients[faceIndex][cellIndex] = -(1/V)*h;
					}
					else{
						coefficients[faceIndex][cellIndex] = -(1/V)*faces(W)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,B)){
					FaceIndex faceIndex = faces(B)->getIndex();
					CellIndex cellIndex(i,j);
					coefficients[faceIndex][cellIndex] = faces(B)->getArea()/V;
				}
			}
		}
	}
}

}
}
