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
					int faceIndex = faces(N)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					if(g->isFaceRegular(i,j,N)){
						coefficients[ind] = (1/V)*h;
					}
					else{
						coefficients[ind] = (1/V)*faces(N)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,S)){
					int faceIndex = faces(S)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					if(g->isFaceRegular(i,j,S)){
						coefficients[ind] = -(1/V)*h;
					}
					else{
						coefficients[ind] = -(1/V)*faces(S)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,E)){
					int faceIndex = faces(E)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					if(g->isFaceRegular(i,j,E)){
						coefficients[ind] = (1/V)*h;
					}
					else{
						coefficients[ind] = (1/V)*faces(E)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,W)){
					int faceIndex = faces(W)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					if(g->isFaceRegular(i,j,W)){
						coefficients[ind] = -(1/V)*h;
					}
					else{
						coefficients[ind] = -(1/V)*faces(W)->getArea();
					}
				}
				if(g->isFaceUncovered(i,j,B)){
					int faceIndex = faces(B)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					coefficients[ind] = faces(B)->getArea()/V;
				}
			}
		}
	}
}

}
}
