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
					coefficients[ind] = faces(N)->getArea()/V;
				}
				if(g->isFaceUncovered(i,j,S)){
					int faceIndex = faces(S)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					coefficients[ind] = -faces(S)->getArea()/V;
				}
				if(g->isFaceUncovered(i,j,E)){
					int faceIndex = faces(E)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					coefficients[ind] = faces(E)->getArea()/V;
				}
				if(g->isFaceUncovered(i,j,W)){
					int faceIndex = faces(W)->getIndex();
					CellToFaceIndex ind(i,j,faceIndex);
					coefficients[ind] = -faces(W)->getArea()/V;
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
