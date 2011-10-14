/*
 * Gradient.cpp
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

void Gradient::interpolateIrregularFace(double alpha, CellToFaceIndex ind1, CellToFaceIndex ind2, CellToFaceIndex ind3, CellToFaceIndex ind4){
	double c1, c2;
	double h = g->getH();
	c1 = -1/h;
	c2 = 1/h;
	coefficients[ind1] = c1;
	coefficients[ind2] = c2;

	/*double c1, c2, c3, c4;
	double h = g->getH();
	c1 = -(h - alpha)/(pow2(h));
	c2 =  (h - alpha)/(pow2(h));
	c3 = -alpha/pow2(h);
	c4 =  alpha/pow2(h);
	coefficients[ind1] = c1;
	coefficients[ind2] = c2;
	coefficients[ind3] = c3;
	coefficients[ind4] = c4;*/
}

Gradient::Gradient(Grid * g){
	this->g = g;
	double h = g->getH();
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->cells(i,j)->isUncovered()){
				TinyVector<Face*,5> faces = g->cells(i,j)->faces;
				if(faces(N) != 0){
					int faceIndex = faces(N)->getIndex();
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(N)->isRegular()){
						coefficients[c2f] = -1/h;
					}
					else if(faces(N)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(N)->getArea()/2;
							if(g->isFaceUncovered(i+1,j,N)){
								CellToFaceIndex ind1(i,j,faceIndex);
								CellToFaceIndex ind2(i,j+1,faceIndex);
								CellToFaceIndex ind3(i+1,j,faceIndex);
								CellToFaceIndex ind4(i+1,j+1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
							else if(g->isFaceUncovered(i-1,j,N)){
								CellToFaceIndex ind1(i,j,faceIndex);
								CellToFaceIndex ind2(i,j+1,faceIndex);
								CellToFaceIndex ind3(i-1,j,faceIndex);
								CellToFaceIndex ind4(i-1,j+1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
						}
					}
				}
				if(faces(S) != 0){
					int faceIndex = faces(S)->getIndex();
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(S)->isRegular()){
						coefficients[c2f] = 1/h;
					}
					else if(faces(S)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(S)->getArea()/2;
							if(g->isFaceUncovered(i+1,j,S)){
								CellToFaceIndex ind1(i,j-1,faceIndex);
								CellToFaceIndex ind2(i,j,faceIndex);
								CellToFaceIndex ind3(i+1,j-1,faceIndex);
								CellToFaceIndex ind4(i+1,j,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
							else if(g->isFaceUncovered(i-1,j,S)){
								CellToFaceIndex ind1(i,j-1,faceIndex);
								CellToFaceIndex ind2(i,j,faceIndex);
								CellToFaceIndex ind3(i-1,j-1,faceIndex);
								CellToFaceIndex ind4(i-1,j,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
						}
					}
				}
				if(faces(E) != 0){
					int faceIndex = faces(E)->getIndex();
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(E)->isRegular()){
						coefficients[c2f] = -1/h;
					}
					else if(faces(E)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(E)->getArea()/2;
							if(g->isFaceUncovered(i,j+1,E)){
								CellToFaceIndex ind1(i,j,faceIndex);
								CellToFaceIndex ind2(i+1,j,faceIndex);
								CellToFaceIndex ind3(i,j+1,faceIndex);
								CellToFaceIndex ind4(i+1,j+1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
							else if(g->isFaceUncovered(i,j-1,E)){
								CellToFaceIndex ind1(i,j,faceIndex);
								CellToFaceIndex ind2(i+1,j,faceIndex);
								CellToFaceIndex ind3(i,j-1,faceIndex);
								CellToFaceIndex ind4(i+1,j-1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
						}
					}
				}
				if(faces(W) != 0){
					int faceIndex = faces(W)->getIndex();
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(W)->isRegular()){
						coefficients[c2f] = 1/h;
					}
					else if(faces(W)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(W)->getArea()/2;
							if(g->isFaceUncovered(i,j+1,W)){
								CellToFaceIndex ind1(i+1,j,faceIndex);
								CellToFaceIndex ind2(i,j,faceIndex);
								CellToFaceIndex ind3(i+1,j+1,faceIndex);
								CellToFaceIndex ind4(i,j+1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
							else if(g->isFaceUncovered(i,j-1,W)){
								CellToFaceIndex ind1(i+1,j,faceIndex);
								CellToFaceIndex ind2(i,j,faceIndex);
								CellToFaceIndex ind3(i+1,j-1,faceIndex);
								CellToFaceIndex ind4(i,j-1,faceIndex);
								interpolateIrregularFace(alpha,ind1,ind2,ind3,ind4);
							}
						}
					}
				}
			}
		}
	}

	//cout << "Resized" << endl;
/*
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			cout << "i = " << i << " j = " << j << endl;
			if(g->cells(i,j)->isUncovered()){
				cout << "Uncovered" << endl;
				double c = 1;
				CellIndex index(i,j);
				TinyVector<Face*,5> faces = g->cells(i,j)->faces;
				cout << "Got faces" << endl;
				for(int k = 0; k < 5; k++){
					cout << "Testing face " << k << endl;
					if(faces(k) != 0 && faces(k)->isUncovered()){
						cout << "Face " << k << " is uncovered ";
						cout << "and has index " << faces(k)->getIndex() << endl;
						cout << coefficients(faces(k)->getIndex()).empty() << endl;
						coefficients(faces(k)->getIndex())[index] = c;
					}
				}
			}
		}
	}*/
}

}
}
