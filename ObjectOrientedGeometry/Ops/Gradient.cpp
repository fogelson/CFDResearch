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
	/*double c1, c2;
	double h = g->getH();
	c1 = -1/h;
	c2 = 1/h;
	coefficients[ind1] = c1;
	coefficients[ind2] = c2;*/

	/*if(coefficients.count(ind1) > 0 || coefficients.count(ind2) > 0
			|| coefficients.count(ind3) > 0
			|| coefficients.count(ind4) > 0){
		cout << "Overwriting a set coefficient." << endl;
	}*/

	double c1, c2, c3, c4;
	double h = g->getH();
	c1 = -(1 + alpha/h)/(2*h);
	c2 =  (1 + alpha/h)/(2*h);
	c3 = -(1 - alpha/h)/(2*h);
	c4 =  (1 - alpha/h)/(2*h);
	coefficients[ind1] += c1;
	coefficients[ind2] += c2;
	coefficients[ind3] += c3;
	coefficients[ind4] += c4;
}

Gradient::Gradient(Grid * g){
	this->g = g;
	double h = g->getH();
	constantTerm.resize(g->faces.size());
	constantTerm = 0;
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->cells(i,j)->isUncovered()){
				TinyVector<Face*,5> faces = g->cells(i,j)->faces;
				if(g->isFaceUncovered(i,j,B)){
					FaceIndex faceIndex = faces(B)->getIndex();
					Coord centroid = faces(B)->getCentroid();
					TinyVector<double,2> n = faces(B)->getNormal();
					double x, y, r;
					x = centroid(0);
					y = centroid(1);
					r = sqrt(pow2(x) + pow2(y));
					TinyVector<double,2> gradU;
					double pi = acos(-1);
					//gradU(0) = 3*pow(x,2)*pow(y-2,2)+sin(2*(x-pow(y,2)));
					//gradU(1) = 2*pow(x,3)*(y-2) - 2*y*sin(2*(x-pow(y,2)));
					gradU(0) = -2*pi*sin(2*pi*r)*x/r;
					gradU(1) = -2*pi*sin(2*pi*r)*y/r;
					constantTerm(faceIndex) += gradU(0)*n(0) + gradU(1)*n(1);
				}
				if(faces(N) != 0){
					CellIndex cellIndex(i,j);
					FaceIndex faceIndex = faces(N)->getIndex();
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(N)->isRegular()){
						coefficients[c2f] = -1/h;
					}
					else if(faces(N)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(N)->getArea();
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
					FaceIndex faceIndex = faces(S)->getIndex();
					CellIndex cellIndex(i,j);
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(S)->isRegular()){
						coefficients[c2f] = 1/h;
					}
					else if(faces(S)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(S)->getArea();
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
					FaceIndex faceIndex = faces(E)->getIndex();
					CellIndex cellIndex(i,j);
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(E)->isRegular()){
						coefficients[c2f] = -1/h;
					}
					else if(faces(E)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(E)->getArea();
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
					FaceIndex faceIndex = faces(W)->getIndex();
					CellIndex cellIndex(i,j);
					CellToFaceIndex c2f(i,j,faceIndex);
					if(faces(W)->isRegular()){
						coefficients[c2f] = 1/h;
					}
					else if(faces(W)->isIrregular()){
						if(coefficients.count(c2f) == 0){
							double alpha = g->cells(i,j)->getFace(W)->getArea();
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