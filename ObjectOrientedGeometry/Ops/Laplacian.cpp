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

Laplacian::Laplacian(Grid * g){
	this->g = g;
	Gradient grad(g);
	Divergence div(g);
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = 0;
	CellDoubleArray divGradConstant = g->makeCellDoubleArray();
	divGradConstant = div(grad.constantTerm);
	constantTerm = divGradConstant + div.constantTerm;

	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			for(int k = 0; k < g->faces.size(); k++){
				CellToFaceIndex c2f(i,j,k);
			}
		}
	}
}

/*void Laplacian::applyIrregularGradient(double & c1, double & c2, double & c3, double & c4, double a, double v){
	double h = g->getH();
	c1 = -(a/h + pow2(a/h))/(2*v);
	c2 = (a/h + pow2(a/h))/(2*v);
	c3 = -(a/h - pow2(a/h))/(2*v);
	c4 = (a/h - pow2(a/h))/(2*v);
}

double Laplacian::getDirichlet(Coord c){
	return getDirichlet(c(0),c(1));
}
double Laplacian::getDirichlet(double x, double y){
	double r = sqrt(pow2(x) + pow2(y));
	double pi = acos(-1);
	return cos(2*pi*r);
}

Laplacian::Laplacian(Grid * g){
	this->g = g;
	constantTerm.resize(g->xRange,g->yRange);
	constantTerm = 0;

	double h = g->getH();
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isRegular(i,j)){
				double cN, cS, cE, cW, cC;
				cN = 1/pow2(h);
				cS = 1/pow2(h);
				cE = 1/pow2(h);
				cW = 1/pow2(h);
				cC = -4/pow2(h);

				CellToCellIndex indN(i,j+1,i,j);
				CellToCellIndex indS(i,j-1,i,j);
				CellToCellIndex indE(i+1,j,i,j);
				CellToCellIndex indW(i-1,j,i,j);
				CellToCellIndex indC(i,j,i,j);

				coefficients[indN] = cN;
				coefficients[indS] = cS;
				coefficients[indE] = cE;
				coefficients[indW] = cW;
				coefficients[indC] = cC;
			}
			else if(g->isIrregular(i,j)){
				double v = g->cells(i,j)->getVolume();
				if(false && g->cells(i,j)->getFace(B) != 0){
					double pi = acos(-1);
					Coord centroid = g->cells(i,j)->getFace(B)->getCentroid();
					double x = centroid(0), y = centroid(1);
					double r = sqrt(pow2(x) + pow2(y));
					TinyVector<double,2> n = g->cells(i,j)->getFace(B)->getNormal();
					double gradUx, gradUy;
					gradUx = -2*pi*sin(2*pi*r)*(x/r);
					gradUy = -2*pi*sin(2*pi*r)*(y/r);
					double gradUdotN = gradUx*n(0) + gradUy*n(1);
					double flux = gradUdotN*g->cells(i,j)->getFace(B)->getArea();
					constantTerm(i,j) += flux/v;
				}
				if(g->isFaceUncovered(i,j,N)){
					if(g->isFaceRegular(i+1,j,N)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j,i,j);
						CellToCellIndex ind2(i,j+1,i,j);
						CellToCellIndex ind3(i+1,j,i,j);
						CellToCellIndex ind4(i+1,j+1,i,j);
						double a = g->cells(i,j)->getFace(N)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += c1;
						coefficients[ind2] += c2;
						coefficients[ind3] += c3;
						coefficients[ind4] += c4;
					}
					else if(g->isFaceRegular(i-1,j,N)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j,i,j);
						CellToCellIndex ind2(i,j+1,i,j);
						CellToCellIndex ind3(i-1,j,i,j);
						CellToCellIndex ind4(i-1,j+1,i,j);
						double a = g->cells(i,j)->getFace(N)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += c1;
						coefficients[ind2] += c2;
						coefficients[ind3] += c3;
						coefficients[ind4] += c4;
					}
				}
				if(g->isFaceUncovered(i,j,S)){
					if(g->isFaceRegular(i+1,j,S)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j-1,i,j);
						CellToCellIndex ind2(i,j,i,j);
						CellToCellIndex ind3(i+1,j-1,i,j);
						CellToCellIndex ind4(i+1,j,i,j);
						double a = g->cells(i,j)->getFace(S)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += -c1;
						coefficients[ind2] += -c2;
						coefficients[ind3] += -c3;
						coefficients[ind4] += -c4;
					}
					else if(g->isFaceRegular(i-1,j,S)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j-1,i,j);
						CellToCellIndex ind2(i,j,i,j);
						CellToCellIndex ind3(i-1,j-1,i,j);
						CellToCellIndex ind4(i-1,j,i,j);
						double a = g->cells(i,j)->getFace(S)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += -c1;
						coefficients[ind2] += -c2;
						coefficients[ind3] += -c3;
						coefficients[ind4] += -c4;
					}
				}
				if(g->isFaceUncovered(i,j,E)){
					if(g->isFaceRegular(i,j+1,E)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j,i,j);
						CellToCellIndex ind2(i+1,j,i,j);
						CellToCellIndex ind3(i,j+1,i,j);
						CellToCellIndex ind4(i+1,j+1,i,j);
						double a = g->cells(i,j)->getFace(E)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += c1;
						coefficients[ind2] += c2;
						coefficients[ind3] += c3;
						coefficients[ind4] += c4;
					}
					else if(g->isFaceRegular(i,j-1,E)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i,j,i,j);
						CellToCellIndex ind2(i+1,j,i,j);
						CellToCellIndex ind3(i,j-1,i,j);
						CellToCellIndex ind4(i+1,j-1,i,j);
						double a = g->cells(i,j)->getFace(E)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += c1;
						coefficients[ind2] += c2;
						coefficients[ind3] += c3;
						coefficients[ind4] += c4;
					}
				}
				if(g->isFaceUncovered(i,j,W)){
					if(g->isFaceRegular(i,j+1,W)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i-1,j,i,j);
						CellToCellIndex ind2(i,j,i,j);
						CellToCellIndex ind3(i-1,j+1,i,j);
						CellToCellIndex ind4(i,j+1,i,j);
						double a = g->cells(i,j)->getFace(W)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += -c1;
						coefficients[ind2] += -c2;
						coefficients[ind3] += -c3;
						coefficients[ind4] += -c4;
					}
					else if(g->isFaceRegular(i,j-1,W)){
						double c1, c2, c3, c4;
						CellToCellIndex ind1(i-1,j,i,j);
						CellToCellIndex ind2(i,j,i,j);
						CellToCellIndex ind3(i-1,j-1,i,j);
						CellToCellIndex ind4(i,j-1,i,j);
						double a = g->cells(i,j)->getFace(W)->getArea();

						applyIrregularGradient(c1,c2,c3,c4,a,v);

						coefficients[ind1] += -c1;
						coefficients[ind2] += -c2;
						coefficients[ind3] += -c3;
						coefficients[ind4] += -c4;
					}
				}
			}
		}
	}
}*/

}
}
