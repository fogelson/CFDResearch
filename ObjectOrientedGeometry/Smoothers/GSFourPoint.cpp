/*
 * GSFourPoint.cpp
 *
 *  Created on: Oct 26, 2011
 *      Author: fogelson
 */


#include "Smoothers.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{

void GSFourPoint::smoothPoint(int & i, int & j, CellDoubleArray & u, CellDoubleArray & u0, CellDoubleArray & f, CellToCellOperator & C){
	double rhs = f(i,j) - C.constantTerm(i,j);
	CellCoefficients::iterator cIt;
	CellIndex current(i,j);
	double currentCoefficient = 0;
	for(cIt = C.coefficients[current].begin(); cIt != C.coefficients[current].end(); cIt++){
		if((*cIt).first != current){
			int iRhs = (*cIt).first.i, jRhs = (*cIt).first.j;
			rhs -= (*cIt).second*u(iRhs,jRhs);
		}
		else{
			currentCoefficient = (*cIt).second;
		}
	}
	u(i,j) = rhs/currentCoefficient;
}

void GSFourPoint::smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, CellToCellOperator & C, int its){
	Grid * g = C.g;
	u = u0;

	for(int n = 0; n < its; n++){
		for(int i = g->iMin; i <= g->iMax; i++){
			for(int j = g->jMin; j <= g->jMax; j++){
				this->smoothPoint(i,j,u,u0,f,C);
			}
		}
		for(int i = g->iMax; i >= g->iMin; i--){
			for(int j = g->jMin; j <= g->jMax; j++){
				this->smoothPoint(i,j,u,u0,f,C);
			}
		}
		for(int i = g->iMin; i <= g->iMax; i++){
			for(int j = g->jMax; j >= g->jMin; j--){
				this->smoothPoint(i,j,u,u0,f,C);
			}
		}
		for(int i = g->iMax; i >= g->iMin; i--){
			for(int j = g->jMax; j >= g->jMin; j--){
				this->smoothPoint(i,j,u,u0,f,C);
			}
		}
	}
}

}
}
