/*
 * PiecewiseConstantInterpolator.cpp
 *
 *  Created on: Oct 25, 2011
 *      Author: fogelson
 */


#include "TransferOperators.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{


void PiecewiseConstantInterpolator::doInterpolate(CellDoubleArray & uC, CellDoubleArray & uF, Grid * coarse, Grid * fine){
	uF = 0;
	for(int iC = coarse->iMin; iC <= coarse->iMax; iC++){
		for(int jC = coarse->jMin; jC <= coarse->jMax; jC++){
			if(coarse->isUncovered(iC,jC)){
				uF(2*iC,2*jC) = uC(iC,jC);
				uF(2*iC-1,2*jC) = uC(iC,jC);
				uF(2*iC,2*jC-1) = uC(iC,jC);
				uF(2*iC-1,2*jC-1) = uC(iC,jC);
			}
		}
	}
}

}
}
