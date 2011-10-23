/*
 * VolumeWeightedRestrictor.cpp
 *
 *  Created on: Oct 21, 2011
 *      Author: fogelson
 */


#include "TransferOperators.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
namespace OOMultigrid{


void VolumeWeightedRestrictor::doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * fine, Grid * coarse){
	uC = 0;
	for(int iC = coarse->iMin; iC <= coarse->iMax; iC++){
		for(int jC = coarse->jMin; jC <= coarse->jMax; jC++){
			if(coarse->isUncovered(iC,jC)){
				double vNE = 0, vNW = 0, vSE = 0, vSW = 0;
				double uNE = 0, uNW = 0, uSE = 0, uSW = 0;

				if(fine->isUncovered(2*iC,2*jC)){
					uNE = uF(2*iC,2*jC);
					vNE = fine->cells(2*iC,2*jC)->getVolume();
				}
				if(fine->isUncovered(2*iC-1,2*jC)){
					uNW = uF(2*iC-1,2*jC);
					vNW = fine->cells(2*iC-1,2*jC)->getVolume();
				}
				if(fine->isUncovered(2*iC,2*jC-1)){
					uSE = uF(2*iC,2*jC-1);
					vSE = fine->cells(2*iC,2*jC-1)->getVolume();
				}
				if(fine->isUncovered(2*iC-1,2*jC-1)){
					uSW = uF(2*iC-1,2*jC-1);
					vSW = fine->cells(2*iC-1,2*jC-1)->getVolume();
				}

				uC(iC,jC) = (vNE*uNE + vNW*uNW + vSE*uSE + vSW*uSW)/(vNE + vNW + vSE + vSW);
			}
		}
	}
}


}
}
