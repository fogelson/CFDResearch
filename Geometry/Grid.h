/*
 * Grid.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef GRID_H_
#define GRID_H_

#include "Geometry.h"

namespace CFD{
	namespace Geometry{
		CellType Grid :: getType(int i, int j){
			return types(i,j);
		}
		Cell Grid :: getCell(int i, int j){
			return cells(i,j);
		}
		double Grid :: getH(){
			return h;
		}
		double Grid :: getX(int i, int j){
			return cellCenters(i,j)(0);
		}
		double Grid :: getY(int i, int j){
			return cellCenters(i,j)(1);
		}
	}
}

#endif /* GRID_H_ */
