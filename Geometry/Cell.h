/*
 * Cell.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef CELL_H_
#define CELL_H_

#include "Geometry.h"

namespace CFD{
	namespace Geometry{
		Cell :: Cell(){

		}
		Cell :: Cell(Grid *grid, int i, int j){
			setGrid(grid);
			setIndices(i,j);
		}
		void Cell :: setGrid(Grid *grid){
			this->grid = grid;
		}
		Grid* Cell :: getGrid(){
			return grid;
		}
		int Cell :: getI(){
			return i;
		}
		int Cell :: getJ(){
			return j;
		}
		void Cell :: setIndices(int i, int j){
			this->i = i;
			this->j = j;
		}
		Coord Cell :: getCenter(){
			Coord c;
			c(0) = grid->getX(i,j);
			c(1) = grid->getY(i,j);
			return c;
		}
	}
}

#endif /* CELL_H_ */
