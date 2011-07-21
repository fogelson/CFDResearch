/*
 * Square.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef SQUARE_H_
#define SQUARE_H_

#include "Grid.h"

namespace CFD{
	namespace Geometry{
		class Square : public Grid{
			double s; // Side length
			Coord offset;
			Range xRange, yRange;
			void make(){
				int M = ceil(s/h);
				xRange.setRange(0,M+1,1);
				yRange.setRange(0,M+1,1);
				cellCenters.resize(xRange,yRange);
				cells.resize(xRange,yRange);
				types.resize(xRange,yRange);
				for(int i = 0; i <= M+1; i++){
					for(int j = 0; j <= M+1; j++){
						cellCenters(i,j)(0) = offset(0) + (i-1)*h + 0.5*h;
						cellCenters(i,j)(1) = offset(1) + (j-1)*h + 0.5*h;
						types = IRREGULARCELL;
						//types(i,j) = discoverType(cells(i,j));
					}
				}
			}
			CellType discoverType(Cell &c){
				CellType type;
				if(contains(c.getCenter())){
					type = REGULARCELL;
				}
				return type;
			}
		public:
			Square(double s, double h){
				this->h = h;
				this->s = s;
				offset(0) = 0;
				offset(1) = 0;
				make();
			}
			Square(double s, double h, Coord offset){
				this->h = h;
				this->s = s;
				this->offset = offset;
			}
			Range getXRange(){
				return xRange;
			}
			Range getYRange(){
				return yRange;
			}
			bool contains(double x, double y){
				return true;
			}
			bool contains(Coord c){
				double x = c(0), y = c(1);
				return contains(x,y);
			}
		};
	}
}

#endif /* SQUARE_H_ */
