/*
 * Stencil.h
 *
 *  Created on: Jun 13, 2011
 *      Author: fogelson
 */

#ifndef STENCIL_H_
#define STENCIL_H_

#include "Geometry.h"

namespace CFD{
	namespace Geometry{
		class Stencil{
		private:
			double apply(Array<Coefficients,2> c, CellDoubleArray &u, int i, int j){
				double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;

				Type type = g->getCellType(i,j);
				if(type != COVERED){
					if(type == REGULAR){
						uC = u(i,j);
						uN = u(i,j+1);
						uS = u(i,j-1);
						uE = u(i+1,j);
						uW = u(i-1,j);
						uNE = u(i+1,j+1);
						uNW = u(i-1,j+1);
						uSE = u(i+1,j-1);
						uSW = u(i-1,j-1);
					}
					else{
						uC = g->isUncovered(i,j) ? u(i,j) : 0;
						uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
						uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
						uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
						uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
						uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
						uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
						uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
						uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;
					}

					double Su = c(i,j)(C)*uC
							+ c(i,j)(N)*uN
							+ c(i,j)(S)*uS
							+ c(i,j)(E)*uE
							+ c(i,j)(W)*uW
							+ c(i,j)(NW)*uNW
							+ c(i,j)(NE)*uNE
							+ c(i,j)(SW)*uSW
							+ c(i,j)(SE)*uSE;
					return Su;
				}
				else{
					return 0;
				}
			}
		protected:
			Array<Coefficients,2> c;
			Grid * g;
			Stencil * sCoarse;
			bool hasCoarsened;
			virtual void init() = 0;
		public:
			virtual ~Stencil(){
				if(hasCoarsened){
					delete sCoarse;
				}
			}
			virtual Stencil * copy() const = 0;
			void setGrid(Grid * g){
				this->g = g;
			}
			Grid * getGrid(){
				return g;
			}
			Stencil * coarsen(){
				if(hasCoarsened){
					return sCoarse;
				}
				else{
					Grid * gCoarse = g->coarsen();
					sCoarse = this->copy();
					sCoarse->setGrid(gCoarse);
					sCoarse->init();
					hasCoarsened = true;
					return sCoarse;
				}
			}
			virtual double apply(CellDoubleArray &u, int i, int j){
				return apply(c,u,i,j);
			}
			virtual CellDoubleArray apply(CellDoubleArray &u){
				CellDoubleArray Su = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						Su(i,j) = apply(u, i, j);
					}
				}
				return Su;
			}
			CellDoubleArray operator () (CellDoubleArray &u){
				return apply(u);
			}
			Coefficients getCoefficients(int i, int j){
				return c(i,j);
			}
			Array<Coefficients,2> getCoefficients(){
				return c;
			}
		};
		class IdentityStencil : public Stencil{
			void init(){
				c.resize(g->xRange,g->yRange);
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isCovered(i,j)){
							c(i,j) = 0;
						}
						else{
							c(i,j) = 0;
							c(i,j)(C) = 1;
						}
					}
				}
			}
		public:
			IdentityStencil(Grid * g){
				this->g = g;
				hasCoarsened = false;
				init();
			}
			Stencil * copy() const{
				return new IdentityStencil(g);
			}
		};
	}
}

#endif /* STENCIL_H_ */
