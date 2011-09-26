/*
 * Stencils.h
 *
 *  Created on: Aug 18, 2011
 *      Author: fogelson
 */

#ifndef MULTIGRIDSTENCILS_H_
#define MULTIGRIDSTENCILS_H_

#include "../Geometry/Stencil.h"

namespace CFD{
	using namespace Geometry;
	namespace Multigrid{
		class PoissonStencil : public Stencil{
		public:
			PoissonStencil(Grid * g){
				this->g = g;
				hasCoarsened = false;
				init();
			}
			void init(){
				c.resize(g->xRange,g->yRange);
				double h = g->h;
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						Coefficients cD;
						cD = 0;
						if(g->isIrregular(i,j)){
							cD += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g);
							cD += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g);
							cD -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g);
							cD -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g);

							cD = (1.0/(h*g->volumeFractions(i,j)))*cD;
						}
						else if(g->isRegular(i,j)){
							cD(C) = -4.0/pow2(h);
							cD(N) = 1.0/pow2(h);
							cD(S) = 1.0/pow2(h);
							cD(E) = 1.0/pow2(h);
							cD(W) = 1.0/pow2(h);
						}
						else{
							cD = 0;
							cD(C) = 1;
						}
						c(i,j) = cD;
					}
				}
			}
			Stencil * copy() const{
				return new PoissonStencil(g);
			}
		};
	}
}

#endif /* MULTIGRIDSTENCILS_H_ */
