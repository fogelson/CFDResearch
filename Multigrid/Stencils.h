/*
 * Stencils.h
 *
 *  Created on: Aug 19, 2011
 *      Author: bfogelson
 */

#ifndef STENCILS_H_
#define STENCILS_H_

#include "../Geometry/Stencil.h"

namespace CFD{
	using namespace Geometry;
	namespace Multigrid{

	class BackwardEulerDiffusionStencil : public Stencil{
		double deltaT, D;
	public:
		void init(){
			c.resize(g->xRange,g->yRange);
			double h = g->h;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					Coefficients I, cL;
					I = 0;
					I(C) = 1;
					cL = 0;
					if(g->isRegular(i,j)){
						cL(N) = 1/pow2(h);
						cL(S) = 1/pow2(h);
						cL(E) = 1/pow2(h);
						cL(W) = 1/pow2(h);
						cL(C) = -4/pow2(h);
					}
					else if(g->isIrregular(i,j)){
						if(g->faceTypes(i,j)(N) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(S) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(E) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(W) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
					}
					else if(g->isCovered(i,j)){
						cL(C) = 1;
					}
					c(i,j) = I - D*(deltaT)*cL;
				}
			}
		}
	public:
		BackwardEulerDiffusionStencil(Grid * g, double deltaT, double D){
			this->g = g;
			this->deltaT = deltaT;
			this->D = D;
			hasCoarsened = false;
			init();
		}
		Stencil * copy() const{
			return new BackwardEulerDiffusionStencil(g,deltaT,D);
		}
	};

	class PoissonStencil : public Stencil{
	public:
		void init(){
			c.resize(g->xRange,g->yRange);
			double h = g->h;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					Coefficients cL;
					cL = 0;
					if(g->isRegular(i,j)){
						cL(N) = 1/pow2(h);
						cL(S) = 1/pow2(h);
						cL(E) = 1/pow2(h);
						cL(W) = 1/pow2(h);
						cL(C) = -4/pow2(h);
					}
					else if(g->isIrregular(i,j)){
						/*cL += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g);
						cL -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g);
						cL += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g);
						cL -= (g->areaFractions(i,j)(W))* EBUtilities::getGradientCoefficients(i,j,W,g);

						cL = cL/(h*g->volumeFractions(i,j));*/

						// This commented block of code is a very
						// rudimentary way to try Dirichlet boundary
						// conditions, with some added error at the
						// boundaries
						/*cL(N) = 1/pow2(h);
						cL(S) = 1/pow2(h);
						cL(E) = 1/pow2(h);
						cL(W) = 1/pow2(h);
						cL(C) = -4/pow2(h);*/

						if(g->faceTypes(i,j)(N) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(S) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(E) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
						if(g->faceTypes(i,j)(W) != COVERED){
							Coefficients cTemp;
							cTemp = 0;

							cTemp -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g)/(h*(g->volumeFractions(i,j)));
							cL += cTemp;
						}
					}
					else if(g->isCovered(i,j)){
						cL(C) = 1;
					}
					c(i,j) = cL;
				}
			}
		}
	public:
		PoissonStencil(Grid * g){
			this->g = g;
			hasCoarsened = false;
			init();
		}
		Stencil * copy() const{
			return new PoissonStencil(g);
		}
	};

	}
}

#endif /* STENCILS_H_ */
