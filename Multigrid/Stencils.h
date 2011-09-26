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
		Coefficients gradient(Direction aP, Direction aM, Direction bP, Direction bM, double alpha, double h){
			Coefficients c;
			c = 0;
			c(bP) = (1-alpha)/2;
			c(bM) = -(1-alpha)/2;
			c(aP) = (1+alpha)/2;
			c(aM) = -(1+alpha)/2;

			return c;
		}
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

						//cL = 0;
						//cL(C) = 1;

						Coefficients fN, fS, fE, fW;
						fN = 0;
						fS = 0;
						fE = 0;
						fW = 0;

						if(g->isFaceRegular(i,j,N)){
							fN(N) = 1;
							fN(C) = -1;
						}
						else if(g->isFaceIrregular(i,j,N)){
							Direction aM, aP, bM, bP;

							double alpha = g->areaFractions(i,j)(N);

							aP = N;
							aM = C;

							if(g->isFaceRegular(i+1,j,N)){
								bP = NE;
								bM = E;
							}
							else if(g->isFaceRegular(i-1,j,N)){
								bP = NW;
								bM = W;
							}
							else{
								cout << "Error near boundary" << endl;
								bP = C;
								bM = C;
							}
							fN = alpha*gradient(aP,aM,bP,bM,alpha,h);
						}
						if(g->isFaceRegular(i,j,E)){
							fE(E) = 1;
							fE(C) = -1;
						}
						else if(g->isFaceIrregular(i,j,E)){
							Direction aM, aP, bM, bP;
							double alpha = g->areaFractions(i,j)(E);

							aP = E;
							aM = C;

							if(g->isFaceRegular(i,j+1,E)){
								bP = NE;
								bM = N;
							}
							else if(g->isFaceRegular(i,j-1,E)){
								bP = SE;
								bM = S;
							}
							else{
								cout << "Error near boundary" << endl;
								bP = C;
								bM = C;
							}
							fE = alpha*gradient(aP,aM,bP,bM,alpha,h);
						}
						if(g->isFaceRegular(i,j,S)){
							fS(C) = 1;
							fS(S) = -1;
						}
						else if(g->isFaceIrregular(i,j,S)){
							Direction aM, aP, bM, bP;
							double alpha = g->areaFractions(i,j)(S);

							aP = C;
							aM = S;

							if(g->isFaceRegular(i+1,j,S)){
								bP = E;
								bM = SE;
							}
							else if(g->isFaceRegular(i-1,j,S)){
								bP = W;
								bM = SW;
							}
							else{
								cout << "Error near boundary" << endl;
								bP = C;
								bM = C;
							}
							fS = alpha*gradient(aP,aM,bP,bM,alpha,h);
						}
						if(g->isFaceRegular(i,j,W)){
							fW(C) = 1;
							fW(W) = -1;
						}
						else if(g->isFaceIrregular(i,j,W)){
							Direction aM, aP, bM, bP;
							double alpha = g->areaFractions(i,j)(W);

							aP = C;
							aM = W;

							if(g->isFaceRegular(i,j+1,W)){
								bP = N;
								bM = NW;
							}
							else if(g->isFaceRegular(i,j-1,W)){
								bP = S;
								bM = SW;
							}
							else{
								cout << "Error near boundary" << endl;
							}
							fW = alpha*gradient(aP,aM,bP,bM,alpha,h);
						}
						Coefficients f;
						//f = (fN - fS + fE - fW)/(h*h);
						f = (fN - fS + fE - fW)/((g->volumeFractions(i,j))*(h*h));
						cL = f;


						/*if(g->faceTypes(i,j)(N) != COVERED){
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
						}*/
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
