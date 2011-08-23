/*
 * Smoothers.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef SMOOTHERS_H_
#define SMOOTHERS_H_

#include <blitz/tinyvec-et.h>
#include "../Geometry/Geometry.h"
#include "../Geometry/Stencil.h"
//#include "../LinearAlgebra/Krylov.h"

namespace CFD{
	using namespace Geometry;
	namespace Multigrid{
		class Smoother;

		/* Abstract base class for iterative solvers the smooth the solution
		 * quickly.
		 */
		class Smoother{
		public:
			virtual ~Smoother(){}
			virtual CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its) = 0;
		};

		class StenciledSmoother{
		public:
			virtual ~StenciledSmoother(){}
			virtual void smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, Stencil * s, int its) = 0;
			virtual CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Stencil * s, int its){
				Grid * g = s->getGrid();
				CellDoubleArray u = g->makeCellDoubleArray();
				smooth(u,u0,f,s,its);
				return u;
			}
		};

		class StenciledGSLex : public StenciledSmoother{
		public:
			void smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, Stencil * s, int its){
				Grid * g = s->getGrid();
				u = u0;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(g->isUncovered(i,j)){
								double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;

								if(g->isRegular(i,j)){
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
								Coefficients c;
								c = s->getCoefficients(i,j);

								u(i,j) = (f(i,j) - (c(N)*uN + c(S)*uS + c(E)*uE + c(W)*uW
										+ c(NE)*uNE + c(NW)*uNW + c(SE)*uSE + c(SW)*uSW))/c(C);
							}
						}
					}
				}
			}
		};

		class StenciledFourPointGS : public StenciledSmoother{
			void smoothPoint(int i, int j, CellDoubleArray &u, CellDoubleArray &f, Stencil * s){
				Grid * g = s->getGrid();
				Type type = g->getCellType(i,j);
				if(type != COVERED){
					double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
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
					Coefficients c;
					c = s->getCoefficients(i,j);

					u(i,j) = (f(i,j) - (c(N)*uN + c(S)*uS + c(E)*uE + c(W)*uW
							+ c(NE)*uNE + c(NW)*uNW + c(SE)*uSE + c(SW)*uSW))/c(C);
				}
			}
		public:
			void smooth(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray f, Stencil * s, int its){
				Grid * g = s->getGrid();
				u = u0;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							smoothPoint(i,j,u,f,s);
						}
					}
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMax; j >= g->jMin; j--){
							smoothPoint(i,j,u,f,s);
						}
					}
					for(int i = g->iMax; i >= g->iMin; i--){
						for(int j = g->jMin; j <= g->jMax; j++){
							smoothPoint(i,j,u,f,s);
						}
					}
					for(int i = g->iMax; i >= g->iMin; i--){
						for(int j = g->jMax; j >= g->jMin; j--){
							smoothPoint(i,j,u,f,s);
						}
					}
				}
			}
		};

		class PoissonGSRB : public Smoother{
		public:
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid * g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				u = u0;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(abs(i + j) % 2 == 0){ // Red cells
								if(g->isRegular(i,j)){
									u(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - pow2(h)*f(i,j))/4;
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*h*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*h*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*h*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*h*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									u(i,j) = (pow2(h)*(g->volumeFractions(i,j))*f(i,j)
											- coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
											- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
											- coefficients(E)*(g->isUncovered(i+1,j) ? u(i+1,j) : 0)
											- coefficients(W)*(g->isUncovered(i-1,j) ? u(i-1,j) : 0)
											- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
											- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
											- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
											- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0))/coefficients(C);
								}
							}
						}
					} // End of black cells
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(abs(i + j) % 2 == 1){ // Black cells
								if(g->isRegular(i,j)){
									u(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - pow2(h)*f(i,j))/4;
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*h*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*h*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*h*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*h*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									u(i,j) = (pow2(h)*(g->volumeFractions(i,j))*f(i,j)
											- coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
											- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
											- coefficients(E)*(g->isUncovered(i+1,j) ? u(i+1,j) : 0)
											- coefficients(W)*(g->isUncovered(i-1,j) ? u(i-1,j) : 0)
											- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
											- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
											- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
											- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0))/coefficients(C);
								}
							}
						}
					} // End of black cells
				}
				return u;
			}
		};

		class PoissonGSLex : public Smoother{
		public:
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				u = u0;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(g->isRegular(i,j)){
								u(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - pow2(h)*f(i,j))/4;
							}
							else if(g->isIrregular(i,j)){
								TinyVector<double,9> coefficients;
								coefficients = 0;

								coefficients += (g->areaFractions(i,j)(N))*h*EBUtilities::getGradientCoefficients(u,i,j,N,g);
								coefficients += (g->areaFractions(i,j)(E))*h*EBUtilities::getGradientCoefficients(u,i,j,E,g);
								coefficients -= (g->areaFractions(i,j)(S))*h*EBUtilities::getGradientCoefficients(u,i,j,S,g);
								coefficients -= (g->areaFractions(i,j)(W))*h*EBUtilities::getGradientCoefficients(u,i,j,W,g);

								u(i,j) = (pow2(h)*(g->volumeFractions(i,j))*f(i,j)
										- coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
										- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
										- coefficients(E)*(g->isUncovered(i+1,j) ? u(i+1,j) : 0)
										- coefficients(W)*(g->isUncovered(i-1,j) ? u(i-1,j) : 0)
										- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
										- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
										- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
										- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0))/coefficients(C);
							}
						}
					}
				}
				return u;
			}
		};
/*
		// Point lexicographic Gauss Seidel smoother for the Poisson problem
		class PoissonGSLex : public Smoother{
		public:
			virtual ~PoissonGSLex(){}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its){
				double h = circ.getH();
				GridScalar u = circ.makeScalar();
				u = u0;
				for(int n = 0; n < its; n++){
					for(int i = u.lbound(0); i <= u.ubound(0); i++){
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(circ.getType(i,j) != EXTERIOR){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
				}
				return u;
			}
		};

		// Symmetric point lexicographic Gauss Seidel smoother for the Poisson problem
		class UpwindSymmetricGSLex : public Smoother{
			double aX, aY;
		public:
			virtual ~UpwindSymmetricGSLex(){}
			UpwindSymmetricGSLex(){}
			UpwindSymmetricGSLex(double aX, double aY){
				this->aX = aX;
				this->aY = aY;
			}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its){
				double h = circ.getH();
				GridScalar u = circ.makeScalar();
				u = u0;
				for(int n = 0; n < its; n++){
					// Forward sweep
					for(int i = u.lbound(0); i <= u.ubound(0); i++){
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(circ.getType(i,j) != EXTERIOR){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = f(i,j);
									double denominator = 0;
									if(aX > 0){
										u(i,j) += (aX/h)*u(i-1,j);
										denominator += (aX/h);
									}
									else{
										u(i,j) += -(aX/h)*u(i+1,j);
										denominator += -(aX/h);
									}
									if(aY > 0){
										u(i,j) += (aY/h)*u(i,j-1);
										denominator += (aY/h);
									}
									else{
										u(i,j) += -(aY/h)*u(i,j+1);
										denominator += -(aY/h);
									}
									u(i,j) = u(i,j)/denominator;
								}
								else if(circ.getType(i,j) == IRREGULAR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = f(i,j);
									double denominator = 0;
									if(aX > 0){
										u(i,j) += (aX/hW)*uW;
										denominator += (aX/hW);
									}
									else{
										u(i,j) += -(aX/hE)*uE;
										denominator += -(aX/hE);
									}
									if(aY > 0){
										u(i,j) += (aY/hS)*uS;
										denominator += (aY/hS);
									}
									else{
										u(i,j) += -(aY/hN)*uN;
										denominator += -(aY/hN);
									}
									u(i,j) = u(i,j)/denominator;
								}
							}
						}
					}
					// Backward sweep
					for(int i = u.ubound(0); i >= u.lbound(0); i--){
						for(int j = u.ubound(1); j >= u.lbound(1); j--){
							if(circ.getType(i,j) != EXTERIOR){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
				}
				return u;
			}
		};

		// Symmetric point lexicographic Gauss Seidel smoother for the Poisson problem
		class PoissonSymmetricGSLex : public Smoother{
		public:
			virtual ~PoissonSymmetricGSLex(){}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its){
				double h = circ.getH();
				GridScalar u = circ.makeScalar();
				u = u0;
				for(int n = 0; n < its; n++){
					// Forward sweep
					for(int i = u.lbound(0); i <= u.ubound(0); i++){
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(circ.getType(i,j) != EXTERIOR){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
					// Backward sweep
					for(int i = u.ubound(0); i >= u.lbound(0); i--){
						for(int j = u.ubound(1); j >= u.lbound(1); j--){
							if(circ.getType(i,j) != EXTERIOR){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
				}
				return u;
			}
		};

		// Point Red Black Gauss Seidel smoother for the Poisson problem
		class PoissonGSRB : public Smoother{
		public:
			virtual ~PoissonGSRB(){}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its){
				double h = circ.getH();
				GridScalar u = circ.makeScalar();
				u = u0;
				for(int n = 0; n < its; n++){
					for(int i = u.lbound(0); i <= u.ubound(0); i++){
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(circ.getType(i,j) != EXTERIOR && ((abs(i) + abs(j)) % 2 == 1)){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									double uN, uS, uE, uW;
									double hN, hS, hE, hW;
									Spacing sp = circ.getSpacing(i,j);
									hN = sp(NORTH);
									hS = sp(SOUTH);
									hE = sp(EAST);
									hW = sp(WEST);
									if(circ.getType(i,j+1) == EXTERIOR){
										uN = circ.boundaryValue(i,j,NORTH);
									}
									else{
										uN = u(i,j+1);
									}
									if(circ.getType(i,j-1) == EXTERIOR){
										uS = circ.boundaryValue(i,j,SOUTH);
									}
									else{
										uS = u(i,j-1);
									}
									if(circ.getType(i+1,j) == EXTERIOR){
										uE = circ.boundaryValue(i,j,EAST);
									}
									else{
										uE = u(i+1,j);
									}
									if(circ.getType(i-1,j) == EXTERIOR){
										uW = circ.boundaryValue(i,j,WEST);
									}
									else{
										uW = u(i-1,j);
									}
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
					for(int i = u.lbound(0); i <= u.ubound(0); i++){
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(circ.getType(i,j) != EXTERIOR && ((abs(i) + abs(j)) % 2 == 0)){
								if(circ.getType(i,j) == REGULAR){
									u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - pow2(h)*f(i,j));
								}
								else if(circ.getType(i,j) == IRREGULAR){
									double uN, uS, uE, uW;
									double hN, hS, hE, hW;
									Spacing sp = circ.getSpacing(i,j);
									hN = sp(NORTH);
									hS = sp(SOUTH);
									hE = sp(EAST);
									hW = sp(WEST);
									if(circ.getType(i,j+1) == EXTERIOR){
										uN = circ.boundaryValue(i,j,NORTH);
									}
									else{
										uN = u(i,j+1);
									}
									if(circ.getType(i,j-1) == EXTERIOR){
										uS = circ.boundaryValue(i,j,SOUTH);
									}
									else{
										uS = u(i,j-1);
									}
									if(circ.getType(i+1,j) == EXTERIOR){
										uE = circ.boundaryValue(i,j,EAST);
									}
									else{
										uE = u(i+1,j);
									}
									if(circ.getType(i-1,j) == EXTERIOR){
										uW = circ.boundaryValue(i,j,WEST);
									}
									else{
										uW = u(i-1,j);
									}
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
								}
							}
						}
					}
				}
				return u;
			}
		};*/
	}
}

#endif /* SMOOTHERS_H_ */
