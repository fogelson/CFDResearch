/*
 * DiffusionSolver.h
 *
 *  Created on: May 4, 2011
 *      Author: fogelson
 */

#ifndef DIFFUSIONSOLVER_H_
#define DIFFUSIONSOLVER_H_

#include "../Geometry/Geometry.h"
#include "../Multigrid/Smoothers.h"
#include "../Multigrid/IntergridOperators.h"
#include "../Multigrid/GridOperators.h"
#include "../Multigrid/MultigridSolvers.h"
#include <lapackpp/lapackpp.h>

using namespace CFD;
using namespace Geometry;
using namespace Multigrid;

using namespace la;

namespace CFD{
	namespace Solvers{
		class DiffusionSolver;

		class DiffusionOperator : public GridOperator{
			double D, deltaT;
		public:
			DiffusionOperator(double D, double deltaT){
				this->D = D;
				this->deltaT = deltaT;
			}
			CellDoubleArray apply(CellDoubleArray u, Grid *g){
				double h = g->h;
				CellDoubleArray Lu = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isRegular(i,j)){
							Lu(i,j) = u(i,j) - (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
						}
						else if(g->isIrregular(i,j)){
							FaceDouble F;
							F = 0;
							F(N) = D*(EBUtilities::getFaceGradient(u,i,j,N,g));
							F(E) = D*(EBUtilities::getFaceGradient(u,i,j,E,g));
							F(S) = D*(EBUtilities::getFaceGradient(u,i,j,S,g));
							F(W) = D*(EBUtilities::getFaceGradient(u,i,j,W,g));
							F(B) = 0; // No flux through boundary

							Lu(i,j) = u(i,j) - (deltaT/(2*h*g->volumeFractions(i,j)))
												*((g->areaFractions(i,j)(N)*F(N))
												-(g->areaFractions(i,j)(S)*F(S))
												+(g->areaFractions(i,j)(E)*F(E))
												-(g->areaFractions(i,j)(W)*F(W))
												);
						}
					}
				}
				return Lu;
			}
			CellDoubleArray applyBackwards(CellDoubleArray u, Grid *g){
				double h = g->h;
				CellDoubleArray Lu = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isRegular(i,j)){
							Lu(i,j) = u(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
						}
						else if(g->isIrregular(i,j)){
							FaceDouble F;
							F = 0;
							F(N) = D*(EBUtilities::getFaceGradient(u,i,j,N,g));
							F(E) = D*(EBUtilities::getFaceGradient(u,i,j,E,g));
							F(S) = D*(EBUtilities::getFaceGradient(u,i,j,S,g));
							F(W) = D*(EBUtilities::getFaceGradient(u,i,j,W,g));
							F(B) = 0; // No flux through boundary

							Lu(i,j) = u(i,j) + (deltaT/(2*h*g->volumeFractions(i,j)))
												*((g->areaFractions(i,j)(N)*F(N))
												-(g->areaFractions(i,j)(S)*F(S))
												+(g->areaFractions(i,j)(E)*F(E))
												-(g->areaFractions(i,j)(W)*F(W))
												);
						}
					}
				}
				return Lu;
			}
		};

		class DiffusionGSRB : public Smoother{
			double D, deltaT;
		public:
			DiffusionGSRB(double D, double deltaT){
				this->D = D;
				this->deltaT = deltaT;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				u = u0;
				rhs = f;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(abs(i + j) % 2 == 0){
								if(g->isRegular(i,j)){
									u(i,j) = (rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)))/(1 + 4*(D*deltaT/(2*pow2(h))));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									u(i,j) = (rhs(i,j)
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
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(abs(i + j) % 2 == 1){
								if(g->isRegular(i,j)){
									u(i,j) = (rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)))/(1 + 4*(D*deltaT/(2*pow2(h))));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									u(i,j) = (rhs(i,j)
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
				}
				return u;
			}
		};


		class DiffusionAltLineZebra : public Smoother{
			double D, deltaT;
		public:
			DiffusionAltLineZebra(double D, double deltaT){
				this->D = D;
				this->deltaT = deltaT;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				u = u0;
				rhs = f;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(abs(j) % 2 == 1){
							Array<double,1> lower(g->xRange), main(g->xRange), upper(g->xRange), right(g->xRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							//LaVectorDouble lower(k-1), center(k), upper(k-1), right(k);
							for(int i = g->iMin; i <= g->iMax; i++){
								if(g->isRegular(i,j)){
									right(i) = rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i,j+1) + u(i,j-1));
									main(i) = 1 - (D*deltaT/(2*pow2(h)))*(-4);
									upper(i) = -(D*deltaT/(2*pow2(h)));
									lower(i) = -(D*deltaT/(2*pow2(h)));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									right(i) = rhs(i,j) - coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
														- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
														- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
														- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
														- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
														- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0);

									main(i) = coefficients(C);
									upper(i) = (g->isUncovered(i+1,j) ? coefficients(E) : 0);
									lower(i) = (g->isUncovered(i-1,j) ? coefficients(W) : 0);
								}
								else if(g->isCovered(i,j)){
									right(i) = 0;
									main(i) = 1;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int i = 0; i < k-1; i++){
								laLower(i) = lower(i + g->iMin + 1);
								laMain(i) = main(i + g->iMin);
								laUpper(i) = upper(i + g->iMin);
								laRight(i) = right(i + g->iMin);
							}
							laMain(k-1) = main(k);
							laRight(k-1) = right(k);
							LaTridiagMatDouble laA(laMain,laLower,laUpper);
							LaTridiagFactDouble laAFactor;
							LaTridiagMatFactorize(laA,laAFactor);
							LaLinearSolve(laAFactor,laU,laRight);
							for(int i = 0; i < k; i++){
								u(i+g->iMin,j) = laU(i);
							}
						}
					}
					for(int j = g->jMin; j <= g->jMax; j++){
						if(abs(j) % 2 == 0){
							Array<double,1> lower(g->xRange), main(g->xRange), upper(g->xRange), right(g->xRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							//LaVectorDouble lower(k-1), center(k), upper(k-1), right(k);
							for(int i = g->iMin; i <= g->iMax; i++){
								if(g->isRegular(i,j)){
									right(i) = rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i,j+1) + u(i,j-1));
									main(i) = 1 - (D*deltaT/(2*pow2(h)))*(-4);
									upper(i) = -(D*deltaT/(2*pow2(h)));
									lower(i) = -(D*deltaT/(2*pow2(h)));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									right(i) = rhs(i,j) - coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
														- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
														- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
														- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
														- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
														- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0);

									main(i) = coefficients(C);
									upper(i) = (g->isUncovered(i+1,j) ? coefficients(E) : 0);
									lower(i) = (g->isUncovered(i-1,j) ? coefficients(W) : 0);
								}
								else if(g->isCovered(i,j)){
									right(i) = 0;
									main(i) = 1;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int i = 0; i < k-1; i++){
								laLower(i) = lower(i + g->iMin + 1);
								laMain(i) = main(i + g->iMin);
								laUpper(i) = upper(i + g->iMin);
								laRight(i) = right(i + g->iMin);
							}
							laMain(k-1) = main(k);
							laRight(k-1) = right(k);
							LaTridiagMatDouble laA(laMain,laLower,laUpper);
							LaTridiagFactDouble laAFactor;
							LaTridiagMatFactorize(laA,laAFactor);
							LaLinearSolve(laAFactor,laU,laRight);
							for(int i = 0; i < k; i++){
								u(i+g->iMin,j) = laU(i);
							}
						}
					}
					for(int i = g->iMin; i <= g->iMax; i++){
						if(abs(i) % 2 == 1){
							Array<double,1> lower(g->yRange), main(g->yRange), upper(g->yRange), right(g->yRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							//LaVectorDouble lower(k-1), center(k), upper(k-1), right(k);
							for(int j = g->jMin; j <= g->jMax; j++){
								if(g->isRegular(i,j)){
									right(j) = rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j));
									main(j) = 1 - (D*deltaT/(2*pow2(h)))*(-4);
									upper(j) = -(D*deltaT/(2*pow2(h)));
									lower(j) = -(D*deltaT/(2*pow2(h)));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									right(j) = rhs(i,j) - coefficients(E)*(g->isUncovered(i+1,j+1) ? u(i+1,j) : 0)
														- coefficients(W)*(g->isUncovered(i-1,j-1) ? u(i-1,j) : 0)
														- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
														- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
														- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
														- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0);

									main(j) = coefficients(C);
									upper(j) = (g->isUncovered(i,j+1) ? coefficients(N) : 0);
									lower(j) = (g->isUncovered(i,j-1) ? coefficients(S) : 0);
								}
								else if(g->isCovered(i,j)){
									right(j) = 0;
									main(j) = 1;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int j = 0; j < k-1; j++){
								laLower(j) = lower(j + g->iMin + 1);
								laMain(j) = main(j + g->iMin);
								laUpper(j) = upper(j + g->iMin);
								laRight(j) = right(j + g->iMin);
							}
							laMain(k-1) = main(k);
							laRight(k-1) = right(k);
							LaTridiagMatDouble laA(laMain,laLower,laUpper);
							LaTridiagFactDouble laAFactor;
							LaTridiagMatFactorize(laA,laAFactor);
							LaLinearSolve(laAFactor,laU,laRight);
							for(int j = 0; j < k; j++){
								u(i,j+g->jMin) = laU(j);
							}
						}
					}
					for(int i = g->iMin; i <= g->iMax; i++){
						if(abs(i) % 2 == 0){
							Array<double,1> lower(g->yRange), main(g->yRange), upper(g->yRange), right(g->yRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							//LaVectorDouble lower(k-1), center(k), upper(k-1), right(k);
							for(int j = g->jMin; j <= g->jMax; j++){
								if(g->isRegular(i,j)){
									right(j) = rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j));
									main(j) = 1 - (D*deltaT/(2*pow2(h)))*(-4);
									upper(j) = -(D*deltaT/(2*pow2(h)));
									lower(j) = -(D*deltaT/(2*pow2(h)));
								}
								else if(g->isIrregular(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;
									coefficients = 0;

									coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
									coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
									coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
									coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

									coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

									right(j) = rhs(i,j) - coefficients(E)*(g->isUncovered(i+1,j+1) ? u(i+1,j) : 0)
														- coefficients(W)*(g->isUncovered(i-1,j-1) ? u(i-1,j) : 0)
														- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
														- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
														- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
														- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0);

									main(j) = coefficients(C);
									upper(j) = (g->isUncovered(i,j+1) ? coefficients(N) : 0);
									lower(j) = (g->isUncovered(i,j-1) ? coefficients(S) : 0);
								}
								else if(g->isCovered(i,j)){
									right(j) = 0;
									main(j) = 1;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int j = 0; j < k-1; j++){
								laLower(j) = lower(j + g->iMin + 1);
								laMain(j) = main(j + g->iMin);
								laUpper(j) = upper(j + g->iMin);
								laRight(j) = right(j + g->iMin);
							}
							laMain(k-1) = main(k);
							laRight(k-1) = right(k);
							LaTridiagMatDouble laA(laMain,laLower,laUpper);
							LaTridiagFactDouble laAFactor;
							LaTridiagMatFactorize(laA,laAFactor);
							LaLinearSolve(laAFactor,laU,laRight);
							for(int j = 0; j < k; j++){
								u(i,j+g->jMin) = laU(j);
							}
						}
					}
				}
				return u;
			}
		};

		class DiffusionXLineLex : public Smoother{
			double D, deltaT;
		public:
			DiffusionXLineLex(double D, double deltaT){
				this->D = D;
				this->deltaT = deltaT;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				u = u0;
				rhs = f;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int j = g->jMin; j <= g->jMax; j++){
						Array<double,1> lower(g->xRange), main(g->xRange), upper(g->xRange), right(g->xRange);
						lower = 0;
						main = 0;
						upper = 0;
						right = 0;
						//LaVectorDouble lower(k-1), center(k), upper(k-1), right(k);
						for(int i = g->iMin; i <= g->iMax; i++){
							if(g->isRegular(i,j)){
								right(i) = rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i,j+1) + u(i,j-1));
								main(i) = 1 - (D*deltaT/(2*pow2(h)))*(-4);
								upper(i) = -(D*deltaT/(2*pow2(h)));
								lower(i) = -(D*deltaT/(2*pow2(h)));
							}
							else if(g->isIrregular(i,j)){
								TinyVector<double,9> coefficients, I;
								I = 0;
								I(C) = 1;
								coefficients = 0;

								coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
								coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
								coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
								coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

								coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

								right(i) = rhs(i,j) - coefficients(N)*(g->isUncovered(i,j+1) ? u(i,j+1) : 0)
													- coefficients(S)*(g->isUncovered(i,j-1) ? u(i,j-1) : 0)
													- coefficients(NE)*(g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0)
													- coefficients(NW)*(g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0)
													- coefficients(SE)*(g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0)
													- coefficients(SW)*(g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0);

								main(i) = coefficients(C);
								upper(i) = (g->isUncovered(i+1,j) ? coefficients(E) : 0);
								lower(i) = (g->isUncovered(i-1,j) ? coefficients(W) : 0);
							}
							else if(g->isCovered(i,j)){
								right(i) = 0;
								main(i) = 1;
							}
						}
						int k = main.size();
						LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
						laU = 0;
						for(int i = 0; i < k-1; i++){
							laLower(i) = lower(i + g->iMin + 1);
							laMain(i) = main(i + g->iMin);
							laUpper(i) = upper(i + g->iMin);
							laRight(i) = right(i + g->iMin);
						}
						laMain(k-1) = main(k);
						laRight(k-1) = right(k);
						LaTridiagMatDouble laA(laMain,laLower,laUpper);
						LaTridiagFactDouble laAFactor;
						LaTridiagMatFactorize(laA,laAFactor);
						LaLinearSolve(laAFactor,laU,laRight);
						for(int i = 0; i < k; i++){
							u(i+g->iMin,j) = laU(i);
						}
					}
				}
				return u;
			}
		};

		class DiffusionGSLex : public Smoother{
			double D, deltaT;
		public:
			DiffusionGSLex(double D, double deltaT){
				this->D = D;
				this->deltaT = deltaT;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				u = u0;
				rhs = f;
				double h = g->h;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(g->isRegular(i,j)){
								//rhs(i,j) = u0(i,j) + (D*deltaT/(2*pow2(h)))*(u0(i+1,j) + u0(i-1,j) + u0(i,j+1) + u0(i,j-1) - 4*u0(i,j));
								u(i,j) = (rhs(i,j) + (D*deltaT/(2*pow2(h)))*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)))/(1 + 4*(D*deltaT/(2*pow2(h))));
							}
							else if(g->isIrregular(i,j)){
//								FaceDouble F0;
//								F0 = 0;
//								F0(N) = D*(EBUtilities::getFaceGradient(u0,i,j,N,g));
//								F0(E) = D*(EBUtilities::getFaceGradient(u0,i,j,E,g));
//								F0(S) = D*(EBUtilities::getFaceGradient(u0,i,j,S,g));
//								F0(W) = D*(EBUtilities::getFaceGradient(u0,i,j,W,g));
//								F0(B) = 0; // No flux through boundary
//
//								rhs(i,j) = u0(i,j) + (deltaT/(2*h*g->volumeFractions(i,j)))
//										*((g->areaFractions(i,j)(N)*F0(N))
//										-(g->areaFractions(i,j)(S)*F0(S))
//										+(g->areaFractions(i,j)(E)*F0(E))
//										-(g->areaFractions(i,j)(W)*F0(W))
//										);
								TinyVector<double,9> coefficients, I;
								I = 0;
								I(C) = 1;
								coefficients = 0;

								coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
								coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
								coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
								coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

								coefficients = I - (D*deltaT/(2*h*g->volumeFractions(i,j)))*coefficients;

								u(i,j) = (rhs(i,j)
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

		class DiffusionSolver{
			double D, deltaT;
			DiffusionOperator * differenceOperator;
			Smoother * smoother;
			MultigridSolver * multigridSolver;
			Interpolator * interpolator;
			Restrictor * restrictor;
			int v1, v2;
		public:
			DiffusionSolver(double D, double deltaT){
				differenceOperator = new DiffusionOperator(D,deltaT);
				smoother = new DiffusionAltLineZebra(D,deltaT);
				interpolator = new BilinearInterpolator();
				restrictor = new VolumeWeightedRestrictor();

				multigridSolver = new MultigridSolver(differenceOperator,smoother,interpolator,restrictor);

				v1 = 2;
				v2 = 2;
			}
			~DiffusionSolver(){
				delete differenceOperator;
				delete smoother;
				delete multigridSolver;
				delete interpolator;
				delete restrictor;
			}
			CellDoubleArray step(CellDoubleArray u0, CellDoubleArray f, Grid * g){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				rhs = differenceOperator->applyBackwards(u0,g);
				double aX = 1, aY = -2;

				u = multigridSolver->solve(u0,rhs,g,v1,v2,8);
				//u = smoother->smooth(u0,f,g,1500);
				return u;
			}
		};
	}
}

#endif /* DIFFUSIONSOLVER_H_ */
