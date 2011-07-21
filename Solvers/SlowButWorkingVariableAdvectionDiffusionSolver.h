/*
 * VariableAdvectionDiffusionSolver.h
 *
 *  Created on: May 10, 2011
 *      Author: fogelson
 */

#ifndef VARIABLEADVECTIONDIFFUSIONSOLVER_H_
#define VARIABLEADVECTIONDIFFUSIONSOLVER_H_

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
	/*bool isXDir(){
		return xDir;
	}
	void setXDir(bool xDir){
		this->xDir = xDir;
	}
	void toggle(){
		xDir = !xDir;
	}*/

		class AdvectionData{
		public:
			double H, Q0, alpha, u11, u12, u21, u22, D;
		};


		double applyCoefficients(TinyVector<double,9> coefficients, CellDoubleArray &u, int i, int j, Grid * g){
			double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
			uC = g->isUncovered(i,j) ? u(i,j) : 0;
			uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
			uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
			uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
			uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
			uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
			uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
			uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
			uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

			double Lu = coefficients(C)*uC
					+ coefficients(N)*uN
					+ coefficients(S)*uS
					+ coefficients(E)*uE
					+ coefficients(W)*uW
					+ coefficients(NW)*uNW
					+ coefficients(NE)*uNE
					+ coefficients(SW)*uSW
					+ coefficients(SE)*uSE;
			return Lu;
		}
		class DiffusionStencil{
			static double D(double x, double y, AdvectionData &A){
				return A.D;
			}
			static double D(Coord c, AdvectionData &A){
				return D(c(0),c(1),A);
			}
		public:
			static TinyVector<double,9> getCoefficients(CellDoubleArray &u, int i, int j, Grid * g, AdvectionData &A){
				double h = g->h;

				TinyVector<double,9> coefficients;
				coefficients = 0;

				coefficients += (D(g->centroids(i,j)(N),A))*(g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
				coefficients += (D(g->centroids(i,j)(E),A))*(g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
				coefficients -= (D(g->centroids(i,j)(S),A))*(g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
				coefficients -= (D(g->centroids(i,j)(W),A))*(g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

				coefficients = (1.0/(h*g->volumeFractions(i,j)))*coefficients;

				return coefficients;
			}
		};
		class AdvectionStencil{
			static double cX(double x, double y, AdvectionData &A){
				double r = sqrt(pow2(x) + pow2(y));
				return A.alpha*(A.u11*x + A.u12*y) - A.H*x/(1 - pow2(r/A.Q0));
			}
			static double cX(Coord c, AdvectionData &A){
				return cX(c(0),c(1),A);
			}
			static double cY(double x, double y, AdvectionData &A){
				double r = sqrt(pow2(x) + pow2(y));
				return A.alpha*(A.u21*x + A.u22*y) - A.H*y/(1 - pow2(r/A.Q0));
			}
			static double cY(Coord c, AdvectionData &A){
				return cY(c(0),c(1),A);
			}
		public:
			static TinyVector<double,9> getCoefficients(bool xDir, CellDoubleArray &u, int i, int j, Grid * g, AdvectionData &A){
				double h = g->h;

				TinyVector<double,9> coefficients;
				coefficients = 0;
				if(g->isUncovered(i,j)){
					if(xDir){
						double cE = 0, cW = 0;
						if(g->faceTypes(i,j)(E) != COVERED){
							cE = cX(g->centroids(i,j)(E),A);
						}
						if(g->faceTypes(i,j)(W) != COVERED){
							cW = cX(g->centroids(i,j)(W),A);
						}
						if(cE >= 0){
							coefficients(C) = (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(E) = (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
						}
						if(cW >= 0){
							coefficients(W) = -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(C) = -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
						}
					}
					else{
						double cN = 0, cS = 0;
						if(g->faceTypes(i,j)(N) != COVERED){
							cN = cY(g->centroids(i,j)(N),A);
						}
						if(g->faceTypes(i,j)(S) != COVERED){
							cS = cY(g->centroids(i,j)(S),A);
						}
						if(cN >= 0){
							coefficients(C) = (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(N) = (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
						}
						if(cS >= 0){
							coefficients(S) = -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(C) = -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
						}
					}
				}
				return coefficients;
			}
		};
		class AdvectionDiffusionStencil{
		public:
			static TinyVector<double,9> implicitCoefficients(bool xDir, CellDoubleArray &u, int i, int j, Grid * g, AdvectionData &A){
				TinyVector<double,9> coefficients;
				/*// FE advection
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g))/4;*/

				/*// BE advection
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g))/4 - (AdvectionStencil::getCoefficients(xDir,u,i,j,g));*/

				// Unsplit BE
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g,A))/4
						- (AdvectionStencil::getCoefficients(xDir,u,i,j,g,A))/2
						- (AdvectionStencil::getCoefficients(!xDir,u,i,j,g,A))/2;
				return coefficients;
			}
			static TinyVector<double,9> explicitCoefficients(bool xDir, CellDoubleArray &u, int i, int j, Grid * g, AdvectionData &A){
				TinyVector<double,9> coefficients;
				/*// FE advection
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g))/4 - (AdvectionStencil::getCoefficients(xDir,u,i,j,g));*/

				/*// BE advection
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g))/4;*/

				// Unsplit BE
				coefficients = (DiffusionStencil::getCoefficients(u,i,j,g,A))/4;
				return coefficients;
			}
		};
		class AdvectionDiffusionOperator : public GridOperator{
			double deltaT;
			bool xDir;
			AdvectionData A;
		public:
			AdvectionDiffusionOperator(double deltaT, AdvectionData &A){
				this->deltaT = deltaT;
				this->A = A;
			}
			void toggle(){
				xDir = !xDir;
			}
			CellDoubleArray apply(CellDoubleArray u, Grid * g){
				CellDoubleArray Lu = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isUncovered(i,j)){
							TinyVector<double,9> coefficients, I;
							I = 0;
							I(C) = 1;
							coefficients = 0;

							coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));
							Lu(i,j) = applyCoefficients(coefficients,u,i,j,g);
						}
					}
				}
				return Lu;
			}
			CellDoubleArray applyBackwards(CellDoubleArray u, Grid * g){
				CellDoubleArray Lu = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isUncovered(i,j)){
							TinyVector<double,9> coefficients, I;
							I = 0;
							I(C) = 1;
							coefficients = 0;

							coefficients = I + deltaT*(AdvectionDiffusionStencil::explicitCoefficients(xDir,u,i,j,g,A));
							Lu(i,j) = applyCoefficients(coefficients,u,i,j,g);
						}
					}
				}
				return Lu;
			}
		};
		class AdvectionDiffusionFourPointGSLex : public Smoother{
			double deltaT;
			bool xDir;
			AdvectionData A;
		public:
			void toggle(){
				xDir = !xDir;
			}
			AdvectionDiffusionFourPointGSLex(double deltaT, AdvectionData &A){
				this->deltaT = deltaT;
				xDir = true;
				this->A = A;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its){
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();
				u = u0;
				rhs = f;
				double h = g->h;
				TinyVector<double,9> coefficients, I;
				coefficients = 0;
				I = 0;
				I(C) = 1;
				for(int n = 1; n <= its; n++){
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(g->isUncovered(i,j)){
								coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));
								double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
								uC = g->isUncovered(i,j) ? u(i,j) : 0;
								uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
								uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
								uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
								uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
								uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
								uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
								uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
								uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

								u(i,j) = rhs(i,j)
										- coefficients(N)*uN
										- coefficients(S)*uS
										- coefficients(E)*uE
										- coefficients(W)*uW
										- coefficients(NE)*uNE
										- coefficients(NW)*uNW
										- coefficients(SE)*uSE
										- coefficients(SW)*uSW;
								u(i,j) = u(i,j)/coefficients(C);
							}
						}
					}
					for(int i = g->iMax; i >= g->iMin; i--){
						for(int j = g->jMin; j <= g->jMax; j++){
							if(g->isUncovered(i,j)){
								coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));
								double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
								uC = g->isUncovered(i,j) ? u(i,j) : 0;
								uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
								uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
								uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
								uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
								uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
								uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
								uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
								uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

								u(i,j) = rhs(i,j)
										- coefficients(N)*uN
										- coefficients(S)*uS
										- coefficients(E)*uE
										- coefficients(W)*uW
										- coefficients(NE)*uNE
										- coefficients(NW)*uNW
										- coefficients(SE)*uSE
										- coefficients(SW)*uSW;
								u(i,j) = u(i,j)/coefficients(C);
							}
						}
					}
					for(int i = g->iMin; i <= g->iMax; i++){
						for(int j = g->jMax; j >= g->jMin; j--){
							if(g->isUncovered(i,j)){
								coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));
								double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
								uC = g->isUncovered(i,j) ? u(i,j) : 0;
								uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
								uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
								uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
								uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
								uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
								uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
								uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
								uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

								u(i,j) = rhs(i,j)
										- coefficients(N)*uN
										- coefficients(S)*uS
										- coefficients(E)*uE
										- coefficients(W)*uW
										- coefficients(NE)*uNE
										- coefficients(NW)*uNW
										- coefficients(SE)*uSE
										- coefficients(SW)*uSW;
								u(i,j) = u(i,j)/coefficients(C);
							}
						}
					}
					for(int i = g->iMax; i >= g->iMin; i--){
						for(int j = g->jMax; j >= g->jMin; j--){
							if(g->isUncovered(i,j)){
								coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));
								double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
								uC = g->isUncovered(i,j) ? u(i,j) : 0;
								uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
								uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
								uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
								uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
								uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
								uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
								uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
								uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

								u(i,j) = rhs(i,j)
										- coefficients(N)*uN
										- coefficients(S)*uS
										- coefficients(E)*uE
										- coefficients(W)*uW
										- coefficients(NE)*uNE
										- coefficients(NW)*uNW
										- coefficients(SE)*uSE
										- coefficients(SW)*uSW;
								u(i,j) = u(i,j)/coefficients(C);
							}
						}
					}
				}
				return u;
			}
		};
		class AdvectionDiffusionAltLineZebra : public Smoother{
			double deltaT;
			bool xDir;
			AdvectionData A;
		public:
			void toggle(){
				xDir = !xDir;
			}
			AdvectionDiffusionAltLineZebra(double deltaT, AdvectionData &A){
				this->deltaT = deltaT;
				xDir = true;
				this->A = A;
			}
			CellDoubleArray smooth(CellDoubleArray u0, CellDoubleArray f, Grid *g, int its, AdvectionData &A){
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
							for(int i = g->iMin; i <= g->iMax; i++){
								if(g->isUncovered(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;

									coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));

									double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
									uC = g->isUncovered(i,j) ? u(i,j) : 0;
									uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
									uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
									uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
									uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
									uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
									uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
									uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
									uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

									right(i) = rhs(i,j)
											- coefficients(N)*uN
											- coefficients(S)*uS
											- coefficients(NE)*uNE
											- coefficients(NW)*uNW
											- coefficients(SE)*uSE
											- coefficients(SW)*uSW;

									lower(i) = coefficients(W);
									main(i) = coefficients(C);
									upper(i) = coefficients(E);
								}
								else{
									main(i) = 1;
									right(i) = 0;
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
							for(int i = g->iMin; i <= g->iMax; i++){
								if(g->isUncovered(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;

									coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));

									double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
									uC = g->isUncovered(i,j) ? u(i,j) : 0;
									uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
									uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
									uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
									uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
									uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
									uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
									uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
									uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

									right(i) = rhs(i,j)
											- coefficients(N)*uN
											- coefficients(S)*uS
											- coefficients(NE)*uNE
											- coefficients(NW)*uNW
											- coefficients(SE)*uSE
											- coefficients(SW)*uSW;

									lower(i) = coefficients(W);
									main(i) = coefficients(C);
									upper(i) = coefficients(E);
								}
								else{
									main(i) = 1;
									right(i) = 0;
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
							Array<double,1> lower(g->xRange), main(g->xRange), upper(g->xRange), right(g->xRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							for(int j = g->jMin; j <= g->jMax; j++){
								if(g->isUncovered(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;

									coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));

									double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
									uC = g->isUncovered(i,j) ? u(i,j) : 0;
									uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
									uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
									uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
									uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
									uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
									uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
									uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
									uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

									right(j) = rhs(i,j)
											- coefficients(E)*uE
											- coefficients(W)*uW
											- coefficients(NE)*uNE
											- coefficients(NW)*uNW
											- coefficients(SE)*uSE
											- coefficients(SW)*uSW;

									lower(j) = coefficients(S);
									main(j) = coefficients(C);
									upper(j) = coefficients(N);
								}
								else{
									main(j) = 1;
									right(j) = 0;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int j = 0; j < k-1; j++){
								laLower(j) = lower(j + g->jMin + 1);
								laMain(j) = main(j + g->jMin);
								laUpper(j) = upper(j + g->jMin);
								laRight(j) = right(j + g->jMin);
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
							Array<double,1> lower(g->xRange), main(g->xRange), upper(g->xRange), right(g->xRange);
							lower = 0;
							main = 0;
							upper = 0;
							right = 0;
							for(int j = g->jMin; j <= g->jMax; j++){
								if(g->isUncovered(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;

									coefficients = I - deltaT*(AdvectionDiffusionStencil::implicitCoefficients(xDir,u,i,j,g,A));

									double uC, uN, uS, uE, uW, uNE, uNW, uSE, uSW;
									uC = g->isUncovered(i,j) ? u(i,j) : 0;
									uN = g->isUncovered(i,j+1) ? u(i,j+1) : 0;
									uS = g->isUncovered(i,j-1) ? u(i,j-1) : 0;
									uE = g->isUncovered(i+1,j) ? u(i+1,j) : 0;
									uW = g->isUncovered(i-1,j) ? u(i-1,j) : 0;
									uNW = g->isUncovered(i-1,j+1) ? u(i-1,j+1) : 0;
									uNE = g->isUncovered(i+1,j+1) ? u(i+1,j+1) : 0;
									uSW = g->isUncovered(i-1,j-1) ? u(i-1,j-1) : 0;
									uSE = g->isUncovered(i+1,j-1) ? u(i+1,j-1) : 0;

									right(j) = rhs(i,j)
											- coefficients(E)*uE
											- coefficients(W)*uW
											- coefficients(NE)*uNE
											- coefficients(NW)*uNW
											- coefficients(SE)*uSE
											- coefficients(SW)*uSW;

									lower(j) = coefficients(S);
									main(j) = coefficients(C);
									upper(j) = coefficients(N);
								}
								else{
									main(j) = 1;
									right(j) = 0;
								}
							}
							int k = main.size();
							LaVectorDouble laLower(k-1), laMain(k), laUpper(k-1), laRight(k), laU(k);
							laU = 0;
							for(int j = 0; j < k-1; j++){
								laLower(j) = lower(j + g->jMin + 1);
								laMain(j) = main(j + g->jMin);
								laUpper(j) = upper(j + g->jMin);
								laRight(j) = right(j + g->jMin);
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

		class AdvectionDiffusionSolver{
			AdvectionData A;
			double deltaT;
			AdvectionDiffusionOperator * differenceOperator;
			//AdvectionDiffusionAltLineZebra * smoother;
			AdvectionDiffusionFourPointGSLex * smoother;
			MultigridSolver * multigridSolver;
			Interpolator * interpolator;
			Restrictor * restrictor;
			int v1, v2;
		public:
			AdvectionDiffusionSolver(double deltaT, AdvectionData &A){
				this->A = A;
				differenceOperator = new AdvectionDiffusionOperator(deltaT,A);
				//smoother = new AdvectionDiffusionAltLineZebra(deltaT,A);
				smoother = new AdvectionDiffusionFourPointGSLex(deltaT,A);
				interpolator = new BilinearInterpolator();
				restrictor = new VolumeWeightedRestrictor();

				multigridSolver = new MultigridSolver(differenceOperator,smoother,interpolator,restrictor);

				v1 = 2;
				v2 = 1;

				differenceOperator->toggle();
				smoother->toggle();
			}
			~AdvectionDiffusionSolver(){
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
				u = multigridSolver->solve(u0,rhs,g,v1,v2,4);

				differenceOperator->toggle();
				smoother->toggle();

				rhs = differenceOperator->applyBackwards(u0,g);
				u = multigridSolver->solve(u,rhs,g,v1,v2,4);

				//u = smoother->smooth(u0,f,g,1500);
				return u;
			}
		};
	}
}

#endif /* VARIABLEADVECTIONDIFFUSIONSOLVER_H_ */
