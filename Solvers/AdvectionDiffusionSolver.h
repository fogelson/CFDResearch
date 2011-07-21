/*
 * AdvectionDiffusionSolver.h
 *
 *  Created on: May 9, 2011
 *      Author: fogelson
 */

#ifndef ADVECTIONDIFFUSIONSOLVER_H_
#define ADVECTIONDIFFUSIONSOLVER_H_

#include "../Geometry/Geometry.h"
#include "../Multigrid/Smoothers.h"
#include "../Multigrid/IntergridOperators.h"
#include "../Multigrid/GridOperators.h"
#include "../Multigrid/MultigridSolvers.h"
#include "../Solvers/DiffusionSolver.h"
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
		public:
			static TinyVector<double,9> getCoefficients(double D, CellDoubleArray &u, int i, int j, Grid * g){
				double h = g->h;

				TinyVector<double,9> coefficients;
				coefficients = 0;

				coefficients += (g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(u,i,j,N,g);
				coefficients += (g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(u,i,j,E,g);
				coefficients -= (g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(u,i,j,S,g);
				coefficients -= (g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(u,i,j,W,g);

				coefficients = (D/(h*g->volumeFractions(i,j)))*coefficients;

				return coefficients;
			}
		};
		class AdvectionStencil{
		public:
			static TinyVector<double,9> getCoefficients(double cX, double cY, bool xDir, CellDoubleArray &u, int i, int j, Grid * g){
				double h = g->h;

				TinyVector<double,9> coefficients;
				coefficients = 0;
				if(g->isUncovered(i,j)){
					if(xDir){
						if(cX >= 0){
							coefficients(C) = (g->areaFractions(i,j)(E))*cX/(h*g->volumeFractions(i,j));
							coefficients(W) = -(g->areaFractions(i,j)(W))*cX/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(E) = (g->areaFractions(i,j)(E))*cX/(h*g->volumeFractions(i,j));
							coefficients(C) = -(g->areaFractions(i,j)(W))*cX/(h*g->volumeFractions(i,j));
						}
					}
					else{
						if(cY >= 0){
							coefficients(C) = (g->areaFractions(i,j)(N))*cY/(h*g->volumeFractions(i,j));
							coefficients(S) = -(g->areaFractions(i,j)(S))*cY/(h*g->volumeFractions(i,j));
						}
						else{
							coefficients(N) = (g->areaFractions(i,j)(N))*cY/(h*g->volumeFractions(i,j));
							coefficients(C) = -(g->areaFractions(i,j)(S))*cY/(h*g->volumeFractions(i,j));
						}
					}
				}
				return coefficients;
			}
		};
		class AdvectionDiffusionOperator : public GridOperator{
			double D, deltaT, cX, cY;
			bool xDir;
		public:
			AdvectionDiffusionOperator(double D, double deltaT, double cX, double cY){
				this->D = D;
				this->deltaT = deltaT;
				this->cX = cX;
				this->cY = cY;
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

							coefficients = I + (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
									- (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));
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

							coefficients = I - (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
									+ (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));
							Lu(i,j) = applyCoefficients(coefficients,u,i,j,g);
						}
					}
				}
				return Lu;
			}
		};
		class AdvectionDiffusionAltLineZebra : public Smoother{
			double D, deltaT, cX, cY;
			bool xDir;
		public:
			void toggle(){
				xDir = !xDir;
			}
			AdvectionDiffusionAltLineZebra(double D, double deltaT, double cX, double cY){
				this->D = D;
				this->deltaT = deltaT;
				this->cX = cX;
				this->cY = cY;
				xDir = true;
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
							for(int i = g->iMin; i <= g->iMax; i++){
								if(g->isUncovered(i,j)){
									TinyVector<double,9> coefficients, I;
									I = 0;
									I(C) = 1;

									coefficients = I + (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
											- (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));

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

									coefficients = I + (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
											- (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));

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

									coefficients = I + (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
											- (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));

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

									coefficients = I + (deltaT/2)*(AdvectionStencil::getCoefficients(cX,cY,xDir,u,i,j,g))
											- (deltaT/4)*(DiffusionStencil::getCoefficients(D,u,i,j,g));

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
			double D, deltaT;
			AdvectionDiffusionOperator * differenceOperator;
			AdvectionDiffusionAltLineZebra * smoother;
			MultigridSolver * multigridSolver;
			Interpolator * interpolator;
			Restrictor * restrictor;
			int v1, v2;
		public:
			AdvectionDiffusionSolver(double D, double deltaT, double cX, double cY){
				differenceOperator = new AdvectionDiffusionOperator(D,deltaT,cX,cY);
				smoother = new AdvectionDiffusionAltLineZebra(D,deltaT,cX,cY);
				interpolator = new BilinearInterpolator();
				restrictor = new VolumeWeightedRestrictor();

				multigridSolver = new MultigridSolver(differenceOperator,smoother,interpolator,restrictor);

				v1 = 2;
				v2 = 2;
			}
			~AdvectionDiffusionSolver(){
				delete differenceOperator;
				delete smoother;
				delete multigridSolver;
				delete interpolator;
				delete restrictor;
			}
			CellDoubleArray step(CellDoubleArray u0, CellDoubleArray f, Grid * g){
				differenceOperator->toggle();
				smoother->toggle();


				//cout << "Initializing" << endl;
				CellDoubleArray u = g->makeCellDoubleArray();
				CellDoubleArray rhs = g->makeCellDoubleArray();

				//cout << "Solving x direction 1" << endl;
				rhs = differenceOperator->applyBackwards(u0,g);
				//cout << "Got rhs" << endl;
				u = multigridSolver->solve(u0,rhs,g,v1,v2,2);

				//cout << "Toggling" << endl;
				differenceOperator->toggle();
				smoother->toggle();


				//cout << "Solving y direction 1" << endl;
				rhs = differenceOperator->applyBackwards(u0,g);
				u = multigridSolver->solve(u,rhs,g,v1,v2,2);
				//cout << "Solving y direction 2" << endl;
				rhs = differenceOperator->applyBackwards(u0,g);
				u = multigridSolver->solve(u,rhs,g,v1,v2,2);

				//cout << "Toggling" << endl;
				differenceOperator->toggle();
				smoother->toggle();


				//cout << "Solving x direction 2" << endl;
				rhs = differenceOperator->applyBackwards(u0,g);
				u = multigridSolver->solve(u,rhs,g,v1,v2,2);

				//u = smoother->smooth(u0,f,g,1500);
				return u;
			}
		};
	}
}


#endif /* ADVECTIONDIFFUSIONSOLVER_H_ */
