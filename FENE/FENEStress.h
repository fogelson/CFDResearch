/*
 * FENEStress.h
 *
 *  Created on: Jun 21, 2011
 *      Author: fogelson
 */

#ifndef FENESTRESS_H_
#define FENESTRESS_H_

#include "../Solvers/StenciledAdvectionDiffusionSolver.h"
#include "../Multigrid/Smoothers.h"
#include "MacroscaleObjects.h"
#include <list>
#include <map>
#include <vector>

namespace CFD{
	using namespace Solvers;
	using namespace Geometry;
	using namespace Multigrid;

	class FokkerPlanckSolver{
		Grid * g;
		FENEStencil * stencil;
		Interpolator * interpolator;
		Restrictor * restrictor;
		StenciledSmoother * smoother;
		StenciledMultigridSolver * solver;

		int n;
		double deltaX, deltaY;
		double deltaT, h, offset, D, H, Q0, lambda;

		CellDoubleArray F;

	public:
		Array<double,4> f;
		~FokkerPlanckSolver(){
			delete g;
			delete stencil;
			delete interpolator;
			delete restrictor;
			delete solver;
		}
		FokkerPlanckSolver(int n, double deltaX, double deltaT, double h, double offset, double D, double H, double Q0, double lambda){
			this->n = n;
			this->deltaX = deltaX;
			this->deltaY = deltaX;

			this->deltaT = deltaT;
			this->h = h;
			this->offset = offset;
			this->D = D;
			this->H = H;
			this->Q0 = Q0;
			this->lambda = lambda;

			g = new Circle(h,Q0,offset);
			stencil = new FENEStencil(g,deltaT,D,H,Q0);
			interpolator = new BilinearInterpolator();
			restrictor = new VolumeWeightedRestrictor();
			smoother = new StenciledFourPointGS();
			solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

			F.resize(g->xRange,g->yRange);
			calculateSpringForce(F);

			Range nRange(0,n-1);
			f.resize(nRange,nRange,g->xRange,g->yRange);
			double A = pow2(g->h)*sum(where(g->cellTypes != COVERED, g->volumeFractions, 0.0));
			f = 1.0/A;

			cout << "f has size " << f.shape() << endl;
		}
		void solveFokkerPlanck(Array<double,3> & U){
			Array<double,2> gradU(shape(2,2));
			TinyVector<int,2> gradUIndex;
			gradUIndex(0) = 1;
			gradUIndex(1) = 1;
			gradU.reindexSelf(gradUIndex);

			int iP, iM, jP, jM;
			for(int i = 0; i < n; i++){
				iP = (i + 1 + n) % n;
				iM = (i - 1 + n) % n;
				for(int j = 0; j < n; j++){
					jP = (j + 1 + n) % n;
					jM = (j - 1 + n) % n;

					gradU(1,1) = (U(iP,j,0) - U(iM,j,0))/(2*deltaX); // Du_11
					gradU(1,2) = (U(i,jP,0) - U(i,jM,0))/(2*deltaY); // Du_12
					gradU(2,1) = (U(iP,j,1) - U(iM,j,1))/(2*deltaX); // Du_21
					gradU(2,2) = (U(i,jP,1) - U(i,jM,1))/(2*deltaY); // Du_22

					stencil->setGradU(gradU);
					CellDoubleArray fij = f(i,j,g->xRange,g->yRange);

					//void solve(CellDoubleArray & u, CellDoubleArray & f, Stencil * stencil, int v1, int v2, int its){

					CellDoubleArray rhs = g->makeCellDoubleArray();
					solver->solve(fij,rhs,stencil,2,2,1);
					//fij = solver->solve(fij,fij,stencil,2,2,1);
					f(i,j,g->xRange,g->yRange) = fij;
				}
			}
		}

		void calculateStress(Array<double,3> & S){
			double S11, S12, S22;

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					CellDoubleArray fij = f(i,j,g->xRange,g->yRange);

					stressAtPoint(fij,S11,S12,S22);
					S(i,j,0) = S11;
					S(i,j,1) = S12;
					S(i,j,2) = S22;
				}
			}
		}

		void stressAtPoint(CellDoubleArray &f, double &S11, double &S12, double &S22){
#ifdef FeneTiming
			CFD::Timing::stressAtPoint++;
			time_t begin, end;
			time(&begin);
#endif
			S11 = 0;
			S12 = 0;
			S22 = 0;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					if(g->isUncovered(i,j)){
						//Coord Q = g->centers(i,j);
						Coord Q = g->cellCentroids(i,j);
						double dQ = (g->volumeFractions(i,j))*pow2(g->h);

						S11 += f(i,j)*F(i,j)*pow2(Q(0))*dQ;
						S22 += f(i,j)*F(i,j)*pow2(Q(1))*dQ;
						S12 += f(i,j)*F(i,j)*Q(0)*Q(1)*dQ;
					}
				}
			}
			S11 *= lambda;
			S12 *= lambda;
			S22 *= lambda;
#ifdef FeneTiming
			time(&end);
			CFD::Timing::stressAtPointTime += difftime(end,begin);
#endif
		}
		double magnitude(Coord c){
			return sqrt(pow2(c(0)) + pow2(c(1)));
		}
		void calculateSpringForce(CellDoubleArray &F){
			F = 0;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					if(g->isUncovered(i,j)){
						//double Q = magnitude(g->centers(i,j));//magnitude(g->cellCentroids(i,j));
						//Q = min(Q0-h/10.0,Q);
						double Q = magnitude(g->cellCentroids(i,j));
						F(i,j) = H*Q;//H/(1 - pow2(Q/Q0));
						//F(i,j) = H*magnitude(g->centers(i,j));
					}
				}
			}
		}
		Grid * getGrid(){
			return g;
		}
	};

}

#endif /* FENESTRESS_H_ */
