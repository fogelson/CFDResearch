/*
 * StenciledAdvectionDiffusionSolver.h
 *
 *  Created on: Jun 13, 2011
 *      Author: fogelson
 */

#ifndef STENCILEDADVECTIONDIFFUSIONSOLVER_H_
#define STENCILEDADVECTIONDIFFUSIONSOLVER_H_

#include "../Geometry/Geometry.h"
#include "../Geometry/Stencil.h"
#include "../Multigrid/Smoothers.h"
#include "../Multigrid/IntergridOperators.h"
#include "../Multigrid/GridOperators.h"
#include "../Multigrid/MultigridSolvers.h"
#include <lapackpp/lapackpp.h>

namespace CFD{
	using namespace Geometry;
	using namespace Multigrid;
	namespace Solvers{
		class FENEStencil : public Stencil{
			double D, H, Q0, deltaT;
			// Coefficient arrays for diffusion, advection
			Array<Coefficients,2> cDiff, cAdv;

			// Wavespeeds from the FENE spring advective term and fluid velocity
			// advective terms at centroids of cells
			FaceDoubleArray aSpring;

			// Gradient of fluid velocity tensor
			Array<double,2> gradU;

			double magnitude(double q1, double q2){
				return sqrt(pow2(q1) + pow2(q2));
			}

			// Advective speed due to spring forces in the q_1 direction
			double aSq1(double q1, double q2){
				double Q = magnitude(q1,q2);
				return -H*q1/(1 - pow2(Q/Q0));
			}
			double aSq1(Coord c){
				return aSq1(c(0),c(1));
			}
			// And the q_2 direction
			double aSq2(double q1, double q2){
				double Q = magnitude(q1,q2);
				return -H*q2/(1 - pow2(Q/Q0));
			}
			double aSq2(Coord c){
				return aSq2(c(0),c(1));
			}

			double aFq1(double q1, double q2){
				return (gradU(1,1)*q1 + gradU(1,2)*q2);
			}
			double aFq1(Coord c){
				return aFq1(c(0),c(1));
			}
			double aFq2(double q1, double q2){
				return (gradU(2,1)*q1 + gradU(2,2)*q2);
			}
			double aFq2(Coord c){
				return aFq2(c(0),c(1));
			}

			void init(){
				g->resizeArray<Coefficients>(c);
				g->resizeArray<Coefficients>(cDiff);
				g->resizeArray<Coefficients>(cAdv);
				g->resizeArray<FaceDouble>(aSpring);

				double h = g->h;

				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						Coefficients cD, cA;
						FaceDouble aS;

						cD = 0;
						cA = 0;
						aS = 0;
						if(g->isUncovered(i,j)){
							cD += D*(g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g);
							cD += D*(g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g);
							cD -= D*(g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g);
							cD -= D*(g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g);

							cD = (1.0/(h*g->volumeFractions(i,j)))*cD;

							aS(N) = aSq2(g->centroids(i,j)(N));
							aS(S) = aSq2(g->centroids(i,j)(S));
							aS(E) = aSq1(g->centroids(i,j)(E));
							aS(W) = aSq1(g->centroids(i,j)(W));
							aS(B) = 0;
						}
						cDiff(i,j) = cD;
						cAdv(i,j) = 0;
						aSpring(i,j) = aS;
					}
				}
				initAdvection();
			}

			void initAdvection(){
				double h = g->h;
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isUncovered(i,j)){
							Coefficients cA;
							cA = 0;

							FaceDouble a, aF;

							aF(N) = aFq2(g->centroids(i,j)(N));
							aF(S) = aFq2(g->centroids(i,j)(S));
							aF(E) = aFq1(g->centroids(i,j)(E));
							aF(W) = aFq1(g->centroids(i,j)(W));
							aF(B) = 0;

							a = aSpring(i,j) + aF;

							double cE = 0, cW = 0;
							if(g->faceTypes(i,j)(E) != COVERED){
								cE = a(E);
							}
							if(g->faceTypes(i,j)(W) != COVERED){
								cW = a(W);
							}
							if(cE >= 0){
								cA(C) += (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
							}
							else{
								cA(E) += (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
							}
							if(cW >= 0){
								cA(W) += -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
							}
							else{
								cA(C) += -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
							}

							double cN = 0, cS = 0;
							if(g->faceTypes(i,j)(N) != COVERED){
								cN = a(N);
							}
							if(g->faceTypes(i,j)(S) != COVERED){
								cS = a(S);
							}
							if(cN >= 0){
								cA(C) += (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
							}
							else{
								cA(N) += (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
							}
							if(cS >= 0){
								cA(S) += -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
							}
							else{
								cA(C) += -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
							}
							cAdv(i,j) = cA;

							c(i,j) = -deltaT*(cDiff(i,j) - cAdv(i,j));

							c(i,j)(C) = 1 + c(i,j)(C);
						}
						else{
							c(i,j) = 0;
						}
					}
				}
				//cout << "Ran initAdvection()." << endl;
				//c = cDiff + cAdv;
			}
		public:
			FENEStencil(Grid * g, double deltaT, double D, double H, double Q0){
				this->g = g;
				this->deltaT = deltaT;
				this->D = D;
				this->H = H;
				this->Q0 = Q0;
				hasCoarsened = false;

				//cout << "Trivial part of initializing stencil." << endl;

				Range trivialRange(1,2);
				gradU.resize(trivialRange,trivialRange);
				gradU = 0;

				//cout << "Set gradU to zero." << endl;

				init();

				//cout << "Ran init()." << endl;
			}
			void setGradU(Array<double,2> gradU){
				this->gradU = gradU;
				initAdvection();
			}
			Stencil * copy() const{
				FENEStencil * copyPointer = new FENEStencil(g, deltaT, D, H, Q0);
				copyPointer->setGradU(gradU);
				return copyPointer;
			}
		};

		class AdvectionDiffusionStencil : public Stencil{
			double deltaT;
			double D;
			double alpha, H, u11, u12, u21, u22, Q0;

			double cX(double x, double y){
				double r = sqrt(pow2(x) + pow2(y));
				return alpha*(u11*x + u12*y) - H*x/(1 - pow2(r/Q0));
			}
			double cX(Coord c){
				return cX(c(0),c(1));
			}
			double cY(double x, double y){
				double r = sqrt(pow2(x) + pow2(y));
				return alpha*(u21*x + u22*y) - H*y/(1 - pow2(r/Q0));
			}
			double cY(Coord c){
				return cY(c(0),c(1));
			}
			void init(){
				c.resize(g->xRange,g->yRange);
				double h = g->h;
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isUncovered(i,j)){
							Coefficients cD, cA;
							cD = 0;
							cA = 0;

							cD += D*(g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g);
							cD += D*(g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g);
							cD -= D*(g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g);
							cD -= D*(g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g);

							cD = (1.0/(h*g->volumeFractions(i,j)))*cD;

							double cE = 0, cW = 0;
							if(g->faceTypes(i,j)(E) != COVERED){
								cE = cX(g->centroids(i,j)(E));
							}
							if(g->faceTypes(i,j)(W) != COVERED){
								cW = cX(g->centroids(i,j)(W));
							}
							if(cE >= 0){
								cA(C) += (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
							}
							else{
								cA(E) += (g->areaFractions(i,j)(E))*cE/(h*g->volumeFractions(i,j));
							}
							if(cW >= 0){
								cA(W) += -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
							}
							else{
								cA(C) += -(g->areaFractions(i,j)(W))*cW/(h*g->volumeFractions(i,j));
							}

							double cN = 0, cS = 0;
							if(g->faceTypes(i,j)(N) != COVERED){
								cN = cY(g->centroids(i,j)(N));
							}
							if(g->faceTypes(i,j)(S) != COVERED){
								cS = cY(g->centroids(i,j)(S));
							}
							if(cN >= 0){
								cA(C) += (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
							}
							else{
								cA(N) += (g->areaFractions(i,j)(N))*cN/(h*g->volumeFractions(i,j));
							}
							if(cS >= 0){
								cA(S) += -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
							}
							else{
								cA(C) += -(g->areaFractions(i,j)(S))*cS/(h*g->volumeFractions(i,j));
							}
							c(i,j) = -deltaT*(cD - cA);
							c(i,j)(C) = 1 + c(i,j)(C);
						}
						else{
							c(i,j) = 0;
						}
					}
				}
			}
		public:
			AdvectionDiffusionStencil(Grid * g, double deltaT, double D, double alpha, double H, double Q0, double u11, double u12, double u21, double u22){
				this->g = g;
				this->deltaT = deltaT;
				this->D = D;
				this->alpha = alpha;
				this->H = H;
				this->Q0 = Q0;
				this->u11 = u11;
				this->u12 = u12;
				this->u21 = u21;
				this->u22 = u22;
				hasCoarsened = false;
				init();
			}
			Stencil * copy() const{
				return new AdvectionDiffusionStencil(g,deltaT,D,alpha,H,Q0,u11,u12,u21,u22);
			}
		};
		class DiffusionStencil : public Stencil{
			double D, deltaT;
			void init(){
				c.resize(g->xRange,g->yRange);
				double h = g->h;
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						c(i,j) = 0;
						if(g->isUncovered(i,j)){
							c(i,j) += D*(g->areaFractions(i,j)(N))*EBUtilities::getGradientCoefficients(i,j,N,g);
							c(i,j) += D*(g->areaFractions(i,j)(E))*EBUtilities::getGradientCoefficients(i,j,E,g);
							c(i,j) -= D*(g->areaFractions(i,j)(S))*EBUtilities::getGradientCoefficients(i,j,S,g);
							c(i,j) -= D*(g->areaFractions(i,j)(W))*EBUtilities::getGradientCoefficients(i,j,W,g);

							c(i,j) = (1.0/(h*g->volumeFractions(i,j)))*c(i,j);

							c(i,j) = -deltaT*c(i,j);
							c(i,j)(C) = 1 + c(i,j)(C);
						}
					}
				}
			}
		public:
			DiffusionStencil(double D, double deltaT, Grid * g){
				this->D = D;
				this->deltaT = deltaT;
				this->g = g;
				hasCoarsened = false;
				init();
			}
			Stencil * copy() const{
				return new DiffusionStencil(D,deltaT,g);
			}
		};
	}
}
#endif /* STENCILEDADVECTIONDIFFUSIONSOLVER_H_ */
