/*
 * Smoothers.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef SMOOTHERS_H_
#define SMOOTHERS_H_

#include "../Geometry.h"
#include "../LinearAlgebra/Krylov.h"

//#include <lapackpp/lapackpp.h>

//using namespace la;
using namespace CFD;
using namespace Geometry;

namespace CFD{
	namespace Multigrid{
		class Smoother;

		/* Abstract base class for iterative solvers the smooth the solution
		 * quickly.
		 */
		class Smoother{
		public:
			virtual ~Smoother(){}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its) = 0;
		};

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
		};
	}
}

#endif /* SMOOTHERS_H_ */
