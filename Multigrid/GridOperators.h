/*
 * GridOperators.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef GRIDOPERATORS_H_
#define GRIDOPERATORS_H_

#include "../Geometry/Geometry.h"
//#include "../LinearAlgebra/LinearAlgebra.h"
//using namespace LinearAlgebra;

/* Creates a number of GridOperator child classes for difference
 * operators associated with discretized PDEs. The GridOperator
 * abstract parent class is declared in Geometry.h
 */

namespace CFD{
	using namespace Geometry;
	namespace Multigrid{
		class GridOperator;
		class DifferenceOperator;

		class LaplaceOperator;

		class GridOperator{
		public:
			virtual ~GridOperator(){}
			virtual CellDoubleArray apply(CellDoubleArray u, Grid *g) = 0;
		};

		class LaplaceOperator : public GridOperator{
		public:
			CellDoubleArray apply(CellDoubleArray u, Grid *g){
				double h = g->h;
				CellDoubleArray Lu = g->makeCellDoubleArray();
				for(int i = g->iMin; i <= g->iMax; i++){
					for(int j = g->jMin; j <= g->jMax; j++){
						if(g->isRegular(i,j)){
							Lu(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j))/pow2(h);
						}
						else if(g->isIrregular(i,j)){
							FaceDouble F;
							F = 0;
							F(N) = (g->areaFractions(i,j)(N))*h*(EBUtilities::getFaceGradient(u,i,j,N,g));
							F(E) = (g->areaFractions(i,j)(E))*h*(EBUtilities::getFaceGradient(u,i,j,E,g));
							F(S) = (g->areaFractions(i,j)(S))*h*(EBUtilities::getFaceGradient(u,i,j,S,g));
							F(W) = (g->areaFractions(i,j)(W))*h*(EBUtilities::getFaceGradient(u,i,j,W,g));
							F(B) = 0; // No flux through boundary

							Lu(i,j) = (F(E) - F(W) + F(N) - F(S) - F(B))/(pow2(h)*(g->volumeFractions(i,j)));
						}
					}
				}
				return Lu;
			}
		};


		/* Applies the 2D discrete Laplacian with Dirichlet boundary values
		 * given by the function pointed to by boundaryFunction to GridScalar
		 * u.
		 */

		// Discrete Laplacian with standard stencil.
/*		class LaplaceOperator : public GridOperator{
		public:
			virtual ~LaplaceOperator(){}
			virtual GridScalar apply(GridScalar u, Circle circ){
				double h = circ.getH();
				Coord center = circ.getCenter();
				GridScalar Lu = circ.makeScalar();
				for(int i = Lu.lbound(0); i <= Lu.ubound(0); i++){
					for(int j = Lu.lbound(1); j <= Lu.ubound(1); j++){
						if(circ.getType(i,j) == REGULAR){
							Lu(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4.0*u(i,j))/pow2(h);
						}
						else if(circ.getType(i,j) == IRREGULAR){
							Spacing sp = circ.getSpacing(i,j);
							double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
							NeighborScalar nbr = circ.getNeighbors(u,i,j);
							double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
							Lu(i,j) = 2.0*((uN/hN + uS/hS)/(hN + hS) + (uE/hE + uW/hW)/(hE + hW) - u(i,j)/(hW*hE) - u(i,j)/(hN*hS));
						}
					}
				}
				return Lu;
			}
		};*/
/*
		class LaxWendroffOperator : public GridOperator{
			double aX, aY;
		public:
			virtual ~LaxWendroffOperator(){}
			LaxWendroffOperator(){}
			LaxWendroffOperator(double aX, double aY){
				this->aX = aX;
				this->aY = aY;
			}
			virtual GridScalar apply(GridScalar u, Circle circ){
				double h = circ.getH();
				Coord center = circ.getCenter();
				GridScalar Lu = circ.makeScalar();
				for(int i = Lu.lbound(0); i <= Lu.ubound(0); i++){
					for(int j = Lu.lbound(1); j <= Lu.ubound(1); j++){
						if(circ.getType(i,j) == REGULAR){
							Lu(i,j) = aX*(u(i+1,j) - u(i-1,j))/(2.0*h) + pow2(aX/h)*(u(i-1,j) - 2.0*u(i,j) + u(i+1,j))/(2.0)
									+ aY*(u(i,j+1) - u(i,j-1))/(2.0*h) + pow2(aY/h)*(u(i,j-1) - 2.0*u(i,j) + u(i,j+1))/(2.0);
						}
						else if(circ.getType(i,j) == IRREGULAR){
							Spacing sp = circ.getSpacing(i,j);
							double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
							NeighborScalar nbr = circ.getNeighbors(u,i,j);
							double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
							Lu(i,j) = aX*(uE - uW)/(hE + hW) + pow2(aX)*((uE/hE + uW/hW)/(hE + hW) - u(i,j)/(hE*hW))
									+ aY*(uN - uS)/(hN + hS) + pow2(aY)*((uN/hN + uS/hS)/(hN + hS) - u(i,j)/(hN*hS));
						}
					}
				}
				return Lu;
			}
		};

		class UpwindOperator : public GridOperator{
			double aX, aY;
		public:
			virtual ~UpwindOperator(){}
			UpwindOperator(){}
			UpwindOperator(double aX, double aY){
				this->aX = aX;
				this->aY = aY;
			}
			virtual GridScalar apply(GridScalar u, Circle circ){
				double h = circ.getH();
				Coord center = circ.getCenter();
				GridScalar Lu = circ.makeScalar();
				for(int i = Lu.lbound(0); i <= Lu.ubound(0); i++){
					for(int j = Lu.lbound(1); j <= Lu.ubound(1); j++){
						if(circ.getType(i,j) == REGULAR){
							if(aX > 0){
								Lu(i,j) += aX*(u(i,j) - u(i-1,j))/h;
							}
							else{
								Lu(i,j) += aX*(u(i+1,j) - u(i,j))/h;
							}
							if(aY > 0){
								Lu(i,j) += aY*(u(i,j) - u(i,j-1))/h;
							}
							else{
								Lu(i,j) += aY*(u(i,j+1) - u(i,j))/h;
							}
						}
						else if(circ.getType(i,j) == IRREGULAR){
							Spacing sp = circ.getSpacing(i,j);
							double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
							NeighborScalar nbr = circ.getNeighbors(u,i,j);
							double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
							if(aX > 0){
								Lu(i,j) += aX*(u(i,j) - uW)/hW;
							}
							else{
								Lu(i,j) += aX*(uE - u(i,j))/hE;
							}
							if(aY > 0){
								Lu(i,j) += aY*(u(i,j) - uS)/hS;
							}
							else{
								Lu(i,j) += aY*(uN - u(i,j))/hN;
							}
						}
					}
				}
				return Lu;
			}
		};



		class AdvectionDiffusionOperator : public GridOperator{
			double D, aX, aY;
			UpwindOperator UW;
			LaplaceOperator L;
		public:
			virtual ~AdvectionDiffusionOperator(){}
			AdvectionDiffusionOperator(double D, double aX, double aY){
				this->D = D;
				this->aX = aX;
				this->aY = aY;
				UW = UpwindOperator(aX,aY);
			}
			virtual GridScalar apply(GridScalar u, Circle circ){
				GridScalar Lu = circ.makeScalar();
				Lu = UW.apply(u,circ) + D*L.apply(u,circ);
				return Lu;
			}
		};*/
	}
}

#endif /* GRIDOPERATORS_H_ */
