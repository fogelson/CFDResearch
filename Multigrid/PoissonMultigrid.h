/*
 * PoissonMultigrid.h
 *
 *  Created on: Feb 6, 2011
 *      Author: fogelson
 */

/* Methods for solving Poisson's equation on a circle using
 * a Multigrid solver.
 */

#ifndef POISSONMULTIGRID_H_
#define POISSONMULTIGRID_H_

#include "../Geometry.h"

using namespace CFD::Geometry;

namespace CFD{
	namespace Multigrid{
		/* Many of the later methods require knowing the values on
		 * the boundary. For all except the top level in a Dirichlet
		 * Poisson problem, the function will be zero on the boundary,
		 * so we just build in a function to always return 0.
		 */
		double zeroBoundary(double x, double y){
			return 0;
		}

		// Norms are useful, and the maxNorm is really easy to code.
		double maxNorm(GridScalar u, Circle circ){
			return max(where(circ.getTypes() != EXTERIOR, abs(u), 0));
		}

		/* Use bilinear interpolation to go from a GridScalar uC on a grid
		 * a coarse grid with spacing hC to a fine grid with spacing hF = hC/2.
		 */
		GridScalar BilinearInterpolation(GridScalar uC, Circle fineCircle, Circle coarseCircle, double (*boundary)(double,double)){
			double h = fineCircle.getH();
			GridScalar uF = fineCircle.makeScalar();

			/* Iterate over the fine grid. This is different than
			 * when we are on a normal rectangular grid, because here
			 * we have to treat the irregular points on the fine grid
			 * extremely carefully, so it makes a lot more sense to
			 * just iterate over that grid.
			 */
			for(int iF = uF.lbound(0); iF <= uF.ubound(0); iF++){
				for(int jF = uF.lbound(1); jF <= uF.ubound(1); jF++){
					if(fineCircle.getType(iF,jF) != EXTERIOR){
						// If the fine point overlaps a coarse point, set equal
						if(iF % 2 == 0 && jF % 2 == 0){
							uF(iF,jF) = uC(iF/2,jF/2);
						}
						// Interpolate in just the y-direction for points along a coarse
						// grid line in the x-direction.
						else if(iF % 2 == 0){
							if(fineCircle.getType(iF,jF) == REGULAR){
								uF(iF,jF) = 0.5*(uC(iF/2,(jF+1)/2) + uC(iF/2,(jF-1)/2));
							}
							/* If the point is irregular, we need to get the value of the function
							 * on the boundary of the domain, and do a linear interpolation using
							 * the boundary point instead of an adjacent grid point.
							 */
							else{
								// UN and US are the north and south values we will
								// use to interpolate. One of these values should just
								// be a grid value from the coarse grid, and one of them
								// will be a boundary value from the Dirichlet data.
								double UN, US, hN, hS;
								hN = fineCircle.getSpacing(iF,jF)(NORTH);
								hS = fineCircle.getSpacing(iF,jF)(SOUTH);
								if(fineCircle.getType(iF,jF+1) != EXTERIOR){
									UN = uC(iF/2,(jF+1)/2);
								}
								else{
									Coord b = fineCircle.nearestBoundary(iF,jF,NORTH);
									UN = (*boundary)(b(0),b(1));
								}
								if(fineCircle.getType(iF,jF-1) != EXTERIOR){
									US = uC(iF/2,(jF-1)/2);
								}
								else{
									Coord b = fineCircle.nearestBoundary(iF,jF,SOUTH);
									US = (*boundary)(b(0),b(1));
								}
								uF(iF,jF) = (hS*UN + hN*US)/(hS + hN);
							}
						}
						// Interpolate in just the x-direction for points along a coarse
						// grid line in the y-direction. See the above code for the
						// analogous case if you want detailed comments.
						else if(jF % 2 == 0){
							if(fineCircle.getType(iF,jF) == REGULAR){
								uF(iF,jF) = 0.5*(uC((iF+1)/2,jF/2) + uC((iF-1)/2,jF/2));
							}
							else{
								double UE, UW, hE, hW;
								hE = fineCircle.getSpacing(iF,jF)(EAST);
								hW = fineCircle.getSpacing(iF,jF)(WEST);
								if(fineCircle.getType(iF+1,jF) != EXTERIOR){
									UE = uC((iF+1)/2,jF/2);
								}
								else{
									Coord b = fineCircle.nearestBoundary(iF,jF,EAST);
									UE = (*boundary)(b(0),b(1));
								}
								if(fineCircle.getType(iF-1,jF) != EXTERIOR){
									UW = uC((iF-1)/2,jF/2);
								}
								else{
									Coord b = fineCircle.nearestBoundary(iF,jF,WEST);
									UW = (*boundary)(b(0),b(1));
								}
								uF(iF,jF) = (hW*UE + hE*UW)/(hW + hE);
							}
						}
						/* Interpolate in both directions. We treat both regular and irregular
						 * gridpoints the same here, because it is possible for a point on the
						 * fine grid to be regular while one of its diagonal neighbors on the
						 * coarse grid is outside the domain. So just be equally careful for
						 * all of them.
						 *
						 * This interpolation is built on the fact that bilinear interpolation
						 * can be thought of as the tensor product of interpolation in the x
						 * and y directions.
						 *
						 * 						 uNW  -------  uNE
						 * 						  |				|
						 * 						  |				|
						 * 						  |		 u		|
						 * 						  |				|
						 * 						  |				|
						 * 						uSW  -------   uSE
						 *
						 * In the above diagram, we are trying to interpolate the value u
						 * on the fine grid, when we know the values uNW, uNE, uSW, and
						 * uSE on the coarse grid. If everything is regularly spaced, all
						 * we have to do is linear interpolation in the x-direction above and
						 * below u. This means interpolating uNW and uNE to get a north
						 * value, and uSW and uSE to get a south value. Then we just
						 * interpolate between the north and south values in the y-direction
						 * to get u(F).
						 *
						 * For the regular, rectangular case, this bilinear interpolation is
						 * independent of whether we interpolate in the x-direction then
						 * the y-direction, or vice-versa. When we are near a domain boundary,
						 * however, this is not the case, and our algorithm requires some
						 * modification.
						 *
						 * First, we interpolate in the x-direction above and below u.
						 * If both the SW and SE points are interior, we interpolate as
						 * normal. If one of them is exterior, we just have to make sure
						 * our interpolation is weighted correctly. If both of them are
						 * exterior, rather than interpolating we just use the Dirichlet
						 * data to come up with the southern point (which will be closer
						 * to u than a full grid space). We do similarly for the NW
						 * and NE points to get a northern point, and then interpolate
						 * those to get an estimate for u.
						 *
						 * Note, however, that this estimate for u is actually different
						 * than if we had first interpolated in the y-direction to get
						 * east and west points, and then interpolated those in the x-direction
						 * to get u. So we actually do both. We get one value uXY by
						 * interpolating in x then y, and another uYX by interpolating in
						 * y then x, and we average these.
						 */
						else{
							double hN, hS, hE, hW;
							double uN, uS, uE, uW;
							double uXY, uYX;

							// Interpolate in the x-direction to get north and south points

							// Get uS

							// If we run into the boundary before going a full grid space
							// south, just use Dirichlet data to get uS.
							if(fineCircle.getType(iF,jF-1) == EXTERIOR){
								uS = fineCircle.evaluateOnNearestBoundary(iF,jF,SOUTH,(*boundary));
								hS = fineCircle.getSpacing(iF,jF)(SOUTH);
							}
							// Interpolate between the SW and SE points, using boundary
							// points as needed.
							else{
								double hSW, hSE;
								double uSW, uSE;
								if(fineCircle.getType(iF-1,jF-1) == EXTERIOR){
									uSW = fineCircle.evaluateOnNearestBoundary(iF,jF-1,WEST,(*boundary));
									hSW = fineCircle.getSpacing(iF,jF-1)(WEST);
								}
								else{
									uSW = uC((iF-1)/2,(jF-1)/2);
									hSW = h;
								}
								if(fineCircle.getType(iF+1,jF-1) == EXTERIOR){
									uSE = fineCircle.evaluateOnNearestBoundary(iF,jF-1,EAST,(*boundary));
									hSE = fineCircle.getSpacing(iF,jF-1)(EAST);
								}
								else{
									uSE = uC((iF+1)/2,(jF-1)/2);
									hSE = h;
								}
								uS = uSE*hSW/(hSE + hSW) + uSW*hSE/(hSE + hSW);
								hS = h;
							}

							// Get uN
							if(fineCircle.getType(iF,jF+1) == EXTERIOR){
								uN = fineCircle.evaluateOnNearestBoundary(iF,jF,NORTH,(*boundary));
								hN = fineCircle.getSpacing(iF,jF)(NORTH);
							}
							else{
								double hNW, hNE;
								double uNW, uNE;
								if(fineCircle.getType(iF-1,jF+1) == EXTERIOR){
									uNW = fineCircle.evaluateOnNearestBoundary(iF,jF+1,WEST,(*boundary));
									hNW = fineCircle.getSpacing(iF,jF+1)(WEST);
								}
								else{
									uNW = uC((iF-1)/2,(jF+1)/2);
									hNW = h;
								}
								if(fineCircle.getType(iF+1,jF+1) == EXTERIOR){
									uNE = fineCircle.evaluateOnNearestBoundary(iF,jF+1,EAST,(*boundary));
									hNE = fineCircle.getSpacing(iF,jF+1)(EAST);
								}
								else{
									uNE = uC((iF+1)/2,(jF+1)/2);
									hNE = h;
								}
								uN = uNE*hNW/(hNE + hNW) + uNW*hNE/(hNE + hNW);
								hN = h;
							}

							// Get uW
							if(fineCircle.getType(iF-1,jF) == EXTERIOR){
								uW = fineCircle.evaluateOnNearestBoundary(iF,jF,WEST,(*boundary));
								hW = fineCircle.getSpacing(iF,jF)(WEST);
							}
							else{
								double hWN, hWS;
								double uWN, uWS;
								if(fineCircle.getType(iF-1,jF-1) == EXTERIOR){
									uWS = fineCircle.evaluateOnNearestBoundary(iF-1,jF,SOUTH,(*boundary));
									hWS = fineCircle.getSpacing(iF-1,jF)(SOUTH);
								}
								else{
									uWS = uC((iF-1)/2,(jF-1)/2);
									hWS = h;
								}
								if(fineCircle.getType(iF-1,jF+1) == EXTERIOR){
									uWN = fineCircle.evaluateOnNearestBoundary(iF-1,jF,NORTH,(*boundary));
									hWN = fineCircle.getSpacing(iF-1,jF)(NORTH);
								}
								else{
									uWN = uC((iF-1)/2,(jF+1)/2);
									hWN = h;
								}
								uW = uWS*hWN/(hWN + hWS) + uWN*hWS/(hWN + hWS);
								hW = h;
							}

							// Get uE
							if(fineCircle.getType(iF+1,jF) == EXTERIOR){
								uE = fineCircle.evaluateOnNearestBoundary(iF,jF,EAST,(*boundary));
								hE = fineCircle.getSpacing(iF,jF)(EAST);
							}
							else{
								double hEN, hES;
								double uEN, uES;
								if(fineCircle.getType(iF+1,jF-1) == EXTERIOR){
									uES = fineCircle.evaluateOnNearestBoundary(iF+1,jF,SOUTH,(*boundary));
									hES = fineCircle.getSpacing(iF+1,jF)(SOUTH);
								}
								else{
									uES = uC((iF+1)/2,(jF-1)/2);
									hES = h;
								}
								if(fineCircle.getType(iF+1,jF+1) == EXTERIOR){
									uEN = fineCircle.evaluateOnNearestBoundary(iF+1,jF,NORTH,(*boundary));
									hEN = fineCircle.getSpacing(iF+1,jF)(NORTH);
								}
								else{
									uEN = uC((iF+1)/2,(jF+1)/2);
									hEN = h;
								}
								uE = uES*hEN/(hEN + hES) + uEN*hES/(hEN + hES);
								hE = h;
							}
							uXY = uS*hN/(hN + hS) + uN*hS/(hN + hS);
							uYX = uW*hE/(hE + hW) + uE*hW/(hE + hW);

							uF(iF,jF) = 0.5*(uXY + uYX);
						}
					}
				}
			}

			return uF;
		}

		/* Performs lexicographic Gauss-Seidel point smoothing with
		 * initial solution u0, right hand side f, Dirichlet boundary
		 * data given by the function (*boundaryFunction) on
		 * Circle circ for N iterations.
		 */
		GridScalar GSLex(GridScalar u0, GridScalar f, double (*boundaryFunction)(double,double), Circle circ, int N){
			double h = circ.getH();
			GridScalar u = circ.makeScalar();
			u = u0;
			for(int n = 0; n < N; n++){
				for(int i = u.lbound(0); i <= u.ubound(0); i++){
					for(int j = u.lbound(1); j <= u.ubound(1); j++){
						if(circ.getType(i,j) != EXTERIOR){
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
									Coord b = circ.nearestBoundary(i,j,NORTH);
									uN = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uN = u(i,j+1);
								}
								if(circ.getType(i,j-1) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,SOUTH);
									uS = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uS = u(i,j-1);
								}
								if(circ.getType(i+1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,EAST);
									uE = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uE = u(i+1,j);
								}
								if(circ.getType(i-1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,WEST);
									uW = (*boundaryFunction)(b(0),b(1));
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

		/* Performs red-black Gauss-Seidel point smoothing with
		 * initial solution u0, right hand side f, Dirichlet boundary
		 * data given by the function (*boundaryFunction) on
		 * Circle circ for N iterations.
		 */
		GridScalar GSRB(GridScalar u0, GridScalar f, double (*boundaryFunction)(double,double), Circle circ, int N){
			double h = circ.getH();
			GridScalar u = circ.makeScalar();
			u = u0;
			for(int n = 0; n < N; n++){
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
									Coord b = circ.nearestBoundary(i,j,NORTH);
									uN = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uN = u(i,j+1);
								}
								if(circ.getType(i,j-1) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,SOUTH);
									uS = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uS = u(i,j-1);
								}
								if(circ.getType(i+1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,EAST);
									uE = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uE = u(i+1,j);
								}
								if(circ.getType(i-1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,WEST);
									uW = (*boundaryFunction)(b(0),b(1));
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
									Coord b = circ.nearestBoundary(i,j,NORTH);
									uN = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uN = u(i,j+1);
								}
								if(circ.getType(i,j-1) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,SOUTH);
									uS = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uS = u(i,j-1);
								}
								if(circ.getType(i+1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,EAST);
									uE = (*boundaryFunction)(b(0),b(1));
								}
								else{
									uE = u(i+1,j);
								}
								if(circ.getType(i-1,j) == EXTERIOR){
									Coord b = circ.nearestBoundary(i,j,WEST);
									uW = (*boundaryFunction)(b(0),b(1));
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

		/* Half weighting restriction operator restricts a GridScalar uF
		 * on a finely meshed circle with spacing h to a coarser circle
		 * with spacing 2*h.
		 */
		GridScalar HalfWeighting(GridScalar uF, Circle fine, Circle coarse, double (*boundaryFunction)(double,double)){
			GridScalar uC = coarse.makeScalar();
			for(int iC = uC.lbound(0); iC <= uC.ubound(0); iC++){
				for(int jC = uC.lbound(1); jC <= uC.ubound(1); jC++){
					int iF = 2*iC, jF = 2*jC;
					if(fine.getType(iF,jF) == REGULAR){
						int iF = 2*iC, jF = 2*jC;
						uC(iC,jC) = 4.0*uF(iF,jF)
								+ uF(iF+1,jF) + uF(iF-1,jF) + uF(iF,jF+1) + uF(iF,jF-1);
						uC(iC,jC) = uC(iC,jC)/8.0;
					}
					else if(fine.getType(iF,jF) == IRREGULAR){
						double uN = 0, uS = 0, uE = 0, uW = 0;
						double hN, hS, hE, hW;
						Spacing sp = fine.getSpacing(iF,jF);
						hN = sp(NORTH);
						hS = sp(SOUTH);
						hE = sp(EAST);
						hW = sp(WEST);
						if(fine.getType(iF,jF+1) == EXTERIOR){
							Coord b = fine.nearestBoundary(iF,jF,NORTH);
							uN = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uN = uF(iF,jF+1);
						}
						if(fine.getType(iF,jF-1) == EXTERIOR){
							Coord b = fine.nearestBoundary(iF,jF,SOUTH);
							uS = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uS = uF(iF,jF-1);
						}
						if(fine.getType(iF+1,jF) == EXTERIOR){
							Coord b = fine.nearestBoundary(iF,jF,EAST);
							uE = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uE = uF(iF+1,jF);
						}
						if(fine.getType(iF-1,jF) == EXTERIOR){
							Coord b = fine.nearestBoundary(iF,jF,WEST);
							uW = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uW = uF(iF-1,jF);
						}
						uC(iC,jC) = uF(iF,jF)
								  + (hN*uN + hS*uS)/(hN+hS)
								  + (hE*uE + hW*uW)/(hE+hW);
						uC(iC,jC) = uC(iC,jC)/2.0;
					}
				}
			}
			return uC;
		}

		/* Applies the 2D discrete Laplacian with Dirichlet boundary values
		 * given by the function pointed to by boundaryFunction to GridScalar
		 * u.
		 */
		GridScalar LaplaceOperator(GridScalar u, Circle circ, double (*boundaryFunction)(double,double)){
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
						double uN, uS, uE, uW;
						if(circ.getType(i,j+1) == EXTERIOR){
							Coord b = circ.nearestBoundary(i,j,NORTH);
							uN = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uN = u(i,j+1);
						}
						if(circ.getType(i,j-1) == EXTERIOR){
							Coord b = circ.nearestBoundary(i,j,SOUTH);
							uS = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uS = u(i,j-1);
						}
						if(circ.getType(i+1,j) == EXTERIOR){
							Coord b = circ.nearestBoundary(i,j,EAST);
							uE = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uE = u(i+1,j);
						}
						if(circ.getType(i-1,j) == EXTERIOR){
							Coord b = circ.nearestBoundary(i,j,WEST);
							uW = (*boundaryFunction)(b(0),b(1));
						}
						else{
							uW = u(i-1,j);
						}
						Lu(i,j) = 2.0*((uN/hN + uS/hS)/(hN + hS) + (uE/hE + uW/hW)/(hE + hW) - u(i,j)/(hW*hE) - u(i,j)/(hN*hS));
					}
				}
			}
			return Lu;
		}

		/* Performs one V-Cycle of a multigrid solver with finest grid fineCircle,
		 * v1 presmooths, v2 postsmooths, initial guess u0, right hand side f,
		 * and Dirichlet boundary data function (*boundary).
		 */
		GridScalar PoissonVCycle(GridScalar u0, GridScalar f, Circle fineCircle, double (*boundary)(double,double), int v1, int v2){
			double hF = fineCircle.getH();
			double hC = 2.0*hF;
			Circle coarseCircle(fineCircle.getR(),hC,fineCircle.getCenter());
			GridScalar u = fineCircle.makeScalar();
			u = GSLex(u0, f, (*boundary), fineCircle, v1);
			GridScalar Lu = fineCircle.makeScalar();
			Lu = LaplaceOperator(u,fineCircle,(*boundary));
			GridScalar rF = fineCircle.makeScalar();
			rF = f - Lu;
			GridScalar rC = coarseCircle.makeScalar();
			rC = HalfWeighting(rF, fineCircle, coarseCircle, (*boundary));
			GridScalar eC = coarseCircle.makeScalar();
			if(coarseCircle.interiorPoints() < 5){
				eC = GSLex(eC, rC, (*zeroBoundary), coarseCircle, 3);
			}
			else{
				eC = PoissonVCycle(eC, rC, coarseCircle, (*zeroBoundary), v1, v2);
			}
			GridScalar eF = fineCircle.makeScalar();
			eF = BilinearInterpolation(eC, fineCircle, coarseCircle, (*boundary));
			GridScalar uCorrected = fineCircle.makeScalar();
			uCorrected = u + eF;
			u = GSLex(uCorrected, f, (*boundary), fineCircle, v2);
			return u;
		}
	}
}

#endif /* POISSONMULTIGRID_H_ */
