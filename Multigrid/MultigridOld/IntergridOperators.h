/*
 * IntergridOperators.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef INTERGRIDOPERATORS_H_
#define INTERGRIDOPERATORS_H_

#include "../Geometry.h"

using namespace CFD;
using namespace Geometry;

namespace CFD{
	namespace Multigrid{
		class Interpolator;
		class Restrictor;


		class Interpolator{
		public:
			virtual ~Interpolator(){};
			// doInterpolate is really clumsy wording, but I want to have parallel
			// word choice with doRestrict, and restrict is a reserved word so it
			// can't just be interpolate() and restrict().
			virtual GridScalar doInterpolate(GridScalar uC, Circle fine, Circle coarse) = 0;
		};
		/*class NewBilinearInterpolator : public Interpolator{
			virtual ~BilinearInterpolator(){}
			virtual GridScalar doInterpolate(GridScalar uC, Circle fine, Circle coarse){
				GridScalar uF = fine.makeScalar();
			}

		};*/
		class BilinearInterpolatorExtrapolated : public Interpolator{
		public:
			virtual ~BilinearInterpolatorExtrapolated(){}
			virtual GridScalar doInterpolate(GridScalar uC, Circle fine, Circle coarse){
				//double h = fine.getH();
				GridScalar uF = fine.makeScalar();

				for(int iF = uF.lbound(0); iF <= uF.ubound(0); iF++){
					for(int jF = uF.lbound(1); jF <= uF.ubound(1); jF++){
						if(fine.getType(iF,jF) != EXTERIOR){
							if((iF % 2) == 0 && (jF % 2) == 0){
								uF(iF,jF) = uC(iF/2,jF/2);
							}
							else if((iF % 2) == 0){
								uF(iF,jF) = 0.5*(uC(iF/2,(jF+1)/2) + uC(iF/2,(jF-1)/2));
							}
							else if((jF % 2) == 0){
								uF(iF,jF) = 0.5*(uC((iF+1)/2,jF/2) + uC((iF-1)/2,jF/2));
							}
							else{
								uF(iF,jF) = 0.25*(uC((iF+1)/2,(jF+1)/2) + uC((iF+1)/2,(jF-1)/2) + uC((iF-1)/2,(jF+1)/2) + uC((iF-1)/2,(jF-1)/2));
							}
						}
					}
				}
				return uF;

				return uF;
			}
		};
		class BilinearInterpolator : public Interpolator{
		public:
			virtual ~BilinearInterpolator(){}
			virtual GridScalar doInterpolate(GridScalar uC, Circle fine, Circle coarse){
				double h = fine.getH();
				GridScalar uF = fine.makeScalar();

				/* Iterate over the fine grid. This is different than
				 * when we are on a normal rectangular grid, because here
				 * we have to treat the irregular points on the fine grid
				 * extremely carefully, so it makes a lot more sense to
				 * just iterate over that grid.
				 */
				for(int iF = uF.lbound(0); iF <= uF.ubound(0); iF++){
					for(int jF = uF.lbound(1); jF <= uF.ubound(1); jF++){
						if(fine.getType(iF,jF) != EXTERIOR){
							// If the fine point overlaps a coarse point, set equal
							if(iF % 2 == 0 && jF % 2 == 0){
								uF(iF,jF) = uC(iF/2,jF/2);
							}
							// Interpolate in just the y-direction for points along a coarse
							// grid line in the x-direction.
							else if(iF % 2 == 0){
								if(fine.getType(iF,jF) == REGULAR){
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
									hN = fine.getSpacing(iF,jF)(NORTH);
									hS = fine.getSpacing(iF,jF)(SOUTH);
									if(fine.getType(iF,jF+1) != EXTERIOR){
										UN = uC(iF/2,(jF+1)/2);
									}
									else{
										UN = fine.boundaryValue(iF,jF,NORTH);
									}
									if(fine.getType(iF,jF-1) != EXTERIOR){
										US = uC(iF/2,(jF-1)/2);
									}
									else{
										US = fine.boundaryValue(iF,jF,SOUTH);
									}
									uF(iF,jF) = (hS*UN + hN*US)/(hS + hN);
								}
							}
							// Interpolate in just the x-direction for points along a coarse
							// grid line in the y-direction. See the above code for the
							// analogous case if you want detailed comments.
							else if(jF % 2 == 0){
								if(fine.getType(iF,jF) == REGULAR){
									uF(iF,jF) = 0.5*(uC((iF+1)/2,jF/2) + uC((iF-1)/2,jF/2));
								}
								else{
									double UE, UW, hE, hW;
									hE = fine.getSpacing(iF,jF)(EAST);
									hW = fine.getSpacing(iF,jF)(WEST);
									if(fine.getType(iF+1,jF) != EXTERIOR){
										UE = uC((iF+1)/2,jF/2);
									}
									else{
										UE = fine.boundaryValue(iF,jF,EAST);
									}
									if(fine.getType(iF-1,jF) != EXTERIOR){
										UW = uC((iF-1)/2,jF/2);
									}
									else{
										UW = fine.boundaryValue(iF,jF,WEST);
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
								if(fine.getType(iF,jF-1) == EXTERIOR){
									uS = fine.boundaryValue(iF,jF,SOUTH);
									hS = fine.getSpacing(iF,jF)(SOUTH);
								}
								// Interpolate between the SW and SE points, using boundary
								// points as needed.
								else{
									double hSW, hSE;
									double uSW, uSE;
									if(fine.getType(iF-1,jF-1) == EXTERIOR){
										uSW = fine.boundaryValue(iF,jF-1,WEST);
										hSW = fine.getSpacing(iF,jF-1)(WEST);
									}
									else{
										uSW = uC((iF-1)/2,(jF-1)/2);
										hSW = h;
									}
									if(fine.getType(iF+1,jF-1) == EXTERIOR){
										uSE = fine.boundaryValue(iF,jF-1,EAST);
										hSE = fine.getSpacing(iF,jF-1)(EAST);
									}
									else{
										uSE = uC((iF+1)/2,(jF-1)/2);
										hSE = h;
									}
									uS = uSE*hSW/(hSE + hSW) + uSW*hSE/(hSE + hSW);
									hS = h;
								}

								// Get uN
								if(fine.getType(iF,jF+1) == EXTERIOR){
									uN = fine.boundaryValue(iF,jF,NORTH);
									hN = fine.getSpacing(iF,jF)(NORTH);
								}
								else{
									double hNW, hNE;
									double uNW, uNE;
									if(fine.getType(iF-1,jF+1) == EXTERIOR){
										uNW = fine.boundaryValue(iF,jF+1,WEST);
										hNW = fine.getSpacing(iF,jF+1)(WEST);
									}
									else{
										uNW = uC((iF-1)/2,(jF+1)/2);
										hNW = h;
									}
									if(fine.getType(iF+1,jF+1) == EXTERIOR){
										uNE = fine.boundaryValue(iF,jF+1,EAST);
										hNE = fine.getSpacing(iF,jF+1)(EAST);
									}
									else{
										uNE = uC((iF+1)/2,(jF+1)/2);
										hNE = h;
									}
									uN = uNE*hNW/(hNE + hNW) + uNW*hNE/(hNE + hNW);
									hN = h;
								}

								// Get uW
								if(fine.getType(iF-1,jF) == EXTERIOR){
									uW = fine.boundaryValue(iF,jF,WEST);
									hW = fine.getSpacing(iF,jF)(WEST);
								}
								else{
									double hWN, hWS;
									double uWN, uWS;
									if(fine.getType(iF-1,jF-1) == EXTERIOR){
										uWS = fine.boundaryValue(iF-1,jF,SOUTH);
										hWS = fine.getSpacing(iF-1,jF)(SOUTH);
									}
									else{
										uWS = uC((iF-1)/2,(jF-1)/2);
										hWS = h;
									}
									if(fine.getType(iF-1,jF+1) == EXTERIOR){
										uWN = fine.boundaryValue(iF-1,jF,NORTH);
										hWN = fine.getSpacing(iF-1,jF)(NORTH);
									}
									else{
										uWN = uC((iF-1)/2,(jF+1)/2);
										hWN = h;
									}
									uW = uWS*hWN/(hWN + hWS) + uWN*hWS/(hWN + hWS);
									hW = h;
								}

								// Get uE
								if(fine.getType(iF+1,jF) == EXTERIOR){
									uE = fine.boundaryValue(iF,jF,EAST);
									hE = fine.getSpacing(iF,jF)(EAST);
								}
								else{
									double hEN, hES;
									double uEN, uES;
									if(fine.getType(iF+1,jF-1) == EXTERIOR){
										uES = fine.boundaryValue(iF+1,jF,SOUTH);
										hES = fine.getSpacing(iF+1,jF)(SOUTH);
									}
									else{
										uES = uC((iF+1)/2,(jF-1)/2);
										hES = h;
									}
									if(fine.getType(iF+1,jF+1) == EXTERIOR){
										uEN = fine.boundaryValue(iF+1,jF,NORTH);
										hEN = fine.getSpacing(iF+1,jF)(NORTH);
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
		};

		class Restrictor{
		public:
			virtual ~Restrictor(){};
			// restrict is a reserved word in C++, so we can't just have this
			// function called restrict(). So we call it doRestrict() instead.
			virtual GridScalar doRestrict(GridScalar uF, Circle fine, Circle coarse) = 0;
		};

		class HalfWeighter : public Restrictor{
		public:
			virtual ~HalfWeighter(){}
			virtual GridScalar doRestrict(GridScalar uF, Circle fine, Circle coarse){
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
							/*double uN = 0, uS = 0, uE = 0, uW = 0;
							double hN, hS, hE, hW;
							Spacing sp = fine.getSpacing(iF,jF);
							hN = sp(NORTH);
							hS = sp(SOUTH);
							hE = sp(EAST);
							hW = sp(WEST);
							if(fine.getType(iF,jF+1) == EXTERIOR){
								uN = fine.boundaryValue(iF,jF,NORTH);
							}
							else{
								uN = uF(iF,jF+1);
							}
							if(fine.getType(iF,jF-1) == EXTERIOR){
								uS = fine.boundaryValue(iF,jF,SOUTH);
							}
							else{
								uS = uF(iF,jF-1);
							}
							if(fine.getType(iF+1,jF) == EXTERIOR){
								uE = fine.boundaryValue(iF,jF,EAST);
							}
							else{
								uE = uF(iF+1,jF);
							}
							if(fine.getType(iF-1,jF) == EXTERIOR){
								uW = fine.boundaryValue(iF,jF,WEST);
							}
							else{
								uW = uF(iF-1,jF);
							}*/
							Spacing sp = fine.getSpacing(iF,jF);
							double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
							NeighborScalar nbr = fine.getNeighbors(uF,iF,jF);
							double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
							uC(iC,jC) = uF(iF,jF)
									  + (hN*uN + hS*uS)/(hN+hS)
									  + (hE*uE + hW*uW)/(hE+hW);
							uC(iC,jC) = uC(iC,jC)/2.0;
						}
					}
				}
				return uC;
			}
		};

		class FullWeighter : public Restrictor{
		public:
			virtual ~FullWeighter(){}
			virtual GridScalar doRestrict(GridScalar uF, Circle fine, Circle coarse){
				GridScalar uC = coarse.makeScalar();
				for(int iC = uC.lbound(0); iC <= uC.ubound(0); iC++){
					for(int jC = uC.lbound(1); jC <= uC.ubound(1); jC++){
						int iF = 2*iC, jF = 2*jC;
						double uN, uS, uE, uW, uNE, uSE, uNW, uSW;
						if(fine.getType(iF,jF) == REGULAR){
							uN = uF(iF,jF+1);
							uS = uF(iF,jF-1);
							uE = uF(iF+1,jF);
							uW = uF(iF-1,jF);
							uNE = uF(iF+1,jF+1);
							uNW = uF(iF-1,jF+1);
							uSE = uF(iF+1,jF-1);
							uSW = uF(iF-1,jF-1);
						}
						else if(fine.getType(iF,jF) == IRREGULAR){
							TinyVector<double,8> nbrExtrap = fine.extrapolateNeighbors(uF,iF,jF);
							uN = nbrExtrap(NORTH);
							uS = nbrExtrap(SOUTH);
							uE = nbrExtrap(EAST);
							uW = nbrExtrap(WEST);
							uNE = nbrExtrap(NORTHEAST);
							uSE = nbrExtrap(SOUTHEAST);
							uNW = nbrExtrap(NORTHWEST);
							uSW = nbrExtrap(SOUTHWEST);
						}
						uC(iC,jC) = (4.0*uF(iF,jF) + 2.0*(uN + uS + uE + uW) + (uNE + uSE + uNW + uSW))/16.0;
					}
				}
				return uC;
			}
		};
	}
}

#endif /* INTERGRIDOPERATORS_H_ */
