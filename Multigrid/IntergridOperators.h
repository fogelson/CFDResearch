/*
 * IntergridOperators.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef INTERGRIDOPERATORS_H_
#define INTERGRIDOPERATORS_H_

#include "../Geometry/Geometry.h"

namespace CFD{
	using namespace Geometry;
	namespace Multigrid{
		class Interpolator;
		class Restrictor;


		class Interpolator{
		public:
			virtual ~Interpolator(){};
			// doInterpolate is really clumsy wording, but I want to have parallel
			// word choice with doRestrict, and restrict is a reserved word so it
			// can't just be interpolate() and restrict().
			virtual CellDoubleArray doInterpolate(CellDoubleArray &uC, Grid *fine, Grid *coarse){
				CellDoubleArray uF = fine->makeCellDoubleArray();
				doInterpolate(uF,uC,fine,coarse);
				return uF;
			}
			virtual void doInterpolate(CellDoubleArray & uF, CellDoubleArray & uC, Grid * fine, Grid * coarse) = 0;
		};

		class InjectiveInterpolator : public Interpolator{
		public:
			void doInterpolate(CellDoubleArray &uF, CellDoubleArray &uC, Grid *fine, Grid *coarse){
				//int iFMin = fine->iMin, iFMax = fine->iMax, jFMin = fine->jMin, jFMax = fine->jMax;
				int iCMin = coarse->iMin, iCMax = coarse->iMax, jCMin = coarse->jMin, jCMax = coarse->jMax;
				for(int iC = iCMin; iC <= iCMax; iC++){
					for(int jC = jCMin; jC <= jCMax; jC++){
						if(coarse->getCellType(iC,jC) != COVERED){
							uF(2*iC,2*jC) = uC(iC,jC);
							uF(2*iC,2*jC-1) = uC(iC,jC);
							uF(2*iC-1,2*jC) = uC(iC,jC);
							uF(2*iC-1,2*jC-1) = uC(iC,jC);
						}
					}
				}
			}
		};

		class BilinearInterpolator : public Interpolator{
			double regularStencil(double near, double middle1, double middle2, double far){
				return (9*near + 3*middle1 + 3*middle2 + far)/16;
			}
			double bilinearExtrapolation(double a, double b, double c, double d){
				return (15*a + 5*b - 3*c - d)/16;
			}
			double linearExtrapolation(double p, double q){
				return (5*p - q)/4;
			}
		public:
			void doInterpolate(CellDoubleArray &uF, CellDoubleArray &uC, Grid *fine, Grid *coarse){
				for(int iC = 1; iC <= coarse->iMax; iC++){
					for(int jC = 1; jC <= coarse->jMax; jC++){
						int iF = iC*2, jF = jC*2;
						if(coarse->isUncovered(iC,jC)){
							// NE
							if(fine->isRegular(iF,jF)){
								uF(iF,jF) = regularStencil(uC(iC,jC),uC(iC,jC+1),uC(iC+1,jC),uC(iC+1,jC+1));
							}
							else if(fine->isIrregular(iF,jF)){
								if(coarse->isUncovered(iC,jC+1)
										&& coarse->isUncovered(iC+1,jC)
										&& coarse->isUncovered(iC+1,jC+1)){
									uF(iF,jF) = regularStencil(uC(iC,jC),uC(iC,jC+1),uC(iC+1,jC),uC(iC+1,jC+1));
								}
								else if(coarse->isUncovered(iC,jC-1)
										&& coarse->isUncovered(iC+1,jC)
										&& coarse->isUncovered(iC+1,jC-1)){
									uF(iF,jF) = bilinearExtrapolation(uC(iC,jC),uC(iC+1,jC),uC(iC,jC-1),uC(iC+1,jC-1));
								}
								else if(coarse->isUncovered(iC-1,jC)
										&& coarse->isUncovered(iC,jC+1)
										&& coarse->isUncovered(iC-1,jC+1)){
									uF(iF,jF) = bilinearExtrapolation(uC(iC,jC),uC(iC,jC+1),uC(iC-1,jC),uC(iC-1,jC+1));
								}
								else if(coarse->isUncovered(iC-1,jC-1)){
									uF(iF,jF) = linearExtrapolation(uC(iC,jC),uC(iC-1,jC-1));
								}
								else{
									uF(iF,jF) = uC(iC,jC);
								}
							}
							// NW
							if((fine->isRegular(iF-1,jF))){
								uF(iF-1,jF) = regularStencil(uC(iC,jC),uC(iC-1,jC),uC(iC,jC+1),uC(iC-1,jC+1));
							}
							else if(fine->isIrregular(iF-1,jF)){
								if(coarse->isUncovered(iC-1,jC) // Look at NW block of coarse cells
										&& coarse->isUncovered(iC,jC+1)
										&& coarse->isUncovered(iC-1,jC+1)){
									uF(iF-1,jF) = regularStencil(uC(iC,jC),uC(iC-1,jC),uC(iC,jC+1),uC(iC-1,jC+1));
								}
								else if(coarse->isUncovered(iC,jC-1) // Look at SW block of coarse cells
										&& coarse->isUncovered(iC-1,jC)
										&& coarse->isUncovered(iC-1,jC-1)){
									uF(iF-1,jF) = bilinearExtrapolation(uC(iC,jC),uC(iC-1,jC),uC(iC,jC-1),uC(iC-1,jC-1));
								}
								else if(coarse->isUncovered(iC,jC+1) // Look at NE block of coarse cells
										&& coarse->isUncovered(iC+1,jC)
										&& coarse->isUncovered(iC+1,jC+1)){
									uF(iF-1,jF) = bilinearExtrapolation(uC(iC,jC),uC(iC,jC+1),uC(iC+1,jC),uC(iC+1,jC+1));
								}
								else if(coarse->isUncovered(iC-1,jC-1)){ // Look at diagonal SE coarse cell
									uF(iF-1,jF) = linearExtrapolation(uC(iC,jC),uC(iC-1,jC-1));
								}
								else{ // Look at current coarse cell
									uF(iF-1,jF) = uC(iC,jC);
								}
							}
							// SE
							if(fine->getCellType(iF,jF-1) == REGULAR){
								uF(iF,jF-1) = regularStencil(uC(iC,jC),uC(iC,jC-1),uC(iC+1,jC),uC(iC+1,jC-1));
							}
							else if(fine->getCellType(iF,jF-1) == IRREGULAR){
								if(coarse->isUncovered(iC+1,jC) // Look at SE block of coarse cells
										&& coarse->isUncovered(iC,jC-1)
										&& coarse->isUncovered(iC+1,jC-1)){
									uF(iF,jF-1) = regularStencil(uC(iC,jC),uC(iC,jC-1),uC(iC+1,jC),uC(iC+1,jC-1));
								}
								else if(coarse->isUncovered(iC+1,jC) // Look at NE block of coarse cells
										&& coarse->isUncovered(iC,jC+1)
										&& coarse->isUncovered(iC+1,jC+1)){
									uF(iF,jF-1) = bilinearExtrapolation(uC(iC,jC),uC(iC+1,jC),uC(iC,jC+1),uC(iC+1,jC+1));
								}
								else if(coarse->isUncovered(iC,jC-1) // Look at SW block of coarse cells
										&& coarse->isUncovered(iC-1,jC)
										&& coarse->isUncovered(iC-1,jC-1)){
									uF(iF,jF-1) = bilinearExtrapolation(uC(iC,jC),uC(iC,jC-1),uC(iC-1,jC),uC(iC-1,jC-1));
								}
								else if(coarse->isUncovered(iC-1,jC+1)){ // Look at diagonal NW coarse cell
									uF(iF,jF-1) = linearExtrapolation(uC(iC,jC),uC(iC-1,jC+1));
								}
								else{
									uF(iF,jF-1) = uC(iC,jC);
								}
							}
							// SW
							if(fine->getCellType(iF-1,jF-1) == REGULAR){
								uF(iF-1,jF-1) = regularStencil(uC(iC,jC),uC(iC-1,jC),uC(iC,jC-1),uC(iC-1,jC-1));
							}
							else if(fine->getCellType(iF-1,jF-1) == IRREGULAR){
								if(coarse->isUncovered(iC,jC-1) // Look at SW block of coarse cells
										&& coarse->isUncovered(iC-1,jC)
										&& coarse->isUncovered(iC-1,jC-1)){
									uF(iF-1,jF-1) = regularStencil(uC(iC,jC),uC(iC,jC-1),uC(iC-1,jC),uC(iC-1,jC-1));
								}
								else if(coarse->isUncovered(iC-1,jC) // Look at NW block of coarse cells
										&& coarse->isUncovered(iC,jC+1)
										&& coarse->isUncovered(iC-1,jC+1)){
									uF(iF-1,jF-1) = bilinearExtrapolation(uC(iC,jC),uC(iC-1,jC),uC(iC,jC+1),uC(iC-1,jC+1));
								}
								else if(coarse->isUncovered(iC+1,jC) // Look at SE block of coarse cells
										&& coarse->isUncovered(iC,jC-1)
										&& coarse->isUncovered(iC+1,jC-1)){
									uF(iF-1,jF-1) = bilinearExtrapolation(uC(iC,jC),uC(iC,jC-1),uC(iC+1,jC),uC(iC+1,jC-1));
								}
								else if(coarse->isUncovered(iC+1,jC+1)){ // Look at diagonal NE coarse cell
									uF(iF-1,jF-1) = linearExtrapolation(uC(iC,jC),uC(iC+1,jC+1));
								}
								else{
									uF(iF-1,jF-1) = uC(iC,jC);
								}
							}
						}
					}
				}
				/*for(int iF = 1; iF <= fine->iMax; iF++){
					for(int jF = 1; jF <= fine->jMax; jF++){
						if(fine->getCellType(iF,jF) == REGULAR){
							uF(iF,jF) = 1;
						}
					}
				}*/
				CellDoubleArray out(fine->xRange,fine->yRange);
				out = 1.0*(fine->cellTypes);
				return uF;
			}
		};



		class BilinearInterpolatorDumbButtface : public Interpolator{
			double bilinearExtrapolation(double a, double b, double c, double d){
				/*			   e
				 * 			   |
				 * 			a------------b
				 * 			|  |		 |
				 *  		|  |		 |
				 *   		|  |		 |
				 *    		|  |		 |
				 *    		c------------d
				 *
				 */

				return (15*a + 5*b - 3*c - d)/16;
			}

			double linearExtrapolation(double p, double q){
				/*
				 * 				*
				 * 			   /
				 * 		0-----p
				 *      |	 /|
				 *      |  /  |
				 *      |/	  |
				 *      q-----0
				 */

				return (5*p - q)/4;
			}

			Direction getQuadrant(int iF, int jF){
				Direction out;
				if((iF % 2) == 0){ // E
					if((jF % 2) == 0){ // N
						out = NE;
					}
					else{ // S
						out = SE;
					}
				}
				else{ // W
					if((jF % 2) == 0){ // N
						out = NW;
					}
					else{
						out = SW;
					}
				}
				return out;
			}
			TinyVector<int,2> getCoarseCellIndices(int iF, int jF){
				TinyVector<int,2> out;
				out(0) = (iF + (iF % 2))/2;
				out(1) = (jF + (jF % 2))/2;
				//cout << "(" << iF << ", " << jF << ") ---> (" << out(0) << ", " << out(1) << ")" << endl;
				return out;
			}
			TinyVector<Type,8> getNeighboringCoarseCellTypes(int iF, int jF, Grid *coarse){
				TinyVector<int,2> coarseCellIndices = getCoarseCellIndices(iF,jF);
				int iC = coarseCellIndices(0), jC = coarseCellIndices(1);

				TinyVector<Type,8> out;

				out(N) = coarse->getCellType(iC,jC+1);
				out(S) = coarse->getCellType(iC,jC-1);
				out(E) = coarse->getCellType(iC+1,jC);
				out(W) = coarse->getCellType(iC-1,jC);
				out(NW) = coarse->getCellType(iC-1,jC+1);
				out(NE) = coarse->getCellType(iC+1,jC+1);
				out(SW) = coarse->getCellType(iC-1,jC-1);
				out(SE) = coarse->getCellType(iC+1,jC-1);

				return out;
			}
//			void applyIrregularStencil(CellDoubleArray &uC, CellDoubleArray &uF, int iC, int jC, int iF, int jF, Grid *coarse){
//				uF(iF,jF) = -10;
//				return;
//				Direction quadrant = getQuadrant(iF,jF);
//				TinyVector<bool,CARDINAL_DIRECTIONS> uncovered;
//				bool applyRegular = false, applyBilinear = false, applyLinear = false;
//				uncovered = true;
//
//				double a, b, c, d;
//				double p, q;
//
//				if(coarse->getCellType(iC,jC+1) == COVERED){
//					uncovered(N) = false;
//				}
//				if(coarse->getCellType(iC+1,jC) == COVERED){
//					uncovered(E) = false;
//				}
//				if(coarse->getCellType(iC,jC-1) == COVERED){
//					uncovered(S) = false;
//				}
//				if(coarse->getCellType(iC-1,jC) == COVERED){
//					uncovered(W) = false;
//				}
//				if(coarse->getCellType(iC+1,jC+1) == COVERED){
//					uncovered(NE) = false;
//				}
//				if(coarse->getCellType(iC+1,jC-1) == COVERED){
//					uncovered(SE) = false;
//				}
//				if(coarse->getCellType(iC-1,jC+1) == COVERED){
//					uncovered(NW) = false;
//				}
//				if(coarse->getCellType(iC-1,jC-1) == COVERED){
//					uncovered(SW) = false;
//				}
//
//				if(quadrant == NE){
//					if(uncovered(N) && uncovered(NE) && uncovered(E)){
//						applyRegular = true;
//					}
//					else if(uncovered(N) && uncovered(NW) && uncovered(W)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC,jC+1);
//						c = uC(iC-1,jC);
//						d = uC(iC-1,jC+1);
//					}
//					else if(uncovered(S) && uncovered(SE) && uncovered(E)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC+1,jC);
//						c = uC(iC,jC-1);
//						d = uC(iC+1,jC-1);
//					}
//					else if(uncovered(SW)){
//						applyLinear = true;
//						p = uC(iC,jC);
//						q = uC(iC-1,jC-1);
//					}
//				}
//				else if(quadrant == NW){
//					if(uncovered(N) && uncovered(NW) && uncovered(W)){
//						applyRegular = true;
//					}
//					else if(uncovered(N) && uncovered(NE) && uncovered(E)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC,jC+1);
//						c = uC(iC+1,jC);
//						d = uC(iC+1,jC+1);
//					}
//					else if(uncovered(S) && uncovered(SW) && uncovered(W)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC-1,jC);
//						c = uC(iC,jC-1);
//						d = uC(iC-1,jC-1);
//					}
//					else if(uncovered(SE)){
//						applyLinear = true;
//						p = uC(iC,jC);
//						q = uC(iC+1,jC-1);
//					}
//				}
//				else if(quadrant == SE){
//					if(uncovered(S) && uncovered(SE) && uncovered(E)){
//						applyRegular = true;
//					}
//					else if(uncovered(N) && uncovered(NE) && uncovered(E)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC+1,jC);
//						c = uC(iC,jC+1);
//						d = uC(iC+1,jC+1);
//					}
//					else if(uncovered(S) && uncovered(SW) && uncovered(W)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC,jC-1);
//						c = uC(iC-1,jC);
//						d = uC(iC-1,jC-1);
//					}
//					else if(uncovered(NW)){
//						applyLinear = true;
//						p = uC(iC,jC);
//						q = uC(iC-1,jC+1);
//					}
//				}
//				else if(quadrant == SW){
//					if(uncovered(S) && uncovered(SW) && uncovered(W)){
//						applyRegular = true;
//					}
//					else if(uncovered(N) && uncovered(NW) && uncovered(W)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC-1,jC);
//						c = uC(iC,jC+1);
//						d = uC(iC-1,jC+1);
//					}
//					else if(uncovered(S) && uncovered(SE) && uncovered(E)){
//						applyBilinear = true;
//						a = uC(iC,jC);
//						b = uC(iC,jC-1);
//						c = uC(iC+1,jC);
//						d = uC(iC+1,jC-1);
//					}
//					else if(uncovered(NE)){
//						applyLinear = true;
//						p = uC(iC,jC);
//						q = uC(iC+1,jC+1);
//					}
//				}
//				if(applyRegular){
//					applyRegularStencil(uC, uF, iC, jC, iF, jF);
//					return;
//				}
//				else if(applyBilinear){
//					uF(iF,jF) = bilinearExtrapolation(a,b,c,d);
//					return;
//				}
//				else if(applyLinear){
//					uF(iF,jF) = linearExtrapolation(p,q);
//					return;
//				}
//				else{
//					uF(iF,jF) = -10;
//					cout << "Uh oh. Something went wrong interpolating near the boundary." << endl;
//					cout << "iF = " << iF << ", jF = " << jF << ", iC = " << iC << ", jC = " << jC << endl;
//				}
//			}
//			void applyRegularStencil(CellDoubleArray &uC, CellDoubleArray &uF, int iC, int jC, int iF, int jF){
//				Direction quadrant = getQuadrant(iF,jF);
//				if(quadrant == NE){
//					uF(iF,jF) = 9*uC(iC,jC)
//							+ 3*uC(iC+1,jC)
//							+ 3*uC(iC,jC+1)
//							+ uC(iC+1,jC+1);
//					uF(iF,jF) = uF(iF,jF)/16;
//				}
//				else if(quadrant == NW) {
//					uF(iF,jF) = 9*uC(iC,jC)
//							+ 3*uC(iC-1,jC)
//							+ 3*uC(iC,jC+1)
//							+ uC(iC-1,jC+1);
//					uF(iF,jF) = uF(iF,jF)/16;
//				}
//				else if(quadrant == SE){
//					uF(iF,jF) = 9*uC(iC,jC)
//							+ 3*uC(iC+1,jC)
//							+ 3*uC(iC,jC-1)
//							+ uC(iC+1,jC-1);
//					uF(iF,jF) = uF(iF,jF)/16;
//				}
//				else if(quadrant == SW){
//					uF(iF,jF) = 9*uC(iC,jC)
//							+ 3*uC(iC-1,jC)
//							+ 3*uC(iC,jC-1)
//							+ uC(iC-1,jC-1);
//					uF(iF,jF) = uF(iF,jF)/16;
//				}
//			}
		public:
			void doInterpolate(CellDoubleArray &uF, CellDoubleArray &uC, Grid *fine, Grid *coarse){
				//uC = coarse->makeCellDoubleArray();
				//uC = 1;
				cout << uC.shape() << endl;
				cout << uF.shape() << endl;
				int iFMin = fine->iMin, iFMax = fine->iMax, jFMin = fine->jMin, jFMax = fine->jMax;
				//int iCMin = coarse->iMin, iCMax = coarse->iMax, jCMin = coarse->jMin, jCMax = coarse->jMax;

				for(int iF = iFMin; iF <= iFMax; iF++){
					for(int jF = jFMin; jF <= jFMax; jF++){
						if(fine->getCellType(iF,jF) == REGULAR){
							TinyVector<int,2> coarseIndices = getCoarseCellIndices(iF,jF);
							int iC = coarseIndices(0), jC = coarseIndices(1);

							TinyVector<double,8> neighboringCoarseCellTypes = getNeighboringCoarseCellTypes(iF,jF,coarse);

							Direction fineQuadrant = getQuadrant(iF,jF);

							if(fineQuadrant == NE){
								if(neighboringCoarseCellTypes(N) != COVERED && neighboringCoarseCellTypes(NE) != COVERED && neighboringCoarseCellTypes(E) != COVERED){
									uF(iF,jF) = 9*uC(iC,jC)
										+ 3*uC(iC+1,jC)
										+ 3*uC(iC,jC+1)
										+ uC(iC+1,jC+1);
									uF(iF,jF) = uF(iF,jF)/16;
									//uF(iF,jF) = uC(iC,jC);
								}
							}
						}
						else{
							uF(iF,jF) = 0;
						}
					}
				}
				//cout << uF << endl;

/*				for(int iC = iCMin; iC <= iCMax; iC++){
					for(int jC = jCMin; jC <= jCMax; jC++){
						uF(2*iC,2*jC) += 9*uC(iC,jC);
						uF(2*iC,2*jC-1) += 9*uC(iC,jC);
						uF(2*iC-1,2*jC) += 9*uC(iC,jC);
						uF(2*iC-1,2*jC-1) += 9*uC(iC,jC);


						uF(2*iC,2*jC+1) += 3*uC(iC,jC);
						uF(2*iC-1,2*jC+1) += 3*uC(iC,jC);

						uF(2*iC+1,2*jC) += 3*uC(iC,jC);
						uF(2*iC+1,2*jC-1) += 3*uC(iC,jC);

						uF(2*iC,2*jC-2) += 3*uC(iC,jC);
						uF(2*iC-1,2*jC-2) += 3*uC(iC,jC);

						uF(2*iC-2,2*jC) += 3*uC(iC,jC);
						uF(2*iC-2,2*jC-1) += 3*uC(iC,jC);


						uF(2*iC+1,2*jC+1) += uC(iC,jC);
						uF(2*iC-2,2*jC+1) += uC(iC,jC);
						uF(2*iC+1,2*jC-2) += uC(iC,jC);
						uF(2*iC-2,2*jC-2) += uC(iC,jC);
					}
				}*/
				return uF;
			}
		};


		class Restrictor{
		public:
			virtual ~Restrictor(){};
			// restrict is a reserved word in C++, so we can't just have this
			// function called restrict(). So we call it doRestrict() instead.
			virtual CellDoubleArray doRestrict(CellDoubleArray &uF, Grid *fine, Grid *coarse){
				CellDoubleArray uC = coarse->makeCellDoubleArray();
				doRestrict(uF,uC,fine,coarse);
				return uC;
			}
			virtual void doRestrict(CellDoubleArray & uC, CellDoubleArray & uF, Grid * fine, Grid * coarse) = 0;
		};

		class VolumeWeightedRestrictor : public Restrictor{
		public:
			void doRestrict(CellDoubleArray &uC, CellDoubleArray &uF, Grid *fine, Grid *coarse){
				//int iFMin = fine->iMin, iFMax = fine->iMax, jFMin = fine->jMin, jFMax = fine->jMax;
				int iCMin = coarse->iMin, iCMax = coarse->iMax, jCMin = coarse->jMin, jCMax = coarse->jMax;
				for(int iC = iCMin; iC <= iCMax; iC++){
					for(int jC = jCMin; jC <= jCMax; jC++){
						if(coarse->isUncovered(iC,jC)){
							double vNE = fine->volumeFractions(2*iC,2*jC);
							double vNW = fine->volumeFractions(2*iC-1,2*jC);
							double vSE = fine->volumeFractions(2*iC,2*jC-1);
							double vSW = fine->volumeFractions(2*iC-1,2*jC-1);
							double uNE = uF(2*iC,2*jC);
							double uNW = uF(2*iC-1,2*jC);
							double uSE = uF(2*iC,2*jC-1);
							double uSW = uF(2*iC-1,2*jC-1);

							uC(iC,jC) = (vNE*uNE + vNW*uNW + vSE*uSE + vSW*uSW)/(vNE + vNW + vSE + vSW);
						}
					}
				}
			}
		};
	}
}

#endif /* INTERGRIDOPERATORS_H_ */
