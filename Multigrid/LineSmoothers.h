/*
 * LineSmoothers.h
 *
 *  Created on: Mar 9, 2011
 *      Author: fogelson
 */

#ifndef LINESMOOTHERS_H_
#define LINESMOOTHERS_H_

#include <lapackpp/lapackpp.h>
#include "Smoothers.h"
//#include <vecLib/veclib.h>

using namespace la;

namespace CFD{
	namespace Multigrid{
		class PoissonXLineLex : public Smoother{
		public:
			virtual ~PoissonXLineLex(){}
			virtual GridScalar smooth(GridScalar u0, GridScalar f, Circle circ, int its){
				GridScalar u = circ.makeScalar();
				u = u0;
				double h = circ.getH();
				Range I(u.lbound(0),u.ubound(0));
				/*cout << "I: " << I << endl;
				cout << u.lbound(0) << endl;
				cout << u.ubound(0) << endl;*/
				if(circ.interiorPoints() == 1){
					//cout << "Got to single interior point level. " << endl; ;
					//cout << "\tu is " << u.length(0) << "x" << u.length(1) << ". " <<endl;
					//cout << "\th is " << circ.getH() << ". " << endl;
					//cout << "\tThere are " << circ.interiorPoints() << " interior points." << endl;
				}
				for(int n = 0; n < its; n++){
					if(circ.interiorPoints() == 1){
						for(int i = u.lbound(0); i <= u.ubound(0); i++){
							for(int j = u.lbound(1); j <= u.ubound(1); j++){
								//cout << "\tTried point " << i << ", " << j << ". It has type " << circ.getType(i,j) << "." << endl;
								if(circ.getType(i,j) != EXTERIOR){
									Spacing sp = circ.getSpacing(i,j);
									double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
									NeighborScalar nbr = circ.getNeighbors(u,i,j);
									double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
									u(i,j) = (uW/hW + uE/hE)/(hW + hE) + (uN/hN + uS/hS)/(hN + hS) - f(i,j)/2.0;
									u(i,j) = u(i,j)*(hN*hS*hE*hW)/(hN*hS + hE*hW);
									//cout << "\t\tThe interior point is at (" << circ.getX(i,j) << ", " << circ.getY(i,j) << "). " << endl;
									//cout << "\t\thN = " << hN << ", hS = " << hS << ", hE = " << hE << ", hW = " << hW << "." << endl;
								}
							}
						}
					}
					else{
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(abs(j) % 2 == 0){
								Array<double,1> lower, main, upper, right;
								lower.resize(I);
								main.resize(I);
								upper.resize(I);
								right.resize(I);
								//cout << "Right size: " << right.size() << endl;
								lower = 0;
								main = 1;
								upper = 0;
								right = 0;
								for(int i = u.lbound(0); i <= u.ubound(0); i++){
									if(circ.getType(i,j) == REGULAR){
										lower(i) = 1.0/pow2(h);
										main(i) = -4.0/pow2(h);
										upper(i) = 1.0/pow2(h);
										right(i) = f(i,j) - (u(i,j-1) + u(i,j+1))/pow2(h);
									}
									else if(circ.getType(i,j) == IRREGULAR){
										Spacing sp = circ.getSpacing(i,j);
										double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
										NeighborScalar nbr = circ.getNeighbors(u,i,j);
										double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
										lower(i) = (2.0/hW)/(hE + hW);
										main(i) = -(2.0/(hN*hS) + 2.0/(hE*hW));
										upper(i) = (2.0/hE)/(hE + hW);
										right(i) = f(i,j) - 2.0*(uN/hN + uS/hS)/(hN + hS);
										if(circ.getType(i-1,j) == EXTERIOR){
											/*cout << "***" << endl;
											cout << "At " << i << ", want to set " << i-1 << " to " << uW << endl;
											cout << "Current right value: " << right(i-1) << endl;
											cout << "In range: " << right.isInRange(i-1) << endl;
											cout << "***" << endl;*/
											right(i-1) = uW;
										}
										if(circ.getType(i+1,j) == EXTERIOR){
											/*cout << "***" << endl;
											cout << "At " << i << ", want to set " << i+1 << " to " << uE << endl;
											cout << "Current right value: " << right(i+1) << endl;
											cout << "In range: " << right.isInRange(i+1) << endl;
											cout << "***" << endl;*/
											right(i+1) = uE;
										}
									}
									//cout << right << endl;
								}


								//cout << main.size() << endl;
								//cout << "---" << endl;
								int K = main.size();
								int offset = main.lbound(0);

								LaVectorDouble laLower(K-1), laMain(K), laUpper(K-1), laRight(K), laU(K);
								laU = 0;
								for(int k = 0; k < K-1; k++){
									laLower(k) = lower(k+offset+1);
									laMain(k) = main(k+offset);
									laUpper(k) = upper(k+offset);
									laRight(k) = right(k+offset);
								}
								laMain(K-1) = main(K-1+offset);
								laRight(K-1) = right(K-1+offset);
								LaTridiagMatDouble laA(laMain,laLower,laUpper);
								LaTridiagFactDouble laAFactor;
								LaTridiagMatFactorize(laA,laAFactor);
								LaLinearSolve(laAFactor,laU,laRight);
								for(int k = 0; k < K; k++){
									u(k+offset,j) = laU(k);
								}
							}
						}
						for(int j = u.lbound(1); j <= u.ubound(1); j++){
							if(abs(j) % 2 == 1){
								Array<double,1> lower, main, upper, right;
								lower.resize(I);
								main.resize(I);
								upper.resize(I);
								right.resize(I);
								//cout << "Right size: " << right.size() << endl;
								lower = 0;
								main = 1;
								upper = 0;
								right = 0;
								for(int i = u.lbound(0); i <= u.ubound(0); i++){
									if(circ.getType(i,j) == REGULAR){
										lower(i) = 1.0/pow2(h);
										main(i) = -4.0/pow2(h);
										upper(i) = 1.0/pow2(h);
										right(i) = f(i,j) - (u(i,j-1) + u(i,j+1))/pow2(h);
									}
									else if(circ.getType(i,j) == IRREGULAR){
										Spacing sp = circ.getSpacing(i,j);
										double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
										NeighborScalar nbr = circ.getNeighbors(u,i,j);
										double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);
										lower(i) = (2.0/hW)/(hE + hW);
										main(i) = -(2.0/(hN*hS) + 2.0/(hE*hW));
										upper(i) = (2.0/hE)/(hE + hW);
										right(i) = f(i,j) - 2.0*(uN/hN + uS/hS)/(hN + hS);
										if(circ.getType(i-1,j) == EXTERIOR){
											/*cout << "***" << endl;
											cout << "At " << i << ", want to set " << i-1 << " to " << uW << endl;
											cout << "Current right value: " << right(i-1) << endl;
											cout << "In range: " << right.isInRange(i-1) << endl;
											cout << "***" << endl;*/
											right(i-1) = uW;
										}
										if(circ.getType(i+1,j) == EXTERIOR){
											/*cout << "***" << endl;
											cout << "At " << i << ", want to set " << i+1 << " to " << uE << endl;
											cout << "Current right value: " << right(i+1) << endl;
											cout << "In range: " << right.isInRange(i+1) << endl;
											cout << "***" << endl;*/
											right(i+1) = uE;
										}
									}
									//cout << right << endl;
								}


								//cout << main.size() << endl;
								//cout << "---" << endl;
								int K = main.size();
								int offset = main.lbound(0);

								LaVectorDouble laLower(K-1), laMain(K), laUpper(K-1), laRight(K), laU(K);
								laU = 0;
								for(int k = 0; k < K-1; k++){
									laLower(k) = lower(k+offset+1);
									laMain(k) = main(k+offset);
									laUpper(k) = upper(k+offset);
									laRight(k) = right(k+offset);
								}
								laMain(K-1) = main(K-1+offset);
								laRight(K-1) = right(K-1+offset);
								LaTridiagMatDouble laA(laMain,laLower,laUpper);
								LaTridiagFactDouble laAFactor;
								LaTridiagMatFactorize(laA,laAFactor);
								LaLinearSolve(laAFactor,laU,laRight);
								for(int k = 0; k < K; k++){
									u(k+offset,j) = laU(k);
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
#endif /* LINESMOOTHERS_H_ */
