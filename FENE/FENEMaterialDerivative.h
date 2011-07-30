/*
 * FENEMaterialDerivative.h
 *
 *  Created on: Jul 24, 2011
 *      Author: fogelson
 */

#ifndef FENEMATERIALDERIVATIVE_H_
#define FENEMATERIALDERIVATIVE_H_

#include "MacroscaleObjects.h"
#include "../Geometry.h"

#define MaterialDerivativeDebug

namespace CFD{
	using namespace Geometry;
	class UniformAdvector{

		int n;
		double deltaX, deltaY;
		Grid * grid;
		Array<double,2> Ap;
		Array<double,2> Am;
		Array<double,2> Bp;
		Array<double,2> Bm;
		Array<double,2> F;
		Array<double,2> G;
		Array<double,2> u;
		Array<double,2> v;
		Array<double,2> uM, uP, vM, vP;
	public:
		UniformAdvector(int n, double deltaX, double deltaY, Grid * grid){
			this->n = n;
			this->deltaX = deltaX;
			this->deltaY = deltaY;
			this->grid = grid;

			Ap.resize(n,n);
			Am.resize(n,n);
			Bp.resize(n,n);
			Bm.resize(n,n);
			F.resize(n,n);
			G.resize(n,n);
			u.resize(n,n);
			v.resize(n,n);
			uM.resize(n,n);
			uP.resize(n,n);
			vM.resize(n,n);
			vP.resize(n,n);

			Ap = 0;
			Am = 0;
			Bp = 0;
			Bm = 0;
			F = 0;
			G = 0;
			u = 0;
			v = 0;

		}
		/* Figure out wave propagation for each point in physical space
		 * then apply that to all the points in configuration space.
		 *
		 * Uses the corner transport upwinding (CTU) algorithm described in
		 * section 20.5 of Leveque (FVM).
		 */
		void advectFlat(double deltaT, const Array<double,3> & U, Array<double,2> & f){

			//cout << a << endl; a++;
			// Compute velocities at cell edges
			setEdgeVelocities(U,u,v);
			//cout << a << endl; a++;
			uM = min(u,0);
			uP = max(u,0);
			vM = min(v,0);
			vP = max(v,0);
			//cout << a << endl; a++;

			Ap = 0;
			Am = 0;
			Bp = 0;
			Bm = 0;
			F = 0;
			G = 0;
			//cout << a << endl; a++;


			Array<double,2> fNew;
			Array<double,2> wX, wY, wXLimited, wYLimited;
			fNew.resize(n,n);
			wX.resize(n,n);
			wY.resize(n,n);
			wXLimited.resize(n,n);
			wYLimited.resize(n,n);

			//cout << a << endl; a++;
			fNew = 0;
			wX = 0;
			wY = 0;
			wXLimited = 0;
			wYLimited = 0;
			//cout << a << endl; a++;

			int iM, iP, jM, jP;

			// Solve each Riemann problem
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					//cout << i << ", " << j << ", " << iP << ", " << iM << ", " << jP << ", " << jM << endl;

					wX(i,j) = f(i,j) - f(iM,j);
					wY(i,j) = f(i,j) - f(i,jM);

					Ap(i,j) = uP(i,j)*wX(i,j);
					Am(i,j) = uM(i,j)*wX(i,j);
					Bp(i,j) = vP(i,j)*wY(i,j);
					Bm(i,j) = vM(i,j)*wY(i,j);

					G(iM,j) += (-0.5*deltaT/deltaX)*vM(iM,j)*uM(i,j)*wX(i,j);
					G(iM,jP) += (-0.5*deltaT/deltaX)*vP(iM,jP)*uM(i,j)*wX(i,j);
					G(i,j) += (-0.5*deltaT/deltaX)*vM(i,j)*uP(i,j)*wX(i,j);
					G(i,jP) += (-0.5*deltaT/deltaX)*vP(i,jP)*uP(i,j)*wX(i,j);

					F(i,jM) += (-0.5*deltaT/deltaY)*uM(i,jM)*vM(i,j)*wY(i,j);
					F(iP,jM) += (-0.5*deltaT/deltaY)*uP(iP,jM)*vM(i,j)*wY(i,j);
					F(i,j) += (-0.5*deltaT/deltaY)*uM(i,j)*vP(i,j)*wY(i,j);
					F(iP,j) += (-0.5*deltaT/deltaY)*uP(iP,j)*vP(i,j)*wY(i,j);
				}
			}
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					double thetaX = 0, thetaY = 0;
					if(u(i,j) >= 0 && wX(i,j) != 0){
						thetaX = wX(iM,j)/wX(i,j);
					}
					else if(wX(i,j) != 0){
						thetaX = wX(iP,j)/wX(i,j);
					}
					if(v(i,j) >= 0 && wY(i,j) != 0){
						thetaY = wY(i,jM)/wY(i,j);
					}
					else if(wY(i,j) != 0){
						thetaY = wY(i,jP)/wY(i,j);
					}
					wXLimited(i,j) = wX(i,j)*limiter(thetaX);
					wYLimited(i,j) = wY(i,j)*limiter(thetaY);

					F(i,j) += (0.5)*abs(u(i,j))*(1.0 - (deltaT/deltaX)*abs(u(i,j)))*wXLimited(i,j);
					G(i,j) += (0.5)*abs(v(i,j))*(1.0 - (deltaT/deltaY)*abs(v(i,j)))*wYLimited(i,j);
				}
			}
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					fNew(i,j) = f(i,j)
												- (deltaT/deltaX)*(Ap(i,j) + Am(iP,j))
												- (deltaT/deltaY)*(Bp(i,j) + Bm(i,jP))
												- (deltaT/deltaX)*(F(iP,j) - F(i,j))
												- (deltaT/deltaY)*(G(i,jP) - G(i,j));
				}
			}

			//cout << a << endl; a++;
			//f = fNew;
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					f(i,j) = fNew(i,j);
					//cout << "(" << i << ", " << j << ")." << endl;
				}
			}
			//cout << a << endl; a++;

		}

		double limiter(double theta){
			return 0.0; // Upwind
			//return 1.0; // Lax-Wendroff
			//return max(max(0.0,min(1.0,2*theta)),min(2.0,theta)); // Superbee
			//return max(0.0,min(1.0,theta));
			//return max(0.0,min(min((1.0 + theta)/2.0,2.0),2.0*theta)); // MC
		}
		Array<double,2> limiter(Array<double,2> theta){
			Array<double,2> out(Range(0,n-1),Range(0,n-1));
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					out(i,j) = limiter(theta(i,j));
				}
			}
			return out;
		}

		Array<double,4> make4DArray(){
			Array<double,4> out(Range(0,n-1),Range(0,n-1),grid->xRange,grid->yRange);
			out = 0;
			return out;
		}

		void advect(double deltaT, const Array<double,3> & U, Array<double,4> & f){
			int iM, iP, jM, jP;

			// Compute velocities at cell edges
			setEdgeVelocities(U,u,v);

			double uMax, vMax;
			uMax = max(u);
			vMax = max(v);
			double uCFL, vCFL, CFL;
			uCFL = abs(uMax*deltaT/deltaX);
			vCFL = abs(vMax*deltaT/deltaY);
			CFL = max(uCFL,vCFL);
			if(CFL >= 0.9){
				advect(deltaT/2,U,f);
				advect(deltaT/2,U,f);
				return;
			}
#ifdef MaterialDerivativeDebug
			cout << "uMax = " << uMax << ". vMax = " << vMax << endl;
			cout << "deltaT = " << deltaT <<". deltaX = " << deltaX << ". deltaY = " << deltaY << endl;
			cout << "Advecting f with a CFL number of " << CFL << endl;
#endif

			uM = min(u,0);
			uP = max(u,0);
			vM = min(v,0);
			vP = max(v,0);

			Range all = Range::all();
			Range nRange(0,n-1);

			Array<double,4> wX = make4DArray();
			Array<double,4> wY = make4DArray();
			Array<double,4> F = make4DArray();
			Array<double,4> G = make4DArray();
			Array<double,4> fNew = make4DArray();

			// Compute waves
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					wX(i,j,all,all) = f(i,j,all,all) - f(iM,j,all,all);
					wY(i,j,all,all) = f(i,j,all,all) - f(i,jM,all,all);
				}
			}

			// Compute corner transport and limiter flux corrections
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					G(iM,j,all,all) += (-0.5*deltaT/deltaX)*vM(iM,j)*uM(i,j)*wX(i,j,all,all);
					G(iM,jP,all,all) += (-0.5*deltaT/deltaX)*vP(iM,jP)*uM(i,j)*wX(i,j,all,all);
					G(i,j,all,all) += (-0.5*deltaT/deltaX)*vM(i,j)*uP(i,j)*wX(i,j,all,all);
					G(i,jP,all,all) += (-0.5*deltaT/deltaX)*vP(i,jP)*uP(i,j)*wX(i,j,all,all);

					F(i,jM,all,all) += (-0.5*deltaT/deltaY)*uM(i,jM)*vM(i,j)*wY(i,j,all,all);
					F(iP,jM,all,all) += (-0.5*deltaT/deltaY)*uP(iP,jM)*vM(i,j)*wY(i,j,all,all);
					F(i,j,all,all) += (-0.5*deltaT/deltaY)*uM(i,j)*vP(i,j)*wY(i,j,all,all);
					F(iP,j,all,all) += (-0.5*deltaT/deltaY)*uP(iP,j)*vP(i,j)*wY(i,j,all,all);


					// theta is ratio of upwind edge jump to current edge jump
					Array<double,2> thetaX = grid->makeCellDoubleArray();
					Array<double,2> thetaY = grid->makeCellDoubleArray();

					if(u(i,j) >= 0){
						thetaX = where(wX(i,j,all,all) != 0, wX(iM,j,all,all)/wX(i,j,all,all), 0);
					}
					else{
						thetaX = where(wX(i,j,all,all) != 0, wX(iP,j,all,all)/wX(i,j,all,all), 0);
					}
					if(v(i,j) >= 0){
						thetaY = where(wY(i,j,all,all) != 0, wY(i,jM,all,all)/wY(i,j,all,all), 0);
					}
					else{
						thetaY = where(wY(i,j,all,all) != 0, wY(i,jP,all,all)/wY(i,j,all,all), 0);
					}

					F(i,j,all,all) += (0.5)*abs(u(i,j))*(1.0 - (deltaT/deltaX)*abs(u(i,j)))*limiter(thetaX)*wX(i,j,all,all);
					G(i,j,all,all) += (0.5)*abs(v(i,j))*(1.0 - (deltaT/deltaY)*abs(v(i,j)))*limiter(thetaY)*wY(i,j,all,all);
				}
			}

			// Advance solution

			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				iP = (i + n + 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					jP = (j + n + 1) % n;

					fNew(i,j,all,all) = f(i,j,all,all)
												- (deltaT/deltaX)*(uP(i,j)*wX(i,j,all,all) + uM(iP,j)*wX(iP,j,all,all))
												- (deltaT/deltaY)*(vP(i,j)*wY(i,j,all,all) + uM(i,jP)*wY(i,jP,all,all))
												- (deltaT/deltaX)*(F(iP,j,all,all) - F(i,j,all,all))
												- (deltaT/deltaY)*(G(i,jP,all,all) - G(i,j,all,all));
				}
			}

			f = fNew;

		}
		/* U is an array of cell-centered, collocated x and y velocities.
		 * This uses averaging to compute the edge velocities at cell
		 * interfaces in the x-direction (array u) and the y direction (array v).
		 */
		void setEdgeVelocities(const Array<double,3> & U, Array<double,2> & u, Array<double,2> & v){
			int iM, jM;
			for(int i = 0; i < n; i++){
				iM = (i + n - 1) % n;
				for(int j = 0; j < n; j++){
					jM = (j + n - 1) % n;
					u(i,j) = 0.5*(U(iM,j,0) + U(i,j,0));
					v(i,j) = 0.5*(U(i,jM,1) + U(i,j,1));
				}
			}
		}
		/*		void riemannSolver(double & u, double & ql, double & qr, double & Am, double & Ap){
				double w = qr - ql;
				double up = max(u,0);
				double um = min(u,0);
				Ap = up*w;
				Am = um*w;
			}
			void transverseRiemannSolver(double & v, double & A, double & BmTrans, double & BpTrans){
				double vp = max(v,0);
				double vm = min(v,0);
				Bm = vm*A;
				Bp = vp*A;
			}*/
	};
}

#endif /* FENEMATERIALDERIVATIVE_H_ */
