/*
 * mxForwardEulerAdvectionDiffusion.cpp
 *
 *  Created on: Feb 1, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "FiniteVolumeGeometry.h"
#include "BlitzMatlab.h"

using namespace blitzmatlab;
using namespace CFD::FiniteVolume::Geometry;

// [X,Y,uInit,uSol] = mxForwardEulerDifference(r,h,a,b,D,T,uInit)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Check for proper number of arguments. */
	if(nrhs > 7){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > 4){
		mexErrMsgTxt("Too many outputs.");
	}
	double r, h, a, b, D, T, deltaT;
	r = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);
	a = getMxDouble(prhs[2]);
	b = getMxDouble(prhs[3]);
	D = getMxDouble(prhs[4]);
	T = getMxDouble(prhs[5]);

	Circle c(r,h);
	CoordArray points = c.getVertexCenteredCoords();
	CoordTypeArray pointTypes = c.getVertexCenteredCoordTypes();
	Array<double,2> x = points.extractComponent(double(),0,2);
	Array<double,2> y = points.extractComponent(double(),1,2);

	Range iRange(points.lbound(0),points.ubound(0));
	Range jRange(points.lbound(1),points.ubound(1));
	Array<double,2> u0(iRange,jRange), u(iRange,jRange);
	if(nrhs == 8){
		u0 = getMxArray(prhs[6]);
	}
	else{
		u0 = 0;
		u0 = where(pointTypes == BOUNDARY, x*y*sin(5*x)*cos(y) + 1.0, u0);
	}
	u = u0;

	double hMin = h;
	for(int i = points.lbound(firstDim); i <= points.ubound(firstDim); i++){
		for(int j = points.lbound(secondDim); j <= points.ubound(secondDim); j++){
			if(pointTypes(i,j) == INTERIOR){
				double hN = y(i,j+1) - y(i,j);
				double hS = y(i,j) - y(i,j-1);
				double hE = x(i+1,j) - x(i,j);
				double hW = x(i,j) - x(i-1,j);
				hMin = min(hMin, hN);
				hMin = min(hMin, hS);
				hMin = min(hMin, hE);
				hMin = min(hMin, hW);
			}
		}
	}
	double dTDiffusion = (1.0/8.0)*h*h/(2.0);
	double mu = dTDiffusion*D/(h*h);
	deltaT = (1.0/8.0)*(h/a)*sqrt(2.0*mu);
	cout << "deltaT = " << deltaT << endl;



	double t = 0;
	for(int n = 0; t < T; n++){
		Array<double,2> uNew(iRange,jRange);
		uNew = where(pointTypes != INTERIOR, u, 0);
		for(int i = points.lbound(firstDim); i <= points.ubound(firstDim); i++){
			for(int j = points.lbound(secondDim); j <= points.ubound(secondDim); j++){
				if(pointTypes(i,j) == INTERIOR){
					double hN = y(i,j+1) - y(i,j);
					double hS = y(i,j) - y(i,j-1);
					double hE = x(i+1,j) - x(i,j);
					double hW = x(i,j) - x(i-1,j);
					uNew(i,j) = u(i,j);
					uNew(i,j) += 2.0*deltaT*D*(u(i-1,j)/(hW*(hE+hW)) + u(i+1,j)/(hE*(hE+hW)) + u(i,j-1)/(hS*(hN+hS)) + u(i,j+1)/(hN*(hN+hS)) - u(i,j)*(hN*hS + hE*hW)/(hN*hS*hE*hW));
					uNew(i,j) += 2.0*deltaT*a*(-u(i-1,j)*(hE/hW)/(hE+hW) + u(i+1,j)*(hW/hE)/(hE+hW) + u(i,j)*(hE-hW)/(hE*hW));
					uNew(i,j) += 2.0*deltaT*b*(-u(i,j-1)*(hN/hS)/(hN+hS) + u(i,j+1)*(hS/hN)/(hN+hS) + u(i,j)*(hN-hS)/(hN*hS));
//					u(i,j) = (u(i-1,j))/(hW*(hE+hW)) + (u(i+1,j))/(hE*(hE+hW)) + (u(i,j-1))/(hS*(hN+hS)) + (u(i,j+1))/(hN*(hN+hS));
//					u(i,j) *= (hN*hS*hE*hW)/(hN*hS + hE*hW);
				}
			}
		}
		t = n*deltaT;
		u = uNew;
	}

	Array<double,2> X = getDomain<double>(x,pointTypes);
	Array<double,2> Y = getDomain<double>(y,pointTypes);
	Array<double,2> uInit = getDomain<double>(u0,pointTypes);
	Array<double,2> uSol = getDomain<double>(u,pointTypes);

	/*cout << X << endl;
	cout << Y << endl;
	cout << uInit << endl;
	cout << uSol << endl;*/
	x = where(pointTypes==EXTERIOR,mxGetNaN(),x);
	y = where(pointTypes==EXTERIOR,mxGetNaN(),y);
	u0 = where(pointTypes==EXTERIOR,mxGetNaN(),u0);
	u = where(pointTypes==EXTERIOR,mxGetNaN(),u);

	plhs[0] = setMxArray(x);
	plhs[1] = setMxArray(y);
	plhs[2] = setMxArray(u0);
	plhs[3] = setMxArray(u);
}
