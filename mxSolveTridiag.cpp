/*
 * mxSolveTridiag.cpp
 *
 *  Created on: Feb 25, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
//#include "Multigrid/MultigridSolvers.h"
//#include "LinearAlgebra/Krylov.h"
#include <blitz/array.h>
#include <lapackpp/lapackpp.h>
#include "Multigrid/MultigridSolvers.h"

using namespace la;
using namespace blitzmatlab;
using namespace blitz;
using namespace CFD;
using namespace Multigrid;
using namespace Geometry;

// A = mxSolveTridiag(r,h)
// x = mxSolveTridiag(lower,center,upper,b)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int inputs = 2;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}
	double r, h;
	r = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);

	Circle circ(r,h);
	GridScalar f = circ.makeScalar();

	GridScalar X = circ.getXes();
	GridScalar Y = circ.getYes();

	f = pow2(X) - pow2(Y);

	PoissonDirectSolver pds;
	Array<double,2> A = pds.solve(f,circ);

	plhs[0] = setMxArray(A);

/*
	Array<double,1> lower = getMxVector(prhs[0]);
	Array<double,1> center = getMxVector(prhs[1]);
	Array<double,1> upper = getMxVector(prhs[2]);
	Array<double,1> b = getMxVector(prhs[3]);

	int N = center.size();

	LaVectorDouble laLower(N-1);
	LaVectorDouble laUpper(N-1);
	LaVectorDouble laCenter(N);
	LaVectorDouble laB(N);
	laLower(0) = lower(0);
	laCenter(0) = center(0);
	laUpper(0) = upper(1);
	laB(0) = b(0);
	for(int n = 1; n <= N-2; n++){
		laLower(n) = lower(n);
		laCenter(n) = center(n);
		laUpper(n) = upper(n+1);
		laB(n) = b(n);
	}
	laLower(N-1) = lower(N-1);
	laCenter(N-1) = center(N-1);
	laB(N-1) = b(N-1);

	LaTridiagMatDouble L(laCenter, laLower, laUpper);

	LaTridiagFactDouble fact;
	LaGenMatDouble laX(N,1);
	LaTridiagMatFactorize(L,fact);
	LaLinearSolve(fact, laX, laB);

	Array<double,2> x(N,1);

	for(int n = 0; n < N; n++){
		x(n,0) = laX(n,0);
	}
	plhs[0] = setMxArray(x);*/

}
