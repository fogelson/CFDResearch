#include "mex.h"
#include "matrix.h"
#include "multigrid2D.h"
#include "BlitzMatlab.h"

using namespace blitzmatlab;
using namespace multigrid2D;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Check for proper number of arguments. */
	if(nrhs != 6){
		mexErrMsgTxt("Six inputs required: initial guess for u, rhs b, v1, v2, tolerance, and max cycles.");
	}
	else if(nlhs > 2){
		mexErrMsgTxt("Number of outputs should be two.");
	}
	double tol;
	int v1, v2, maxIts;
	v1 = getMxInt(prhs[2]);
	v2 = getMxInt(prhs[3]);
	tol = getMxDouble(prhs[4]);
	maxIts = getMxInt(prhs[5]);

	Array<double,2> u(getMxArray(prhs[0]));
	Array<double,2> b(getMxArray(prhs[1]));

	int n = u.length(0)-2;
	Grid g(n);
	g.u = u;
	g.b = b;

	MultigridSolver mg(v1,v2,tol,maxIts);
	int its = mg.solve(g);

	u = g.u;

  	plhs[0] = setMxArray(u);
  	plhs[1] = setMxInt(its);
}
