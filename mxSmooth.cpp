/*
 * mxSmooth.cpp
 *
 *  Created on: Aug 29, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "Multigrid/IntergridOperators.h"
#include "Geometry/Stencil.h"
#include "Multigrid/Smoothers.h"
#include "Multigrid/Stencils.h"
#include "Multigrid/MultigridSolvers.h"

using namespace blitzmatlab;
using namespace CFD;
using namespace Geometry;
using namespace Multigrid;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [uSmooth] = mxInterpolate(u,f,h,offset,its)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 5;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	/* Get data from the right hand side pointers. Note
	 * that we use pointers from prhs[] in the order
	 * we want them to appear in arguments to the MATLAB
	 * function call.
	 */
	double h = getMxDouble(prhs[2]);
	double offset = getMxDouble(prhs[3]);
	int its = getMxInt(prhs[4]);
	Grid * g = new Circle(h,1,offset);
	CellDoubleArray u = g->makeCellDoubleArray();
	u = getMxArray(prhs[0]);
	CellDoubleArray f = g->makeCellDoubleArray();
	f = getMxArray(prhs[1]);

	Stencil * stencil = new PoissonStencil(g);

	StenciledSmoother * smoother = new StenciledGSLex();

	CellDoubleArray uSmooth = g->makeCellDoubleArray();

	smoother->smooth(uSmooth,u,f,stencil,its);

	uSmooth = where(g->cellTypes == COVERED, mxGetNaN(), uSmooth);

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(uSmooth);

	delete g;
	delete stencil;
	delete smoother;
}
