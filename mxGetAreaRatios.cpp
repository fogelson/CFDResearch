/*
 * mxGetAreaRatios.cpp
 *
 *  Created on: Sep 29, 2011
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
 *
 * [rat] = mxGetAreaRatios(h,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int inputs = 2;
	int outputs = 1;

	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}

	if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	/* Get data from the right hand side pointers. Note
	 * that we use pointers from prhs[] in the order
	 * we want them to appear in arguments to the MATLAB
	 * function call.
	 */

	double h = getMxDouble(prhs[0]);
	double offset = getMxDouble(prhs[1]);

	Circle circ(h,1,offset);
	Grid * g = &circ;

	CellDoubleArray u = g->makeCellDoubleArray();
	g->calculateRatios();

	u = g->lengthRatios;
	u = where(g->cellTypes != COVERED, u, mxGetNaN());

	plhs[0] = setMxArray(u);
}
