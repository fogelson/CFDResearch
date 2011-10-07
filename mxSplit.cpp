/*
 * mxSplit.cpp
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
 *
 * [uI,uP,u3,u4,u5] = mxSplit(h,offset,u)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int inputs = 3;
	int outputs = 5;

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
	u = getMxArray(prhs[2]);

	CellDoubleArray uI = g->makeCellDoubleArray(), uP = g->makeCellDoubleArray();
	CellDoubleArray u3 = g->makeCellDoubleArray();
	CellDoubleArray u4 = g->makeCellDoubleArray();
	CellDoubleArray u5 = g->makeCellDoubleArray();

	uI = where(g->cellTypes == REGULAR, u, mxGetNaN());
	uP = where(g->cellTypes == IRREGULAR, u, mxGetNaN());
	u3 = where(g->numberOfVertices == 3, uP, mxGetNaN());
	u4 = where(g->numberOfVertices == 4, uP, mxGetNaN());
	u5 = where(g->numberOfVertices == 5, uP, mxGetNaN());

	plhs[0] = setMxArray(uI);
	plhs[1] = setMxArray(uP);
	plhs[2] = setMxArray(u3);
	plhs[3] = setMxArray(u4);
	plhs[4] = setMxArray(u5);
}
