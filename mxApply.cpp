/*
 * mxApply.cpp
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
 * [Lu] = mxInterpolate(u,h,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 3;
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
	double h = getMxDouble(prhs[1]);
	double offset = getMxDouble(prhs[2]);

	//Grid * g = new Circle(h,1,offset);
	Grid * g = new UnitSquare(h,offset);
	CellDoubleArray u = g->makeCellDoubleArray();
	u = getMxArray(prhs[0]);
	CellDoubleArray Lu = g->makeCellDoubleArray();

	Stencil * stencil = new PoissonStencil(g);

	Lu = stencil->apply(u);

	Lu = where(g->cellTypes == COVERED, mxGetNaN(), Lu);

	//Lu = where(g->numberOfVertices == 3, mxGetNaN(), Lu);

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(Lu);

	delete g;
	delete stencil;
}
/*
 * mxApply.cpp
 *
 *  Created on: Aug 29, 2011
 *      Author: fogelson
 */

