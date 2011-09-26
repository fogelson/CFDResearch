/*
 * mxRestrict.cpp
 *
 *  Created on: Aug 26, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "Multigrid/IntergridOperators.h"

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
 * [uC] = mxInterpolate(uF,hF,offset)
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
	double hF = getMxDouble(prhs[1]);
	double hC = 2*hF;
	double offset = getMxDouble(prhs[2]);
	Grid * gF = new Circle(hF,1,offset);
	CellDoubleArray uF = gF->makeCellDoubleArray();
	uF = getMxArray(prhs[0]);
	Grid * gC = new Circle(hC,1,offset);
	CellDoubleArray uC = gC->makeCellDoubleArray();

	Restrictor * restrictor = new VolumeWeightedRestrictor();

	restrictor->doRestrict(uC,uF,gF,gC);

	uC = where(gC->cellTypes == COVERED, mxGetNaN(), uC);

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(uC);

	delete gF;
	delete gC;
	delete restrictor;
}
