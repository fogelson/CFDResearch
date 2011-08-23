/*
 * mxGetEBComponents.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: bfogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Geometry/Geometry.h"

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [out1,out2,out3] = mxFileTemplate(h,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 2;
	int outputs = 3;

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
	double h = getMxDouble(prhs[0]);
	double offset = getMxDouble(prhs[1]);

	// Do something with the data
	//CFD::Geometry::Grid * grid = new CFD::Geometry::Circle(h,1,offset);
	CFD::Geometry::Grid * grid = new CFD::Geometry::UnitSquare(h,offset);

	CellDoubleArray x = grid->centers.extractComponent(double(),0,2);
	CellDoubleArray y = grid->centers.extractComponent(double(),1,2);
	CellDoubleArray volumeFractions = grid->volumeFractions;

	CellDoubleArray areaFractions = grid->areaFractions.extractComponent(double(),B,5);

	Array<int,2> faceType = grid->faceTypes.extractComponent(int(),W,5);
	CellDoubleArray zero = grid->makeCellDoubleArray();

	CellDoubleArray faceTypeDouble = grid->makeCellDoubleArray();
	faceTypeDouble = zero + faceType;

	CellDoubleArray out1 = x;
	CellDoubleArray out2 = y;
	CellDoubleArray out3 = faceTypeDouble;

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(out1);
	plhs[1] = setMxArray(out2);
	plhs[2] = setMxArray(out3);

	delete grid;
}
