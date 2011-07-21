/*
 * mxSteadyFlowFromFile.cpp
 *
 *  Created on: May 22, 2011
 *      Author: bfogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "Geometry/Geometry.h"
#include "Multigrid/IntergridOperators.h"
#include "Multigrid/GridOperators.h"
#include "Multigrid/Smoothers.h"
#include "Multigrid/MultigridSolvers.h"
#include "Solvers/VariableAdvectionDiffusionSolver.h"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

using namespace CFD;
using namespace Geometry;
using namespace Multigrid;
using namespace Solvers;

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [h,offset,u] = mxFileTemplate()
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

	char *input_buf;
	mwSize buflen;

	/* input must be a string */
	if ( mxIsChar(prhs[0]) != 1)
	  mexErrMsgTxt("Input must be a string.");

	/* input must be a row vector */
	if (mxGetM(prhs[0])!=1)
	  mexErrMsgTxt("Input must be a row vector.");

	/* get the length of the input string */
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

	/* copy the string data from prhs[0] into a C string input_ buf.    */
	input_buf = mxArrayToString(prhs[0]);

	if(input_buf == NULL)
	  mexErrMsgTxt("Could not convert input to string.");

	double h, offset, t;

	h = getMxDouble(prhs[1]);
	offset = h/2;


	ifstream infile;
	infile.open(input_buf, ios::in | ios::binary);

//	infile.read((char *)(&h), sizeof(h));
//	infile.read((char *)(&deltaT), sizeof(deltaT));
//	infile.read((char *)(&its), sizeof(its));

	offset = h/2;

	Grid * grid = new Circle(h,1,offset);
	CellDoubleArray u = grid->makeCellDoubleArray();

	infile.read((char *)(&t),sizeof(t));
	infile.read((char *)(u.data()), ((long)u.size())*sizeof(double));

	infile.close();

	u = where(grid->cellTypes == COVERED, mxGetNaN(), u);

	plhs[0] = setMxDouble(h);
	plhs[1] = setMxArray(u);
	plhs[2] = setMxDouble(t);

	mxFree(input_buf);
	delete grid;
}

