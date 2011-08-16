/*
 * mxFENESolver.cpp
 *
 *  Created on: Aug 16, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "FENE/FENEStress.h"

using namespace CFD;
using namespace Multigrid;
using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [f] = mxFileTemplate(f0,h,deltaT,D,H,Q0,v1,v1,its,timesteps)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 10;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	double h, deltaT, D, H, Q0, v1, v2, its, timesteps;

	CellDoubleArray f0 = getMxArray(prhs[0]);

	h = getMxDouble(prhs[1]);
	deltaT = getMxDouble(prhs[2]);
	D = getMxDouble(prhs[3]);
	Q0 = getMxDouble(prhs[4]);
	v1 = getMxDouble(prhs[5]);
	v2 = getMxDouble(prhs[6]);
	its = getMxDouble(prhs[7]);
	timesteps = getMxDouble(prhs[8]);

	Grid * g = new Circle(h,Q0,h/2);
	Stencil * s = new FENEStencil(g,deltaT,D,H,Q0);
	Interpolator * interpolator = new BilinearInterpolator();
	Restrictor * restrictor = new VolumeWeightedRestrictor();
	StenciledSmoother * smoother = new StenciledFourPointGS();
	StenciledMultigridSolver * solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

	CellDoubleArray f = g->makeCellDoubleArray();
	CellDoubleArray zero = g->makeCellDoubleArray();
	f = f0;

	for(int i = 0; i < timesteps; i++){
		solver->solve(f,zero,s,v1,v2,its);
	}

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(f);

	delete g;
	delete s;
	delete interpolator;
	delete restrictor;
	delete smoother;
	delete solver;
}
