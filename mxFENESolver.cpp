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
#include "Multigrid/Stencils.h"

using namespace CFD;
using namespace Multigrid;
using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [f] = mxFileTemplate(f0,rhs,h,deltaT,D,H,Q0,v1,v1,its,timesteps,offsetFrac)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 12;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	double h, deltaT, D, H, Q0, v1, v2, its, timesteps, offsetFrac;

	h = getMxDouble(prhs[2]);
	deltaT = getMxDouble(prhs[3]);
	D = getMxDouble(prhs[4]);
	H = getMxDouble(prhs[5]);
	Q0 = getMxDouble(prhs[6]);
	v1 = getMxDouble(prhs[7]);
	v2 = getMxDouble(prhs[8]);
	its = getMxDouble(prhs[9]);
	timesteps = getMxDouble(prhs[10]);
	offsetFrac = getMxDouble(prhs[11]);

	//Grid * g = new UnitSquare(h,h*offsetFrac);
	Grid * g = new Circle(h,Q0,h*offsetFrac);
	//Stencil * s = new FENEStencil(g,deltaT,D,H,Q0);
	Stencil * s = new PoissonStencil(g);
	//Stencil * s = new BackwardEulerDiffusionStencil(g,deltaT,D);
	Interpolator * interpolator = new BilinearInterpolator();//QuadraticBoundaries();
	Restrictor * restrictor = new VolumeWeightedRestrictor();
	//StenciledSmoother * smoother = new StenciledGSLex();
	StenciledSmoother * smoother = new StenciledFourPointGS();
	StenciledMultigridSolver * solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

	CellDoubleArray f0 = g->makeCellDoubleArray();
	CellDoubleArray rhs = g->makeCellDoubleArray();
	f0 = getMxArray(prhs[0]);
	rhs = getMxArray(prhs[1]);

	CellDoubleArray f = g->makeCellDoubleArray();
	CellDoubleArray zero = g->makeCellDoubleArray();
	f = f0;

	//PoissonStencil L(g);

	for(int i = 0; i < timesteps; i++){
		//smoother->smooth(f,f0,rhs,s,its);
		CellDoubleArray fullRHS = g->makeCellDoubleArray();
		//CellDoubleArray Lf = g->makeCellDoubleArray();
		//Lf = L(f);
		fullRHS = rhs;// + f;// + (D*deltaT/2)*Lf;
		//smoother->smooth(f,f0,fullRHS,s,its);
		solver->solve(f,fullRHS,s,v1,v2,its);
	}

	// Set output pointers to the desired outputs

	f = where(g->cellTypes == COVERED, mxGetNaN(), f);
	plhs[0] = setMxArray(f);

	delete g;
	delete s;
	delete interpolator;
	delete restrictor;
	delete smoother;
	delete solver;
}
