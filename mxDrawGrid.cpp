/*
 * mxDrawGrid.cpp
 *
 *  Created on: Sep 28, 2011
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
 * mxSplit(h,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int outputs = 0;

	if(nrhs == 2){
		mexEvalString("figure; hold on");
	}
	else if(nrhs == 3){
		mxArray * plhsFig[0];
		mxArray * prhsFig[1];
		int nlhsFig = 0;
		int nrhsFig = 1;
		int figNum = getMxInt(prhs[2]);
		prhsFig[0] = setMxInt(figNum);
		mexCallMATLAB(nlhsFig,plhsFig,nrhsFig,prhsFig,"figure");
		mexEvalString("hold on");
	}
	else{
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

	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isUncovered(i,j)){
				int N = g->isRegular(i,j) ? 4 : g->numberOfVertices(i,j);
				for(int n = 0; n < N; n++){
					int nP = (n+1) % N;
					Coord c1, c2;
					c1 = g->vertices(i,j)(n);
					c2 = g->vertices(i,j)(nP);
					Array<double,2> x(2,1);
					Array<double,2> y(2,1);
					x(0,0) = c1(0);
					x(1,0) = c2(0);
					y(0,0) = c1(1);
					y(1,0) = c2(1);

					int nlhsPlot = 0;
					int nrhsPlot = 2;
					mxArray * plhsPlot[0];
					mxArray * prhsPlot[2];
					prhsPlot[0] = setMxArray(x);
					prhsPlot[1] = setMxArray(y);

					mexCallMATLAB(nlhsPlot,plhsPlot,nrhsPlot,prhsPlot,"plot");
				}
			}
		}
	}

	/*CellDoubleArray u = g->makeCellDoubleArray();
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
	plhs[4] = setMxArray(u5);*/
}
