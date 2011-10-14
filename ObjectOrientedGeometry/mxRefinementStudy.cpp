/*
 * mxDrawGrid.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: fogelson
 */

#include "MexTools/MexTools.h"
#include "Geo/Geometry.h"
#include "Ops/Operators.h"
#include <vector>

using namespace std;

using namespace blitzmatlab;
using namespace CFD::OOGeometry;
using namespace CFD::OOOps;
using namespace CFD::OOMexTools;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 *
 * mxRefinementStudy(hCoarsest,r,offset,levels)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int outputs = 0;
	int inputs = 4;

	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}


	double hCoarsest, r, offset;
	int refinementLevels;

	hCoarsest = getMxDouble(prhs[0]);
	r = getMxDouble(prhs[1]);
	offset = getMxDouble(prhs[2]);
	refinementLevels = getMxInt(prhs[3]);

	Array<double,1> errNorm(refinementLevels), ratios(refinementLevels-1);

	for(int l = 0; l < refinementLevels; l++){
		double h = hCoarsest*pow(2.0,((double)(l)));
		Circle * g = new Circle(h,r,offset);

		Gradient grad(g);
		Divergence div(g);

		CellDoubleArray u = g->makeCellDoubleArray();
		CellDoubleArray LuExact = g->makeCellDoubleArray();
		CellDoubleArray LuNumeric = g->makeCellDoubleArray();

		double pi = acos(-1);

		for(int i = g->iMin; i <= g->iMax; i++){
			for(int j = g->jMin; j <= g->jMax; j++){
				if(g->isUncovered(i,j)){
					double x, y;
					x = g->cells(i,j)->getCenter()(0);
					y = g->cells(i,j)->getCenter()(1);
					double r = sqrt(pow2(x) + pow2(y));
					u(i,j) = cos(2*pi*r);
					LuExact(i,j) = -2.0*pi*r*sin(2.0*pi*r)/r - 4.0*pow2(pi)*cos(2.0*pi*r);
				}
			}
		}

		FaceDoubleArray gradU = grad(u);
		LuNumeric = div(gradU);

		CellDoubleArray error = g->makeCellDoubleArray();
		error = abs(LuNumeric - LuExact);
		errNorm(l) = max(error);

		delete g;
	}
	Range first(0,refinementLevels-2);
	Range second(1,refinementLevels-1);
	ratios = errNorm(first)/errNorm(second);

	cout << errNorm << endl;
	cout << ratios << endl;
}
