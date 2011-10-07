/*
 * mxDrawGrid.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: fogelson
 */

#include "mex.h"
#include "matrix.h"
#include "../BlitzMatlab.h"
#include "Geo/Geometry.h"
#include <vector>

using namespace std;

using namespace blitzmatlab;
using namespace CFD::OOGeometry;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 *
 * mxSplit(xMin,xMax,yMin,yMax,h)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int outputs = 1;
	int inputs = 0;

	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

/*
	double xMin, xMax, yMin, yMax, h;
	xMin = getMxDouble(prhs[0]);
	xMax = getMxDouble(prhs[1]);
	yMin = getMxDouble(prhs[2]);
	yMax = getMxDouble(prhs[3]);
	h = getMxDouble(prhs[4]);

	double r = 1;
	Circle * g = new Circle(h,r,h/2);

	//g->setH(h);

	/*CellDoubleArray u;
	u.resize(g->xRange,g->yRange);
	u = 0;*/

	/*mexEvalString("figure; hold on;");

	vector<Face*>::iterator it;
	for(it = g->faces.begin(); it != g->faces.end(); it++){
		Face * curr = (*it);

		Coord a, b;
		a = curr->getA()->getCoord();
		b = curr->getB()->getCoord();
		Array<double,2> x(2,1);
		Array<double,2> y(2,1);


		x(0,0) = a(0);
		x(1,0) = b(0);
		y(0,0) = a(1);
		y(1,0) = b(1);

		int nlhsPlot = 0;
		int nrhsPlot = 2;
		mxArray * plhsPlot[0];
		mxArray * prhsPlot[2];
		prhsPlot[0] = setMxArray(x);
		prhsPlot[1] = setMxArray(y);

		if(curr->isUncovered()){
			//cout << "Plotting face from " << a << " to " << b << endl;
			mexCallMATLAB(nlhsPlot,plhsPlot,nrhsPlot,prhsPlot,"plot");
		}

	}*/

	plhs[0] = setMxInt(1);

	double r = 1, h = .2, offset = h/2;
	Circle * g = new Circle(h,r,offset);

	FaceDoubleArray fda;
	g->resizeFaceDoubleArray(fda);
	//cout << fda << endl;
	//fda.resize(400000000);
	//fda = 0;
//	cout << fda << endl;

	delete g;
	return;
}
