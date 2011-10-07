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
 * mxSplit(h,r,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int outputs = 0;
	int inputs = 3;

	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}


	double h, r, offset;

	h = getMxDouble(prhs[0]);
	r = getMxDouble(prhs[1]);
	offset = getMxDouble(prhs[2]);

	Circle * g = new Circle(h,r,offset);
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

	CellDoubleArray cda;// = Arrays::makeCellDoubleArray(g);
	cda.resize(g->xRange,g->yRange);
	cda = 0;
	FaceDoubleArray fda;
	fda.resize(g->faces.size());
	fda = 0;
	VertexDoubleArray vda;
	vda.resize(g->vertices.size());
	vda = 0;
	//FaceDoubleArray fda = Arrays::makeFaceDoubleArray(g);
	//VertexDoubleArray vda = Arrays::makeVertexDoubleArray(g);

	delete g;
	return;
}
