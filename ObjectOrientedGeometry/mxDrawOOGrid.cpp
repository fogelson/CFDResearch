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

	MexPlotTool plotter;

	double h, r, offset;

	h = getMxDouble(prhs[0]);
	r = getMxDouble(prhs[1]);
	offset = getMxDouble(prhs[2]);

	double pi = acos(-1);

	Circle * g = new Circle(h,r,offset);

	CellDoubleArray xC = g->makeCellDoubleArray();
	CellDoubleArray yC = g->makeCellDoubleArray();
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			xC(i,j) = g->cells(i,j)->getCenter()(0);
			yC(i,j) = g->cells(i,j)->getCenter()(1);
		}
	}
	CellDoubleArray RC = g->makeCellDoubleArray();
	RC = sqrt(pow2(xC) + pow2(yC));

	CellDoubleArray u = g->makeCellDoubleArray();
	u = cos(2*pi*RC);

	FaceDoubleArray xF = g->makeFaceDoubleArray();
	FaceDoubleArray yF = g->makeFaceDoubleArray();
	for(int k = 0; k < g->faces.size(); k++){
		xF(k) = g->faces[k]->getCentroid()(0);
		yF(k) = g->faces[k]->getCentroid()(1);
	}

	FaceDoubleArray RF = g->makeFaceDoubleArray();
	RF = sqrt(pow2(xF) + pow2(yF));

	FaceDoubleArray gradUExact = g->makeFaceDoubleArray();
	gradUExact = -2*pi*sin(2*pi*RF);

	FaceDoubleArray gradUNum = g->makeFaceDoubleArray();
	Gradient grad(g);
	gradUNum = grad(u);

	FaceDoubleArray gradUErr = g->makeFaceDoubleArray();
	gradUErr = abs(gradUNum - gradUExact);

	for(int k = 0; k < g->faces.size(); k++){
		if(g->faces[k]->isIrregular()){
			gradUErr(k) = 100;
		}
	}

	double gradUErrNorm;
	gradUErrNorm = max(gradUErr);

	cout << gradUErrNorm << endl;

	/*string graph = "graphCellDoubleArray";
	plotter.newFigure();
	plotter.graphCellDoubleArray(u,g,graph);
	plotter.holdOn();
	plotter.drawGrid(g);
	plotter.graphFaceDoubleArray(gradUExact,g);

	plotter.newFigure();
	plotter.graphCellDoubleArray(u,g,graph);
	plotter.holdOn();
	plotter.drawGrid(g);
	plotter.graphFaceDoubleArray(gradUNum,g);*/

	delete g;
}
