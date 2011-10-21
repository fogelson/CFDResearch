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

	CellDoubleArray xC = g->getCellX();
	CellDoubleArray yC = g->getCellY();
	CellDoubleArray rC = g->makeCellDoubleArray();
	rC = sqrt(pow2(xC) + pow2(yC));

	CellDoubleArray u = g->makeCellDoubleArray();
	u = cos(2*pi*rC);//pow(xC,3)*pow(yC-2,2) + pow(sin(xC-pow(yC,2)),2);

	Gradient grad(g);
	Divergence div(g);

	CellDoubleArray Lu = div(grad(u));
	//Lu = where(g->getCellTypes() == REGULAR, Lu, 0);

	CellDoubleArray x = g->getCellCentroidX();
	CellDoubleArray y = g->getCellCentroidY();
	CellDoubleArray rCentroid = g->makeCellDoubleArray();
	rCentroid = sqrt(pow2(x) + pow2(y));

	CellDoubleArray LuExact = g->makeCellDoubleArray();
	LuExact = -2*pi*sin(2*pi*rCentroid)/rCentroid - 4*pow2(pi)*cos(2*pi*rCentroid);
	//LuExact = 2*(x*(pow(x,2)+3*pow(y-2,2)) - sin(2*(x-pow(y,2))) + (4*pow(y,2) + 1)*cos(2*(x-pow(y,2))));

	CellDoubleArray LuErr = g->makeCellDoubleArray();
	CellDoubleArray V = g->getVolumes();
	LuErr = abs(Lu - LuExact)*V/pow2(h);

	LuErr = where(g->getCellTypes() != COVERED, LuErr, 0);

	cout << max(LuErr) << endl;

	/*plotter.newFigure();
	plotter.graphCellDoubleArray(LuExact,g,"mesh");

	plotter.newFigure();
	plotter.graphCellDoubleArray(Lu,g,"mesh");*/

	delete g;
}
