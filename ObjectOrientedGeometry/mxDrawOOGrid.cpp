/*
 * mxDrawGrid.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: fogelson
 */

#include "MexTools/MexTools.h"
#include "Geo/Geometry.h"
#include "Ops/Operators.h"
#include "Ops/OperatorFactory.h"
#include "Smoothers/Smoothers.h"
#include "TransferOperators/TransferOperators.h"
#include <vector>

using namespace std;

using namespace blitzmatlab;
using namespace CFD::OOGeometry;
using namespace CFD::OOOps;
using namespace CFD::OOMexTools;
using namespace CFD::OOMultigrid;

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

	Grid * g = new Circle(h,r,offset);

	CellDoubleArray xCC = g->getCellX();
	CellDoubleArray yCC = g->getCellY();
	CellDoubleArray rCC = g->makeCellDoubleArray();
	rCC = sqrt(pow2(xCC) + pow2(yCC));

	CellDoubleArray u = g->makeCellDoubleArray();
	u = cos(2*pi*rCC);

	CellDoubleArray xCT = g->getCellCentroidX();
	CellDoubleArray yCT = g->getCellCentroidY();
	CellDoubleArray rCT = g->makeCellDoubleArray();
	rCT = sqrt(pow2(xCT) + pow2(yCT));

	/*FaceDoubleArray xFace = g->getFaceX();
	FaceDoubleArray yFace = g->getFaceY();
	FaceDoubleArray theta = g->makeFaceDoubleArray();
	theta = atan2(yFace,xFace);

	FaceDoubleArray beta = g->makeFaceDoubleArray();
	beta = 1;//where(abs(theta) >= pi/2, 3, 1);

	CellToFaceOperator flux;
	flux = beta*grad;

	LaplacianFactory LFac;
	L = *(LFac.get(g));
	//CellToCellOperator L;
	//L = LFac.get(g);*/

	CellDoubleArray LuExact = g->makeCellDoubleArray();
	LuExact = -2*pi*sin(2*pi*rCT)/rCT - 4*pow2(pi)*cos(2*pi*rCT);

	double deltaT = 0.01;

	CrankNicholsonDiffusionFactory CNFactory(deltaT);
	CellToCellOperator * M = CNFactory.getLHS(g), * P = CNFactory.getRHS(g);

	GSLex smoother;

	CellDoubleArray uNew = g->makeCellDoubleArray();

	u = 0;

	mxArray * nullplhs[0], * nullprhs[0];

	plotter.newFigure();
	for(int k = 0; k < 50; k++){
		plotter.graphCellCentroidData(u,g);
		mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
		plotter.drawNow();
		CellDoubleArray rhs = g->makeCellDoubleArray();
		rhs = (*P)(u) - deltaT*LuExact;
		smoother.smooth(uNew,u,rhs,*M,100);
		u = uNew;
	}
	plotter.graphCellCentroidData(u,g);
	mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");




	delete g;
}
