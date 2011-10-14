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
	u = cos(2*pi*rC);

	Gradient grad(g);

	FaceDoubleArray gradU = g->makeFaceDoubleArray();
	gradU = grad(u);

	FaceDoubleArray xF = g->getFaceX();
	FaceDoubleArray yF = g->getFaceY();
	FaceDoubleArray rF = g->makeFaceDoubleArray();
	rF = sqrt(pow2(xF) + pow2(yF));

	FaceDoubleArray gradUExact = g->makeFaceDoubleArray();
	for(int k = 0; k < g->faces.size(); k++){
		double nX, nY;
		nX = g->faces[k]->getNormal()(0);
		nY = g->faces[k]->getNormal()(1);
		gradUExact(k) = (-2*pi*xF(k)*sin(2*pi*rF(k))/rF(k))*nX + (-2*pi*yF(k)*sin(2*pi*rF(k))/rF(k))*nY;
	}

	FaceDoubleArray gradUErr = g->makeFaceDoubleArray();
	gradUErr = abs(gradU - gradUExact);
	gradUErr = where(g->getFaceTypes() != COVERED, gradUErr, 0);

	double gradUErrNorm = max(gradUErr);

	//cout << "gradUErrNorm = " << gradUErrNorm << endl;

	Divergence div(g);

	CellDoubleArray divGradU = g->makeCellDoubleArray();
	divGradU = div(gradU);

	CellDoubleArray divGradUExact = g->makeCellDoubleArray();
	divGradUExact = -4*pow2(pi)*cos(2*pi*rC) - 2*pi*sin(2*pi*rC)/rC;

	CellDoubleArray divGradUErr = g->makeCellDoubleArray();
	divGradUErr = abs(divGradU - divGradUExact);
	divGradUErr = where(g->getCellTypes() != COVERED, divGradUErr, 0);

	double divGradUErrNorm = max(divGradUErr);

	/*divGradUErrNorm = 0;
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isUncovered(i,j)){
				double v = g->cells(i,j)->getVolume();
				divGradUErrNorm += v*pow(divGradUErr(i,j),15.0);
			}
		}
	}
	divGradUErrNorm = pow(divGradUErrNorm,1.0/15.0);*/

	cout << "divGradUErrNorm = " << divGradUErrNorm << endl;

    /*plotter.newFigure();
	plotter.holdOn();
	//plotter.drawGrid(g);
	plotter.graphCellDoubleArray(u,g,"mesh");

	//plotter.newFigure();
	//plotter.holdOn();
	//plotter.drawGrid(g);
	//plotter.graphFaceDoubleArray(gradUErr,g);

	plotter.newFigure();
	plotter.holdOn();
	plotter.graphCellDoubleArray(divGradUErr,g,"mesh");
	//plotter.drawGrid(g);*/

	plotter.newFigure();
	plotter.holdOn();
	plotter.graphFaceDoubleArray(gradU,g);

	plotter.newFigure();
	//for(int i = g->iMin; i <= g->iMax; i++){
	{
		int i = 5;
		plotter.plotFaceDoubleArrayXLine(gradU,g,i,E);
		//mexEvalString("drawnow");
		//cout << i << endl;
		//mexEvalString("pause(0.2)");
	//}
	}


	/* This block of code produces a refinement study showing
	 * that the divergence operator is O(h) on the boundary
	 * and O(h^2) on the interior.
	 */
/*	CellDoubleArray xC = g->makeCellDoubleArray();
	CellDoubleArray yC = g->makeCellDoubleArray();
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			xC(i,j) = g->cells(i,j)->getCenter()(0);
			yC(i,j) = g->cells(i,j)->getCenter()(1);
		}
	}

	CellDoubleArray RC = g->makeCellDoubleArray();
	RC = sqrt(pow2(xC) + pow2(yC));

	FaceDoubleArray xF = g->makeFaceDoubleArray();
	FaceDoubleArray yF = g->makeFaceDoubleArray();
	for(int k = 0; k < g->faces.size(); k++){
		xF(k) = g->faces[k]->getCentroid()(0);
		yF(k) = g->faces[k]->getCentroid()(1);
	}

	FaceDoubleArray RF = g->makeFaceDoubleArray();
	RF = sqrt(pow2(xF) + pow2(yF));

	FaceDoubleArray u = g->makeFaceDoubleArray();

	for(int k = 0; k < g->faces.size(); k++){
		TinyVector<double,2> n = g->faces[k]->getNormal();
		double theta = atan2(g->faces[k]->getCentroid()(1),g->faces[k]->getCentroid()(0));
		u(k) = (-2*pi*sin(2*pi*RF(k))/(RF(k)))*(xF(k)*n(0) + yF(k)*n(1));
	}

	//Array<Type,1> faceTypes = g->getFaceTypes();
	//u = where(faceTypes != IRREGULAR, 0, 1);

	Divergence div(g);

	CellDoubleArray divUNum = g->makeCellDoubleArray();
	divUNum = div(u);

	CellDoubleArray divUExact = g->makeCellDoubleArray();
	divUExact = -2*pi*sin(2*pi*RC)/RC - 4*pow2(pi)*cos(2*pi*RC);

	CellDoubleArray divUErr = g->makeCellDoubleArray();
	divUErr = abs(divUNum - divUExact);

	Array<Type,2> cellTypes = g->getCellTypes();
	divUErr = where(cellTypes != COVERED, divUErr, 0);

	double divUErrNorm = max(divUErr);

	divUErrNorm = 0;

	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isUncovered(i,j)){
				double V = g->cells(i,j)->getVolume();
				divUErrNorm += V*pow(divUErr(i,j),1);
			}
		}
	}

	divUErrNorm = pow(divUErrNorm,1.0/1.0);

	cout << divUErrNorm << endl;*/
/*
	plotter.newFigure();
	plotter.holdOn();
	plotter.graphCellDoubleArray(divUErr,g,"pcolor");
	plotter.drawGrid(g);
//	plotter.graphFaceDoubleArray(u,g);

	plotter.newFigure();
	plotter.graphCellDoubleArray(divUNum,g,"pcolor");

	plotter.newFigure();
	plotter.graphCellDoubleArray(divUExact,g,"pcolor");
*/

	// The code below produces local truncation errors in the
	// gradient that converge at O(h^2) on the interior
	// and O(h) on the boundary
/*	CellDoubleArray xC = g->makeCellDoubleArray();
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
	for(int k = 0; k < g->faces.size(); k++){
		double nX, nY;
		nX = g->faces[k]->getNormal()(0);
		nY = g->faces[k]->getNormal()(1);
		gradUExact(k) = -2*pi*sin(2*pi*RF(k))*(xF(k)*nX/RF(k) + yF(k)*nY/RF(k));
	}

	FaceDoubleArray gradUNum = g->makeFaceDoubleArray();
	Gradient grad(g);
	gradUNum = grad(u);

	FaceDoubleArray gradUErr = g->makeFaceDoubleArray();
	gradUErr = abs(gradUNum - gradUExact);

	for(int k = 0; k < g->faces.size(); k++){
		if(g->faces[k]->isCovered()){
			gradUErr(k) = 0;
		}
		else if(g->faces[k]->isIrregular()){
			gradUErr(k) = 0;
		}
	}

	double gradUErrNorm;
	gradUErrNorm = max(gradUErr);

	cout << gradUErrNorm << endl;
*/


	delete g;
}
