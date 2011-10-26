/*
 * MexPlotTool.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: fogelson
 */


#include "MexTools.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOMexTools{

void MexPlotTool::holdOn(){
	mexEvalString("hold on");
}

void MexPlotTool::holdOff(){
	mexEvalString("hold off");
}

void MexPlotTool::newFigure(){
	mexEvalString("figure");
}

void MexPlotTool::drawNow(){
	mexEvalString("drawnow");
}

void MexPlotTool::drawGrid(Grid * g){
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
			mexCallMATLAB(nlhsPlot,plhsPlot,nrhsPlot,prhsPlot,"plot");
		}
	}
}

void MexPlotTool::plotFaceDoubleArrayXLine(FaceDoubleArray & u, Grid * g, int i, Direction d){
	FaceDoubleArray x = g->getFaceX();
	FaceDoubleArray y = g->getFaceY();

	vector<Face*> faces;

	for(int j = g->jMin; j <= g->jMax; j++){
		if(g->isFaceUncovered(i,j,d)){
			faces.push_back(g->cells(i,j)->getFace(N));
		}
	}

	Array<double,1> yOut(faces.size());
	Array<double,1> uOut(faces.size());

	double pos = 0;

	for(int k = 0; k < faces.size(); k++){
		yOut(k) = pos;
		pos += faces[k]->getArea();
		uOut(k) = u(faces[k]->getIndex());
	}
	int nlhs = 0;
	int nrhs = 2;
	mxArray * plhs[nlhs];
	mxArray * prhs[nrhs];

	prhs[0] = setMxVector(yOut);
	prhs[1] = setMxVector(uOut);

	mexCallMATLAB(nlhs,plhs,nrhs,prhs,"plot");
}

void MexPlotTool::graphFaceDoubleArray(FaceDoubleArray & u, Grid * g){
	FaceDoubleArray x = g->getFaceX();
	FaceDoubleArray y = g->getFaceY();
	FaceDoubleArray nX = g->makeFaceDoubleArray();
	FaceDoubleArray nY = g->makeFaceDoubleArray();
	FaceDoubleArray uPlot = g->makeFaceDoubleArray();
	uPlot = mxGetNaN();
	for(int k = 0; k < g->faces.size(); k++){
		if(g->faces[k] != 0 && g->faces[k]->isUncovered()){
			//x(k) = g->faces[k]->getCentroid()(0);
			//y(k) = g->faces[k]->getCentroid()(1);
			/*Coord A, B;
			A = g->faces[k]->getA()->getCoord();
			B = g->faces[k]->getB()->getCoord();
			double x1, x2, y1, y2;
			x1 = A(0);
			x2 = B(0);
			y1 = A(1);
			y2 = B(1);
			double r = sqrt(pow2(x1 - x2) + pow2(y1 - y2));
			nX(k) = (y2 - y1)/r;
			nY(k) = (x1 - x2)/r;*/
			TinyVector<double,2> normal = g->faces[k]->getNormal();
			nX(k) = normal(0)*u(k);
			nY(k) = normal(1)*u(k);
			uPlot(k) = u(k);
		}
	}
	/*int nlhs = 0;
	int nrhs = 4;
	mxArray * plhs[nlhs];
	mxArray * prhs[nrhs];
	prhs[0] = setMxVector(x);
	prhs[1] = setMxVector(y);
	prhs[2] = setMxInt(20);
	prhs[3] = setMxVector(uPlot);
	mexCallMATLAB(nlhs,plhs,nrhs,prhs,"scatter");*/

	int nlhs = 0;
	int nrhs = 4;
	mxArray * plhs[nlhs];
	mxArray * prhs[nrhs];
	prhs[0] = setMxVector(x);
	prhs[1] = setMxVector(y);
	prhs[2] = setMxVector(nX);
	prhs[3] = setMxVector(nY);
	mexCallMATLAB(nlhs,plhs,nrhs,prhs,"quiver");
}

void MexPlotTool::graphCellCentroidData(CellDoubleArray & u, Grid * g){
	int numUncovered = sum(where(g->getCellTypes() != COVERED, 1, 0));
	cout << "There are " << numUncovered << " uncovered cells." << endl;
	int nlhs = 0;
	int nrhs = 3*numUncovered;
	mxArray * plhs[nlhs];
	mxArray * prhs[nrhs];
	int n = 0;
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			if(g->isUncovered(i,j)){
				//cout << n << ": " << 3*n << ", " << 3*n+1 << ", " << 3*n+2 << endl;

				//cout << "On uncovered cell number " << n << ", (" << i << ", " << j << ")" << endl;
				vector<Vertex*> vertices = g->cells(i,j)->getVertices();
				//cout << "Cell has " << vertices.size() << " uncovered vertices." << endl;
				Array<double,1> x(vertices.size()), y(vertices.size());
				for(int k = 0; k < vertices.size(); k++){
					x(k) = vertices[k]->getCoord()(0);
					y(k) = vertices[k]->getCoord()(1);
				}
				//cout << "\tAdded " << x << ", " << y << ", " << u(i,j) << endl;
				prhs[3*n] = setMxVector(x);
				prhs[3*n+1] = setMxVector(y);
				prhs[3*n+2] = setMxDouble(u(i,j));
				n++;
			}
		}
	}
	//cout << "Calling MATLAB" << endl;
	mexCallMATLAB(nlhs,plhs,nrhs,prhs,"fill");
}

void MexPlotTool::graphCellDoubleArray(CellDoubleArray & u, Grid * g, string graphCommand){
	CellDoubleArray x = g->makeCellDoubleArray();
	CellDoubleArray y = g->makeCellDoubleArray();
	CellDoubleArray uPlot = g->makeCellDoubleArray();
	uPlot = u;
	double h = g->getH();
	for(int i = g->iMin; i <= g->iMax; i++){
		for(int j = g->jMin; j <= g->jMax; j++){
			x(i,j) = g->cells(i,j)->getCenter()(0) - h/2;
			y(i,j) = g->cells(i,j)->getCenter()(1) - h/2;
			uPlot(i,j) = g->cells(i,j)->isUncovered() ? u(i,j) : mxGetNaN();
		}
	}
	int nlhs = 0;
	int nrhs = 3;
	mxArray * plhs[nlhs];
	mxArray * prhs[nrhs];
	prhs[0] = setMxArray(x);
	prhs[1] = setMxArray(y);
	prhs[2] = setMxArray(uPlot);
	mexCallMATLAB(nlhs,plhs,nrhs,prhs,graphCommand.c_str());
	//double plotHandle = getMxDouble(plhs[0]);
}

}
}
