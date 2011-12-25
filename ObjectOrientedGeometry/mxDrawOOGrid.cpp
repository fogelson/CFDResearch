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
#include "Multigrid/Multigrid.h"
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

TinyVector<double,2> solidBody(Coord pos){
	double r = sqrt(pow2(pos(0)) + pow2(pos(1)));
	double theta = atan2(pos(1),pos(0));
	double pi = acos(-1);
	double x = pos(0), y = pos(1);
	TinyVector<double,2> out;

	out(0) = y*(1.03-y)*sin(.97+y);
	out(1) = -x*(1-.99*x)*cos(1+1.01*x);

	/*else{
		out = 0;
	}*/
	return out;
}

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

	double deltaT = .01;
	double D = .1;

	//OperatorFactory<CellToCellOperator> * factory = new AdvectionDiffusionBackwardEulerFactory(3,0,D,deltaT);
	//AdvectionDiffusionBackwardEulerFactory  * factory = new AdvectionDiffusionBackwardEulerFactory(8,0,D,deltaT);
	FENEBackwardEulerFactory * factory = new FENEBackwardEulerFactory(deltaT);
	factory->setD(D);
	factory->setH(0);
	factory->setQmax(1);
	factory->setGradU(0,0,0,0);
	StenciledSmoother * smoother = new GSFourPoint();
	Interpolator * interpolator = new PiecewiseConstantInterpolator();
	Restrictor * restrictor = new VolumeWeightedRestrictor();

	MultigridSolver * solver = new MultigridSolver(smoother,interpolator,restrictor);

	CellToCellOperator * lhs, * rhs;

	u = where(xCC <= .01, 1, 0);


	plotter.newFigure();
	CellDoubleArray uNew = g->makeCellDoubleArray();
	CellDoubleArray f = g->makeCellDoubleArray();

	for(int n = 0; n <= 50; n++){
//		CellDoubleArray s = g->makeCellDoubleArray();
//		s = (g->getVolumes())*u;
//		double total = sum(where(g->getCellTypes() != COVERED, s, 0));
//		cout << total << endl;

		if(n >= 25){
			//factory->setA(-3,0);
		}

		lhs = factory->getLHS(g);
		rhs = factory->getRHS(g);
		f = (*rhs)(u);
		//smoother->smooth(uNew,u,f,*lhs,400);
		for(int k = 0; k < 3; k++){
			solver->vCycle(uNew,u,f,2,2,g,factory);
			//factory->clearMap(factory->lhs);
			//factory->clearMap(factory->rhs);
			//lhs = factory->getLHS(g);
			//rhs = factory->getRHS(g);
		}
		u = uNew;
		if(n % 1 == 0){
			plotter.graphCellCentroidData(u,g);
//			plotter.graphCellDoubleArray(u,g,"mesh");
			plotter.colorbar();
			plotter.drawNow();
		}
	}

	/*OperatorFactory<CellToCellOperator> * upwind = new UpwindFactory(*solidBody);
	CellToCellOperator * uw = upwind->get(g);

	u = where(abs(xCC) <= .2 || abs(yCC) <= .2 , 0,1);

	mxArray * nullplhs[0], * nullprhs[0];

	CellDoubleArray uNew = g->makeCellDoubleArray();
	plotter.initializeMovie();
	plotter.newFigure();

	for(int n = 0; n <= 6000; n++){
		uNew = u + deltaT*((*uw)(u));
		u = uNew;
		if(n % 200 == 0){
			plotter.graphCellCentroidData(u,g);
			mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
			//plotter.captureFrame();
			plotter.drawNow();
		}
	}
	//plotter.playMovie();*/



/*	SplitOperatorFactory<CellToCellOperator> * factory = new CrankNicholsonDiffusionFactory(deltaT);
	CellToCellOperator * P = factory->getRHS(g);

	StenciledSmoother * smoother = new GSFourPoint();
	Interpolator * interpolator = new PiecewiseConstantInterpolator();
	Restrictor * restrictor = new VolumeWeightedRestrictor();

	MultigridSolver solver(smoother,interpolator,restrictor);


	u = 0;

	mxArray * nullplhs[0], * nullprhs[0];

	Grid * fine = g;

	//CellDoubleArray uNew = g->makeCellDoubleArray();
	CellDoubleArray rhs = g->makeCellDoubleArray();

	for(int n = 0; n < 20; n++){
		rhs = (*P)(u) - deltaT*LuExact;
		for(int k = 0; k < 3; k++){
			//smoother->smooth(u,u,rhs,*(factory->getLHS(g)),100);
			solver.vCycle(u,u,rhs,3,3,fine,factory);
		}
		//plotter.graphCellCentroidData(u,fine);
		//mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
		//plotter.drawNow();
	}
	//plotter.graphCellCentroidData(u,fine);
	//mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");*/

/*	for(int k = 0; k < 1; k++){
		CellToCellOperator * Afine = factory->getLHS(g);
		smoother->smooth(uNew,u,rhs,*Afine,3);
		CellDoubleArray Lu = g->makeCellDoubleArray();
		Lu = (*Afine)(u);
		CellDoubleArray rF = g->makeCellDoubleArray();
		rF = rhs - Lu;
		Grid * coarse = g->getCoarse();
		CellDoubleArray rC = coarse->makeCellDoubleArray();
		restrictor->doRestrict(rC,rF,coarse,fine);
		CellDoubleArray eC = coarse->makeCellDoubleArray();
		CellDoubleArray zC = coarse->makeCellDoubleArray();
		CellToCellOperator * Acoarse = factory->getLHS(coarse);
		smoother->smooth(eC,zC,rC,(*Acoarse),500);
		CellDoubleArray eF = fine->makeCellDoubleArray();
		interpolator->doInterpolate(eC,eF,coarse,fine);
		CellDoubleArray uCorrected = fine->makeCellDoubleArray();
		uCorrected = uNew + eF;
		smoother->smooth(uNew,uCorrected,rhs,(*Afine),3);
		u = uNew;
	}*/




	/*plotter.newFigure();
	for(int k = 0; k < 50; k++){
		plotter.graphCellCentroidData(u,g);
		mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
		plotter.drawNow();
		CellDoubleArray rhs = g->makeCellDoubleArray();
		rhs = (*P)(u) - deltaT*LuExact;
		solver.vCycle(uNew,u,rhs,4,4,g,factory);
		u = uNew;
	}*/



	delete g;
	delete factory;
	//delete interpolator;
	//delete restrictor;
	//delete smoother;
	//delete factory;
}
