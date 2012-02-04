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

/* The following command should be invoked from MATLAB:
 *
 * mxMicroRefinementStudy(deltaT,h,r,t0,tF,H,D,v1,v2,vCycles)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs

	int outputs = 0;
	int inputs = 10;

	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	MexPlotTool plotter;

	double deltaT, h, r, t0, tF, H, D;
	int v1, v2, vCycles;

	deltaT = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);
	r = getMxDouble(prhs[2]);
	t0 = getMxDouble(prhs[3]);
	tF = getMxDouble(prhs[4]);
	H = getMxDouble(prhs[5]);
	D = getMxDouble(prhs[6]);
	v1 = getMxInt(prhs[7]);
	v2 = getMxInt(prhs[8]);
	vCycles = getMxInt(prhs[9]);

	double pi = acos(-1);

	double offset = h/2;
	Grid * g = new Circle(h,r,offset);

	CellDoubleArray xCentroid = g->getCellCentroidX();
	CellDoubleArray yCentroid = g->getCellCentroidY();
	CellDoubleArray xCell = g->getCellX();
	CellDoubleArray yCell = g->getCellY();

	CellDoubleArray u = g->makeCellDoubleArray();
	//u = pow2(xCentroid) + pow2(yCentroid);


	DiffusionBackwardEulerFactory * factory = new DiffusionBackwardEulerFactory(D,deltaT);

	//FENEBackwardEulerFactory * factory = new FENEBackwardEulerFactory(deltaT);

	//factory->setD(D);
	//factory->setH(H);
	//factory->setQmax(r);

	Interpolator * interpolator = new PiecewiseConstantInterpolator();
	Restrictor * restrictor = new VolumeWeightedRestrictor();
	StenciledSmoother * smoother = new GSLex();//GSFourPoint();
	MultigridSolver * solver = new MultigridSolver(smoother, interpolator, restrictor);

	CellDoubleArray zed = g->makeCellDoubleArray();
	CellDoubleArray u0 = g->makeCellDoubleArray();

	double t = t0;

	u0 = (pow2(xCell) + pow2(yCell))*exp(t);
	u = u0;

	CellDoubleArray rhs = g->makeCellDoubleArray();
	CellDoubleArray f = g->makeCellDoubleArray();

	CellDoubleArray uExact = g->makeCellDoubleArray();


	for(int n = 0; t < tF; n++){
		t = n*deltaT + t0;

		uExact = (pow2(xCentroid) + pow2(yCentroid))*exp(t);

		/*plotter.newFigure();
		plotter.graphCellCentroidData(uExact,g);
		stringstream titleExact;
		titleExact << "uExact, t = ";
		titleExact << t;
		plotter.title(titleExact.str());
		plotter.colorbar();

		plotter.newFigure();
		plotter.graphCellCentroidData(u,g);
		stringstream titleU;
		titleU << "u, t = ";
		titleU << t;
		plotter.title(titleU.str());
		plotter.colorbar();

		plotter.newFigure();
		plotter.graphCellCentroidData(factory->getLHS(g)->constantTerm,g);
		plotter.colorbar();*/

		factory->setTime(t);
		f = exp(t+deltaT)*(pow2(xCentroid)+pow2(yCentroid)-4);
		rhs = deltaT*f + factory->getRHS(g)->apply(u);
		for(int c = 0; c < vCycles; c++){
			solver->vCycle(u,u0,rhs,v1,v2,g,factory);
		}
		u0 = u;

		/*t = n*deltaT + t0;
		cout << t << endl;
		factory->setTime(t);
		f = exp(t+deltaT)*(pow2(xCentroid)+pow2(yCentroid)-4);
		rhs = deltaT*f + factory->getRHS(g)->apply(u0);
		for(int c = 0; c < vCycles; c++){
			solver->vCycle(u,u0,rhs,v1,v2,g,factory);
		}
		u0 = u;*/
	}

	uExact = (pow2(xCentroid) + pow2(yCentroid))*exp(t+deltaT);


	CellDoubleArray err = g->makeCellDoubleArray();
	err = u - uExact;

	double errl2 = sum(where(g->getCellTypes()!=COVERED, pow2(err)*(g->getVolumes()), 0));
	errl2 = errl2/sum(where(g->getCellTypes()!=COVERED, g->getVolumes(), 0));
	errl2 = sqrt(errl2);

	double maxErr = max(where(g->getCellTypes() != COVERED, abs(err), 0));
	cout << errl2 << endl;

	plotter.newFigure();
	plotter.graphCellCentroidData(err,g);
	plotter.colorbar();

	/*plotter.newFigure();
	plotter.graphCellCentroidData(u,g);
	plotter.colorbar();
	plotter.newFigure();
	plotter.graphCellCentroidData(uExact,g);
	plotter.colorbar();*/


	delete interpolator;
	delete restrictor;
	delete smoother;
	delete solver;
	delete factory;
	delete g;

//	CellDoubleArray xCC = g->getCellX();
//	CellDoubleArray yCC = g->getCellY();
//	CellDoubleArray rCC = g->makeCellDoubleArray();
//	rCC = sqrt(pow2(xCC) + pow2(yCC));
//
//	CellDoubleArray u = g->makeCellDoubleArray();
//	u = cos(2*pi*rCC);
//
//	CellDoubleArray xCT = g->getCellCentroidX();
//	CellDoubleArray yCT = g->getCellCentroidY();
//	CellDoubleArray rCT = g->makeCellDoubleArray();
//	rCT = sqrt(pow2(xCT) + pow2(yCT));
//
//	/*FaceDoubleArray xFace = g->getFaceX();
//	FaceDoubleArray yFace = g->getFaceY();
//	FaceDoubleArray theta = g->makeFaceDoubleArray();
//	theta = atan2(yFace,xFace);
//
//	FaceDoubleArray beta = g->makeFaceDoubleArray();
//	beta = 1;//where(abs(theta) >= pi/2, 3, 1);
//
//	CellToFaceOperator flux;
//	flux = beta*grad;
//
//	LaplacianFactory LFac;
//	L = *(LFac.get(g));
//	//CellToCellOperator L;
//	//L = LFac.get(g);*/
//
//	CellDoubleArray LuExact = g->makeCellDoubleArray();
//	LuExact = -2*pi*sin(2*pi*rCT)/rCT - 4*pow2(pi)*cos(2*pi*rCT);
//
//	double deltaT = .05;
//	double D = .1;
//
//	//OperatorFactory<CellToCellOperator> * factory = new AdvectionDiffusionBackwardEulerFactory(3,0,D,deltaT);
//	//AdvectionDiffusionBackwardEulerFactory  * factory = new AdvectionDiffusionBackwardEulerFactory(8,0,D,deltaT);
//	//FENEBackwardEulerFactory * factory = new FENEBackwardEulerFactory(deltaT);
//	FENESteadyFactory * factory = new FENESteadyFactory();
//	factory->setD(1);
//	factory->setH(1);
//	factory->setQmax(1);
//	factory->setGradU(0,10,0,0);
//	StenciledSmoother * smoother = new GSFourPoint();
//	Interpolator * interpolator = new PiecewiseConstantInterpolator();
//	Restrictor * restrictor = new VolumeWeightedRestrictor();
//
//	MultigridSolver * solver = new MultigridSolver(smoother,interpolator,restrictor);
//
//	CellToCellOperator * lhs, * rhs;
//
//	u = 1;//where(xCC <= .01, 1, 0);
//
//
//	plotter.newFigure();
//	CellDoubleArray uNew = g->makeCellDoubleArray();
//	CellDoubleArray f = g->makeCellDoubleArray();
//
//	lhs = factory->getLHS(g);
//	rhs = factory->getRHS(g);
//
//
//	for(int n = 0; n <= 0; n++){
//		cout << n << endl;
////		CellDoubleArray s = g->makeCellDoubleArray();
////		s = (g->getVolumes())*u;
////		double total = sum(where(g->getCellTypes() != COVERED, s, 0));
////		cout << total << endl;
//
//		lhs = factory->getLHS(g);
//		rhs = factory->getRHS(g);
//		f = (*rhs)(u);
//		//smoother->smooth(uNew,u,f,*lhs,400);
//		for(int k = 0; k < 30; k++){
//			solver->vCycle(uNew,u,f,2,2,g,factory);
//			//factory->clearMap(factory->lhs);
//			//factory->clearMap(factory->rhs);
//			//lhs = factory->getLHS(g);
//			//rhs = factory->getRHS(g);
//		}
//		u = uNew;
//		/*if(n % 1 == 0){
//			plotter.graphCellCentroidData(u,g);
////			plotter.graphCellDoubleArray(u,g,"mesh");
//			plotter.colorbar();
//			plotter.drawNow();
//		}*/
//	}
//	plotter.graphCellCentroidData(u,g);
//	plotter.colorbar();
//	plotter.drawNow();
//
//	/*OperatorFactory<CellToCellOperator> * upwind = new UpwindFactory(*solidBody);
//	CellToCellOperator * uw = upwind->get(g);
//
//	u = where(abs(xCC) <= .2 || abs(yCC) <= .2 , 0,1);
//
//	mxArray * nullplhs[0], * nullprhs[0];
//
//	CellDoubleArray uNew = g->makeCellDoubleArray();
//	plotter.initializeMovie();
//	plotter.newFigure();
//
//	for(int n = 0; n <= 6000; n++){
//		uNew = u + deltaT*((*uw)(u));
//		u = uNew;
//		if(n % 200 == 0){
//			plotter.graphCellCentroidData(u,g);
//			mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
//			//plotter.captureFrame();
//			plotter.drawNow();
//		}
//	}
//	//plotter.playMovie();*/
//
//
//
///*	SplitOperatorFactory<CellToCellOperator> * factory = new CrankNicholsonDiffusionFactory(deltaT);
//	CellToCellOperator * P = factory->getRHS(g);
//
//	StenciledSmoother * smoother = new GSFourPoint();
//	Interpolator * interpolator = new PiecewiseConstantInterpolator();
//	Restrictor * restrictor = new VolumeWeightedRestrictor();
//
//	MultigridSolver solver(smoother,interpolator,restrictor);
//
//
//	u = 0;
//
//	mxArray * nullplhs[0], * nullprhs[0];
//
//	Grid * fine = g;
//
//	//CellDoubleArray uNew = g->makeCellDoubleArray();
//	CellDoubleArray rhs = g->makeCellDoubleArray();
//
//	for(int n = 0; n < 20; n++){
//		rhs = (*P)(u) - deltaT*LuExact;
//		for(int k = 0; k < 3; k++){
//			//smoother->smooth(u,u,rhs,*(factory->getLHS(g)),100);
//			solver.vCycle(u,u,rhs,3,3,fine,factory);
//		}
//		//plotter.graphCellCentroidData(u,fine);
//		//mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
//		//plotter.drawNow();
//	}
//	//plotter.graphCellCentroidData(u,fine);
//	//mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");*/
//
///*	for(int k = 0; k < 1; k++){
//		CellToCellOperator * Afine = factory->getLHS(g);
//		smoother->smooth(uNew,u,rhs,*Afine,3);
//		CellDoubleArray Lu = g->makeCellDoubleArray();
//		Lu = (*Afine)(u);
//		CellDoubleArray rF = g->makeCellDoubleArray();
//		rF = rhs - Lu;
//		Grid * coarse = g->getCoarse();
//		CellDoubleArray rC = coarse->makeCellDoubleArray();
//		restrictor->doRestrict(rC,rF,coarse,fine);
//		CellDoubleArray eC = coarse->makeCellDoubleArray();
//		CellDoubleArray zC = coarse->makeCellDoubleArray();
//		CellToCellOperator * Acoarse = factory->getLHS(coarse);
//		smoother->smooth(eC,zC,rC,(*Acoarse),500);
//		CellDoubleArray eF = fine->makeCellDoubleArray();
//		interpolator->doInterpolate(eC,eF,coarse,fine);
//		CellDoubleArray uCorrected = fine->makeCellDoubleArray();
//		uCorrected = uNew + eF;
//		smoother->smooth(uNew,uCorrected,rhs,(*Afine),3);
//		u = uNew;
//	}*/
//
//
//
//
//	/*plotter.newFigure();
//	for(int k = 0; k < 50; k++){
//		plotter.graphCellCentroidData(u,g);
//		mexCallMATLAB(0,nullplhs,0,nullprhs,"colorbar");
//		plotter.drawNow();
//		CellDoubleArray rhs = g->makeCellDoubleArray();
//		rhs = (*P)(u) - deltaT*LuExact;
//		solver.vCycle(uNew,u,rhs,4,4,g,factory);
//		u = uNew;
//	}*/
//
//
//
//	delete g;
//	delete factory;
	//delete interpolator;
	//delete restrictor;
	//delete smoother;
	//delete factory;
}
