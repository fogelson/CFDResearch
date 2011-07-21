/*
 * StencilTester.cpp
 *
 *  Created on: Jun 13, 2011
 *      Author: fogelson
 */

#include "Solvers/StenciledAdvectionDiffusionSolver.h"
#include "Multigrid/Smoothers.h"

using namespace CFD;
using namespace Solvers;
using namespace Geometry;
using namespace Multigrid;

int mainNope(){
	double h = .5;
	double offset = h/2;
	double r = 1;
	double deltaT = 0.00001;

	double D = 1;
	double alpha = 0;
	double u11 = 0;
	double u21 = 0;
	double u12 = 0;
	double u22 = 0;
	double H = 0;
	double Q0 = r;

	Grid * g = new Circle(h,r,offset);

	StenciledGSLex smoother;

	BilinearInterpolator bi;
	VolumeWeightedRestrictor vw;

	StenciledMultigridSolver multigrid(&smoother,&bi,&vw);

	AdvectionDiffusionStencil S(g,deltaT,D,alpha,H,Q0,u11,u12,u21,u22);
	DiffusionStencil DS(1,deltaT,g);

	CellDoubleArray u = g->makeCellDoubleArray();
	CellDoubleArray f = g->makeCellDoubleArray();

	u = where(g->cellTypes == COVERED, 0, 1);
	u(3,3) = 10;

	cout << u << endl;

	CellDoubleArray uNew = g->makeCellDoubleArray();

	for(int t = 0; t <= 100; t++){
		cout << t << endl;
		uNew = multigrid.solve(u,u,&DS,2,2,1);
		//uNew = smoother.smooth(u,u,&DS,100);
		u = uNew;
	}

	cout << u << endl;



	delete g;

	return 0;
}
