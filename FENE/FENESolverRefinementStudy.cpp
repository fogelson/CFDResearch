/*
 * FENESolverRefinementStudy.cpp
 *
 *  Created on: Aug 15, 2011
 *      Author: fogelson
 */

#include "FENEStress.h"

using namespace CFD;
using namespace Multigrid;

int main(){
	int levels = 8;

	Array<double,1> norms(levels-1);
	Array<double,1> ratios(levels-2);

	norms = 0;
	ratios = 0;

	double hFinest = 0.005;

	double H, Q0, D, deltaT, offsetFrac, offset, h;
	offsetFrac = 0.5;
	H = 0.0;
	Q0 = 1.0;
	D = 1.0;
	deltaT = 0.02;
	int timeSteps = 1;

	Array<double,2> gradU(shape(2,2));
	TinyVector<int,2> gradUIndex;
	gradUIndex(0) = 1;
	gradUIndex(1) = 1;
	gradU.reindexSelf(gradUIndex);

	gradU(1,1) = 0;
	gradU(1,2) = 0;
	gradU(2,1) = 0;
	gradU(2,2) = 0;

	Array<double,2> f, fNew;

	Grid * g, * gCurrent;
	Stencil * stencil, * stencilCurrent;
	Interpolator * interpolator;
	Restrictor * restrictor;
	StenciledSmoother * smoother;
	StenciledMultigridSolver * solver;

	g = new Circle(hFinest,Q0,hFinest/2);
	stencil = new FENEStencil(g,deltaT,D,H,Q0);

	gCurrent = g;
	stencilCurrent = stencil;

	interpolator = new BilinearInterpolator();
	restrictor = new VolumeWeightedRestrictor();
	smoother = new StenciledFourPointGS();
	solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

	fNew.resize(gCurrent->xRange,gCurrent->yRange);
	fNew = 1;

	CellDoubleArray rhs = gCurrent->makeCellDoubleArray();
	rhs = fNew;
	for(int i = 0; i < timeSteps; i++){
		solver->solve(fNew,rhs,stencilCurrent,2,2,1);
	}

	for(int l = 1; l < levels; l++){
		f.resize(gCurrent->xRange,gCurrent->yRange);
		f = fNew.copy();
		CellDoubleArray fRestricted = gCurrent->coarsen()->makeCellDoubleArray();
		fRestricted = restrictor->doRestrict(f,gCurrent,gCurrent->coarsen());

		gCurrent = gCurrent->coarsen();
		stencilCurrent = stencilCurrent->coarsen();

		fNew.resize(gCurrent->xRange,gCurrent->yRange);
		fNew = 1;

		CellDoubleArray rhs2 = gCurrent->makeCellDoubleArray();
		rhs2 = fNew;
		for(int i = 0; i < timeSteps; i++){
			solver->solve(fNew,rhs2,stencilCurrent,2,2,1);
		}

		CellDoubleArray diff = gCurrent->makeCellDoubleArray();
		diff = abs(fNew - fRestricted);

		double maxNorm = max(where(gCurrent->cellTypes == REGULAR, diff, -1));

		for(int q1 = gCurrent->iMin; q1 <= gCurrent->iMax; q1++){
			for(int q2 = gCurrent->jMin; q2 <= gCurrent->jMax; q2++){
				if(diff(q1,q2) == maxNorm)
					cout << q1 << ", " << q2 << endl;
			}
		}

		cout << "---" << endl;
		norms(l-1) = maxNorm;
	}
	for(int l = 0; l < levels-2; l++){
		ratios(l) = norms(l)/norms(l+1);
	}

	cout << norms << endl;
	cout << ratios << endl;
}
