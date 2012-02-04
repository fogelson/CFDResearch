/*
 * Multigrid.cpp
 *
 *  Created on: Oct 20, 2011
 *      Author: fogelson
 */


#include "Multigrid.h"
#include <iostream>

//#define MATLAB_DEBUG


#ifdef MATLAB_DEBUG

#include "../MexTools/MexTools.h"

#endif


using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
#ifdef MATLAB_DEBUG
using namespace CFD::OOMexTools;
#endif
namespace OOMultigrid{

MultigridSolver::MultigridSolver(StenciledSmoother * smoother, Interpolator * interpolator, Restrictor * restrictor){
	this->smoother = smoother;
	this->interpolator = interpolator;
	this->restrictor = restrictor;
}

void MultigridSolver::vCycle(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray rhs, int v1, int v2, Grid * fine, OperatorFactory<CellToCellOperator> * fac){
#ifdef MATLAB_DEBUG
	MexPlotTool plotter;
#endif

	//MexPlotTool plotter;

	//u = u0;

	CellToCellOperator * Afine = fac->getLHS(fine);
	//cout << "Got Afine" << endl;

	Grid * coarse = fine->getCoarse();
	//cout << "Got coarse" << endl;

	CellToCellOperator * Acoarse = fac->getLHS(coarse);
	//cout << "Got Acoarse" << endl;

	CellDoubleArray uNew = fine->makeCellDoubleArray();
	CellDoubleArray Lu = fine->makeCellDoubleArray();
	CellDoubleArray rF = fine->makeCellDoubleArray();
	CellDoubleArray eF = fine->makeCellDoubleArray();
	CellDoubleArray uCorrected = fine->makeCellDoubleArray();

	CellDoubleArray rC = coarse->makeCellDoubleArray();
	CellDoubleArray eC = coarse->makeCellDoubleArray();
	CellDoubleArray zC = coarse->makeCellDoubleArray();
	//cout << "Initialized arrays" << endl;
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(u0,fine);
	plotter.title("u0");
	plotter.drawNow();
#endif
	//plotter.graphCellCentroidData(u,fine);
	smoother->smooth(uNew,u0,rhs,*Afine,v1);
	//cout << "Presmoothed" << endl;
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(uNew,fine);
	plotter.title("presmoothed");
	plotter.drawNow();
#endif
	Lu = (*Afine)(uNew);
	//cout << "Applied operator" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(Lu,fine);
	//plotter.drawNow();
	rF = rhs - Lu;
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(rF,fine);
	plotter.title("fine residual");
	plotter.drawNow();
#endif
	//cout << "Computed residual" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(rF,fine);
	//plotter.drawNow();
	restrictor->doRestrict(rC,rF,coarse,fine);
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(rC,coarse);
	plotter.title("coarse residual");
	plotter.drawNow();
#endif
	//cout << "Restricted residual" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(rC,coarse);
	//plotter.drawNow();
	if(coarse->cells.size() <= 64){
		//cout << "On coarsest grid" << endl;
		smoother->smooth(eC,zC,rC,(*Acoarse),64);
		//cout << "Solved exactly" << endl;
		//plotter.newFigure();
		//plotter.graphCellCentroidData(eC,coarse);
		//plotter.drawNow();
	}
	else{
		//cout << "Not on coarsest grid, doing new vCycle" << endl;
		this->vCycle(eC,zC,rC,v1,v2,coarse,fac);
		//cout << "Solved via vCycle" << endl;
		//plotter.newFigure();
		//plotter.graphCellCentroidData(eC,coarse);
		//plotter.drawNow();
	}

#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(eC,coarse);
	plotter.title("coarse error");
	plotter.drawNow();
#endif
	interpolator->doInterpolate(eC,eF,coarse,fine);
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(eF,fine);
	plotter.title("fine error");
	plotter.drawNow();
#endif
	//cout << "Interpolated error" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(eF,fine);
	//plotter.drawNow();
	uCorrected = uNew + eF;
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(uCorrected,fine);
	plotter.title("corrected solution");
	plotter.drawNow();
#endif
	//cout << "Applied coarse grid correction" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(uCorrected,fine);
	//plotter.drawNow();
	smoother->smooth(uNew,uCorrected,rhs,(*Afine),v2);
#ifdef MATLAB_DEBUG
	plotter.newFigure();
	plotter.graphCellCentroidData(uNew,fine);
	plotter.title("postsmoothed");
	plotter.drawNow();
#endif
	//cout << "Postsmoothed" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(uNew,fine);
	//plotter.drawNow();
	u = uNew;
	//cout << "Set u to uNew" << endl;
}

}
}
