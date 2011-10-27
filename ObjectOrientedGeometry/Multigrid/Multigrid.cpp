/*
 * Multigrid.cpp
 *
 *  Created on: Oct 20, 2011
 *      Author: fogelson
 */


#include "Multigrid.h"
#include "../MexTools/MexTools.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
using namespace OOOps;
using namespace OOMexTools;
namespace OOMultigrid{

MultigridSolver::MultigridSolver(StenciledSmoother * smoother, Interpolator * interpolator, Restrictor * restrictor){
	this->smoother = smoother;
	this->interpolator = interpolator;
	this->restrictor = restrictor;
}

void MultigridSolver::vCycle(CellDoubleArray & u, CellDoubleArray u0, CellDoubleArray & rhs, int v1, int v2, Grid * fine, SplitOperatorFactory<CellToCellOperator> * fac){
	MexPlotTool plotter;
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

	//plotter.graphCellCentroidData(u,fine);
	smoother->smooth(uNew,u,rhs,*Afine,v1);
	//cout << "Presmoothed" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(uNew,fine);
	//plotter.drawNow();
	Lu = (*Afine)(uNew);
	//cout << "Applied operator" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(Lu,fine);
	//plotter.drawNow();
	rF = rhs - Lu;
	//cout << "Computed residual" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(rF,fine);
	//plotter.drawNow();
	restrictor->doRestrict(rC,rF,coarse,fine);
	//cout << "Restricted residual" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(rC,coarse);
	//plotter.drawNow();
	if(coarse->cells.size() <= 100){
		//cout << "On coarsest grid" << endl;
		smoother->smooth(eC,zC,rC,(*Acoarse),100);
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
	interpolator->doInterpolate(eC,eF,coarse,fine);
	//cout << "Interpolated error" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(eF,fine);
	//plotter.drawNow();
	uCorrected = uNew + eF;
	//cout << "Applied coarse grid correction" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(uCorrected,fine);
	//plotter.drawNow();
	smoother->smooth(uNew,uCorrected,rhs,(*Afine),v2);
	//cout << "Postsmoothed" << endl;
	//plotter.newFigure();
	//plotter.graphCellCentroidData(uNew,fine);
	//plotter.drawNow();
	u = uNew;
	//cout << "Set u to uNew" << endl;
}

}
}
