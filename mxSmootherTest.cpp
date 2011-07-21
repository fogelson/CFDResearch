#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "Geometry/Geometry.h"
#include "Multigrid/IntergridOperators.h"
#include "Multigrid/GridOperators.h"
#include "Multigrid/Smoothers.h"
#include "Multigrid/MultigridSolvers.h"
#include "Solvers/VariableAdvectionDiffusionSolver.h"

using namespace CFD;
using namespace Geometry;
using namespace Multigrid;
using namespace Solvers;

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [u] = mxSmootherTest(D,cX,cY,h,offset,u0,f,its,deltaT)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 9;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	/* Get data from the right hand side pointers. Note
	 * that we use pointers from prhs[] in the order
	 * we want them to appear in arguments to the MATLAB
	 * function call.
	 */

	double D = getMxDouble(prhs[0]);
	double cX = getMxDouble(prhs[1]);
	double cY = getMxDouble(prhs[2]);
	double h = getMxDouble(prhs[3]);
	double offset = getMxDouble(prhs[4]);
	CellDoubleArray u0 = getMxArray(prhs[5]);
	TinyVector<int,2> ones;
	ones(0) = 1;
	ones(1) = 1;
	u0.reindexSelf(ones);
	CellDoubleArray f = getMxArray(prhs[6]);
	f.reindexSelf(ones);
	int its = getMxInt(prhs[7]);
	double deltaT = getMxDouble(prhs[8]);

//	Smoother * smoother = new PoissonGSRB();
//	GridOperator * differenceOperator = new LaplaceOperator();
//	Interpolator * interpolator = new BilinearInterpolator();
//	Restrictor * restrictor = new VolumeWeightedRestrictor();
//
//	Grid * grid = new UnitSquare(h,offset);
//
//	MultigridSolver mg(differenceOperator,smoother,interpolator,restrictor);
//
//	CellDoubleArray u = grid->makeCellDoubleArray();
//
////	u = smoother->smooth(u0,f,grid,its);
//
//	int v1 = 2;
//	int v2 = 1;
//
//	u = mg.solve(u0, f, grid, v1, v2, its);
//
//	// Set output pointers to the desired outputs
//
//	u = where(grid->cellTypes == COVERED, mxGetNaN(), u);
//
//	delete grid;
//	delete differenceOperator;
//	delete smoother;
//	delete interpolator;
//	delete restrictor;

	//Grid * grid = new UnitSquare(h,offset);
	Grid * grid = new Circle(h,1,offset);
	CellDoubleArray u = grid->makeCellDoubleArray();


	AdvectionDiffusionSolver diffusion(deltaT);
	//AdvectionDiffusionSolver diffusion(D,deltaT,cX,cY);
	u = u0;
	for(int n = 0; n < its; n++){
		u = diffusion.step(u,f,grid);
	}

	/*for(int i = grid->iMin; i <= grid->iMax; i++){
		for(int j = grid->jMin; j <= grid->jMax; j++){
			if(grid->isUncovered(i,j)){
				TinyVector<double,9> coefficients = AdvectionStencil::getCoefficients(false,u,i,j,grid);
				u(i,j) = coefficients(N);
			}
		}
	}*/


//	BilinearInterpolator bi;
//	VolumeWeightedRestrictor vw;
//
//	//CellDoubleArray uC = vw.doRestrict(u,grid,grid->coarsen());
//	CellDoubleArray uF = grid->makeCellDoubleArray();
//	uF = bi.doInterpolate(u,grid,coarse);


	//Array<Type,2> extractedFaceTypes = grid->faceTypes.extractComponent(int(),W,5);
	//u = 1.0*extractedFaceTypes;

	//u = grid->areaFractions.extractComponent(double(),W,5);

	//u = 1.0*(grid->areaFractions.extractComponent(double(),W,5));


	u = where(grid->cellTypes == COVERED, mxGetNaN(), u);
	plhs[0] = setMxArray(u);

	delete grid;
}
