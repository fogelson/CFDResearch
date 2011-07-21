#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Solvers/StenciledAdvectionDiffusionSolver.h"
#include "Multigrid/Smoothers.h"

using namespace CFD;
using namespace Solvers;
using namespace Geometry;
using namespace Multigrid;

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * u = mxStenciledMultigrid(h,deltaT,T,u0,D,alpha,H,u11,u12,u21,u22)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 11;
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

	double h, deltaT, T, D, alpha, H, Q0, u11, u12, u21, u22;
	h = getMxDouble(prhs[0]);
	deltaT = getMxDouble(prhs[1]);
	T = getMxDouble(prhs[2]);

	D = getMxDouble(prhs[4]);
	alpha = getMxDouble(prhs[5]);
	H = getMxDouble(prhs[6]);

	Q0 = 1;

	u11 = getMxDouble(prhs[7]);
	u12 = getMxDouble(prhs[8]);
	u21 = getMxDouble(prhs[9]);
	u22 = getMxDouble(prhs[10]);

	Array<double,2> gradU;
	Range trivialRange(1,2);
	gradU.resize(trivialRange,trivialRange);
	gradU(1,1) = u11;
	gradU(1,2) = u12;
	gradU(2,1) = u21;
	gradU(2,2) = u22;

	gradU = alpha*gradU;

	Grid * grid = new Circle(h,Q0,h/2);
	CellDoubleArray u0 = grid->makeCellDoubleArray();
	u0 = getMxArray(prhs[3]);

	// Do something with the data

	//StenciledGSLex smoother;
	StenciledFourPointGS smoother;
	BilinearInterpolator bi;
	VolumeWeightedRestrictor vw;

	StenciledMultigridSolver multigrid(&smoother,&bi,&vw);

	//AdvectionDiffusionStencil stencil(grid,deltaT,D,alpha,H,Q0,u11,u12,u21,u22);
	//DiffusionStencil stencil(D,deltaT,grid);
	FENEStencil stencil(grid,deltaT,D,H,Q0);
	stencil.setGradU(gradU);

	/*Array<Coefficients,2> cADS, cFENE;
	grid->resizeArray<Coefficients>(cADS);
	grid->resizeArray<Coefficients>(cFENE);
	cADS = stencilADS.getCoefficients();
	cFENE = stencilFENE.getCoefficients();

	Array<double,2> cADSextracted = cADS.extractComponent(double(),C,9);
	Array<double,2> cFENEextracted = cFENE.extractComponent(double(),C,9);*/

	CellDoubleArray u = grid->makeCellDoubleArray();
	CellDoubleArray uNew = grid->makeCellDoubleArray();
	u = u0;

	//cout << u << endl;

	for(double t = 0; t <= T; t += deltaT){
//		uNew = smoother.smooth(u,u,&DS,100);
		uNew = multigrid.solve(u,u,&stencil,2,2,1);
//		stencil.setGradU(gradU);
		u = uNew;
	}

	//cout << u << endl;

	u = where(grid->cellTypes == COVERED, mxGetNaN(), u);


	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(u);

//	cADSextracted = where(grid->cellTypes == COVERED, mxGetNaN(), cADSextracted);
//	cFENEextracted = where(grid->cellTypes == COVERED, mxGetNaN(), cFENEextracted);
//
//	plhs[0] = setMxArray(cADSextracted);
//	plhs[1] = setMxArray(cFENEextracted);

	// Delete anything created with the new command
	delete grid;
}
