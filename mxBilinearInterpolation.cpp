#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "Geometry/Geometry.h"
#include "Multigrid/IntergridOperators.h"
#include "Multigrid/GridOperators.h"
#include "Multigrid/Smoothers.h"

using namespace CFD;
using namespace Geometry;
using namespace Multigrid;

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [uF,uC,LuF] = mxFileTemplate(hC,offset,uC)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 3;
	int outputs = 3;

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
	double hC = getMxDouble(prhs[0]);
	double offset = getMxDouble(prhs[1]);
	CellDoubleArray uC = getMxArray(prhs[2]);
	TinyVector<int,2> ones;
	ones(0) = 1;
	ones(1) = 1;
	uC.reindexSelf(ones);

	BilinearInterpolator *bi = new BilinearInterpolator();
	Grid *coarse = new UnitSquare(hC,offset);
	double hF = hC/2;
	Grid *fine = new UnitSquare(hF,offset);
	CellDoubleArray uF = bi->doInterpolate(uC,fine,coarse);
	LaplaceOperator *L = new LaplaceOperator();

	CellDoubleArray LuF = L->apply(uF,fine);

	Restrictor *r = new VolumeWeightedRestrictor();
	uC = r->doRestrict(uF,fine,coarse);

	delete bi;
	delete coarse;
	delete fine;
	delete r;
	delete L;
	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(uF);
	plhs[1] = setMxArray(uC);
	plhs[2] = setMxArray(LuF);
}
