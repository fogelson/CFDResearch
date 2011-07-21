#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Geometry/Geometry.h"

using namespace blitzmatlab;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [x,y,volumeFractions] = mxFileTemplate(h,offset)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 2;
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
	double h = getMxDouble(prhs[0]);
	double offset = getMxDouble(prhs[1]);

	// Do something with the data
	CFD::Geometry::Grid * grid = new CFD::Geometry::Circle(h,1,offset);
	//CFD::Geometry::Grid * grid = new CFD::Geometry::UnitSquare(h,offset);

	CellDoubleArray x = grid->centers.extractComponent(double(),0,2);
	CellDoubleArray y = grid->centers.extractComponent(double(),1,2);
	CellDoubleArray volumeFractions = grid->volumeFractions;

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(x);
	plhs[1] = setMxArray(y);
	plhs[2] = setMxArray(volumeFractions);

	delete grid;
}
