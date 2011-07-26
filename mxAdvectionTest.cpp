#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"

#include "FENE/FENEMaterialDerivative.h"
#include "FENE/MacroscaleObjects.h"

using namespace blitzmatlab;

using namespace CFD;
using namespace Geometry;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [q] = mxFileTemplate(q0,u,v,n,deltaT)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 5;
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

	Array<double,2> q = getMxArray(prhs[0]);
	Array<double,2> u = getMxArray(prhs[1]);
	Array<double,2> v = getMxArray(prhs[2]);
	int n = getMxInt(prhs[3]);
	double deltaT = getMxDouble(prhs[4]);

	TinyVector<double,2> zeroBase;
	zeroBase(0) = 0;
	zeroBase(1) = 0;
	q.reindex(zeroBase);

	double deltaX, deltaY;
	deltaX = 1.0/n;
	deltaY = deltaX;

	Grid * grid;

	UniformAdvector advector(n,deltaX,deltaY,deltaT,grid);

	Array<double,3> U;
	VelocityArrayOrder<3> uOrder;
	U.setStorage(uOrder);
	U.resize(n,n,2);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			U(i,j,0) = u(i,j);
			U(i,j,1) = v(i,j);
		}
	}

	//cout << q << endl;

	advector.advectFlat(U,q);

	//cout << q << endl;


	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(q);
}
