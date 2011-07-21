/*
 * mxConjugateGradient
 *
 *  Created on: Feb 22, 2011
 *      Author: fogelson
 */


#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Multigrid/MultigridSolvers.h"
#include "LinearAlgebra/Krylov.h"

using namespace blitzmatlab;
using namespace CFD::Geometry;
using namespace CFD::Multigrid;
using namespace CFD::LinearAlgebra;


class MultigridPreconditioner : public Solver{
	Circle *circ;
	int v1, v2, N;
	UpwindOperator UW;
	UpwindSymmetricGSLex SGSLex;
	BilinearInterpolator bi;
	HalfWeighter hw;
public:
	virtual ~MultigridPreconditioner(){}
	MultigridPreconditioner(Circle *circ, int v1, int v2, int N, double aX, double aY){
		this->circ = circ;
		this->v1 = v1;
		this->v2 = v2;
		this->N = N;
		SGSLex = UpwindSymmetricGSLex(aX,aY);
		UW = UpwindOperator(aX,aY);
	}
	virtual GridScalar solve(GridScalar in){
		MultigridSolver mg(&UW,&SGSLex,&bi,&hw);
		GridScalar u0 = circ->makeScalar();
		return mg.solve(u0,in,*circ,v1,v2,N);
	}
};

double boundary(double x, double y){
	if(x < 0 && y < 0){
		return 1;
	}
	else if(x > 0 && y > 0){
		return 1;
	}
	else{
		return 0;//sin(5.0*pow2(x))*cos(31.0*y);
	}
}

// [x,y,u] = mxPoissonMultigrid(r,h,tol,maxIts)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int inputs = 8;
	int outputs = 3;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	double r, h, tol, aX, aY;
	int maxIts, v1, v2;

	r = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);
	tol = getMxDouble(prhs[2]);
	maxIts = getMxInt(prhs[3]);
	v1 = getMxInt(prhs[4]);
	v2 = getMxInt(prhs[5]);
	aX = getMxDouble(prhs[6]);
	aY = getMxDouble(prhs[7]);

	Circle circ(r,h,(*boundary));
	Circle zCirc(r,h);

	GridScalar X = circ.getXes();
	GridScalar Y = circ.getYes();
	GridScalar f = circ.makeScalar();
	GridScalar u0 = circ.makeScalar();
	GridScalar u = circ.makeScalar();

	u = u0;

	IdentitySolver I;
	UpwindOperator UW(aX,aY);
	PCG pcg;
	MaxNorm mn;

	MultigridPreconditioner MGP(&zCirc,v1,v2,1,aX,aY);

	u = pcg.solve(&MGP,&UW,u0,f,circ,&mn,tol,maxIts);

	X = where(circ.getTypes() == EXTERIOR, mxGetNaN(),X);
	Y = where(circ.getTypes() == EXTERIOR, mxGetNaN(),Y);
	u = where(circ.getTypes() == EXTERIOR, mxGetNaN(),u);

	plhs[0] = setMxArray(X);
	plhs[1] = setMxArray(Y);
	plhs[2] = setMxArray(u);
}
