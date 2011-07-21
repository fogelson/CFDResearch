/*
 * mxPoissonMultigrid.cpp
 *
 *  Created on: Feb 7, 2011
 *      Author: fogelson
 */


#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Multigrid/MultigridSolvers.h"
#include "LinearAlgebra/Krylov.h"
#include <random/uniform.h>
#include <stdio.h>

using namespace blitzmatlab;
using namespace CFD::Geometry;
using namespace CFD::Multigrid;
using namespace std;
using namespace ranlib;

class MultigridPreconditioner : public Solver{
	Circle *circ;
	int v1, v2, N;
	LaplaceOperator L;
	PoissonGSRB gsRB;
	BilinearInterpolator bi;
	HalfWeighter hw;
public:
	virtual ~MultigridPreconditioner(){}
	MultigridPreconditioner(Circle *circ, int v1, int v2, int N){
		this->circ = circ;
		this->v1 = v1;
		this->v2 = v2;
		this->N = N;
	}
	virtual GridScalar solve(GridScalar in){
		MultigridSolver mg(&L,&gsRB,&bi,&hw);
		GridScalar u0 = circ->makeScalar();
		return mg.solve(u0,in,*circ,v1,v2,N);
	}
};

double boundary(double x, double y){
	return 0;
	/*if(x < 0 && y < 0){
		return 1;
	}
	else if(x > 0 && y > 0){
		return 1;
	}
	else{
		return 0;//sin(5.0*pow2(x))*cos(31.0*y);
	}*/
	//return sin(5.0*x)*cos(5.0*y);
}

// [x,y,u,n,ratios] = mxPoissonMultigrid(r,h,v1,v2,tol,N)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int inputs = 6;
	int outputs = 5;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	double r, h, tol;
	int v1, v2, N;

	r = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);
	v1 = getMxInt(prhs[2]);
	v2 = getMxInt(prhs[3]);
	tol = getMxDouble(prhs[4]);
	N = getMxInt(prhs[5]);


	LaplaceOperator L;
	FullWeighter fw;
	HalfWeighter hw;
	BilinearInterpolator bi;
	PoissonGSRB gsrb;
	PoissonGSLex gslex;
	//PoissonXLineLex pxl;
	MultigridSolver mg(&L,&gslex,&bi,&hw);
	//Coord center(0.5,0.5);
	Circle circ(r,h,(*boundary));
	GridScalar u0 = circ.makeScalar();
	GridScalar X = circ.getXes();
	GridScalar Y = circ.getYes();
	GridScalar f = circ.makeScalar();
	Uniform<double> unigen;
	u0 = 0;
	for(int i = u0.lbound(0); i <= u0.ubound(0); i++){
		for(int j = u0.lbound(1); j <= u0. ubound(1); j++){
			if(circ.getType(i,j) != EXTERIOR){
				u0(i,j) = 1;//unigen.random();
			}
		}
	}
	f = 0;
	GridScalar diff = circ.makeScalar();
	GridScalar diffNew = circ.makeScalar();
	GridScalar u = circ.makeScalar();
	GridScalar uNew = circ.makeScalar();
	u = u0;
	MaxNorm mn;
	double e = mn(u,circ);
	double eNew = e;
	int n = 0;
	Array<double,2> ratios(5,1);
	ratios = 0;
	while(e > tol && n < N){
		n++;
		//uNew = pxl.smooth(u,f,circ,1);
		uNew = mg.solve(u,f,circ,v1,v2,1);
		diffNew = uNew - u;
		eNew = mn(diffNew,circ);
		ratios(0,0) = ratios(1,0);
		ratios(1,0) = ratios(2,0);
		ratios(2,0) = ratios(3,0);
		ratios(3,0) = ratios(4,0);
		ratios(4,0) = eNew/e;
		cout << "Iteration " << n << ", error ratio is " << eNew/e << endl;
		diff = diffNew;
		u = uNew;
		e = eNew;
	}

	X = where(circ.getTypes() == EXTERIOR, mxGetNaN(),X);
	Y = where(circ.getTypes() == EXTERIOR, mxGetNaN(),Y);
	u = where(circ.getTypes() == EXTERIOR, mxGetNaN(),u);

	plhs[0] = setMxArray(X);
	plhs[1] = setMxArray(Y);
	plhs[2] = setMxArray(u);
	plhs[3] = setMxInt(n);
	plhs[4] = setMxArray(ratios);
}
