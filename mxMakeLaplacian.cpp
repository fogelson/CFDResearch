/*
 * mxMakeLaplacian.cpp
 *
 *  Created on: Mar 6, 2011
 *      Author: bfogelson
 */

#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Multigrid/MultigridSolvers.h"
#include "LinearAlgebra/Krylov.h"
#include <stdio.h>

using namespace blitzmatlab;
using namespace CFD::Geometry;
using namespace CFD::Multigrid;
using namespace std;

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
/*	if(x < 0 && y < 0){
		return 1;
	}
	else if(x > 0 && y > 0){
		return 1;
	}
	else{
		return 0;//sin(5.0*pow2(x))*cos(31.0*y);
	}*/
}

// L = mxPoissonMultigrid(r,h,I)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int inputs = 3;
	int outputs = 1;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	double r, h;
	int I;

	r = getMxDouble(prhs[0]);
	h = getMxDouble(prhs[1]);
	I = getMxInt(prhs[2]);

	Circle circ(r,h);

	GridScalar u = circ.makeScalar();

	int J = u.length(1);
	int M = circ.getM();
	LaVectorDouble lower(J-1), center(J), upper(J-1);
	LaVectorDouble b(J);
	b = 0;
	lower = 0;
	center = 1;
	upper = 0;
	for(int j = u.lbound(1); j <= u.ubound(1); j++){
		if(circ.getType(I,j) == REGULAR){
			lower(j+M-1) = 1.0/pow2(h);
			upper(j+M) = 1.0/pow2(h);
			center(j+M) = -4.0/pow2(h);
			//b(j+M) = -u(I+1,j)/pow2(h) - u(I-1,j)/pow2(h) + f(I,j);
		}
			else if(circ.getType(I,j) == IRREGULAR){
			Spacing sp = circ.getSpacing(I,j);
			double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST);
			NeighborScalar nbr = circ.getNeighbors(u,I,j);
			double uN = nbr(NORTH), uS = nbr(SOUTH), uE = nbr(EAST), uW = nbr(WEST);

			lower(j+M-1) = 2.0/(hN*(hN + hS));
			upper(j+M) = 2.0/(hS*(hN + hS));
			center(j+M) = -2.0/(hN*hS) - 2.0/(hE*hW);

			/*b(j+M) = -2.0*uE/(hE*(hE + hW)) - 2.0*uW/(hW*(hE + hW)) + f(i,j);
			if(circ.getType(i,j+1) == EXTERIOR){
				b(j+M+1) = uN;
			}*/
		}

	}

	LaTridiagMatDouble L(center, lower, upper);
	/*LaGenMatDouble X(J,1);
	LaTridiagFactDouble fact;
	LaTridiagMatFactorize(L,fact);
	LaLinearSolve(fact, X, b);
	for(int j = u.lbound(1); j <= u.ubound(1); j++){
		if(circ.getType(i,j) != EXTERIOR){
			u(i,j) = X(j+M,0);
		}
	}*/
	Array<double,2> LBlitz(J,J);
	for(int i = 0; i < J; i++){
		for(int j = 0; j < J; j++){
			LBlitz(i,j) = L(i,j);
		}
	}

	plhs[0] = setMxArray(LBlitz);
}
