/*
 * Interscale.h
 *
 *  Created on: Jan 12, 2012
 *      Author: fogelson
 */

#ifndef INTERSCALE_H_
#define INTERSCALE_H_

#include "../Geo/Geometry.h"
#include "../Multigrid/Multigrid.h"

namespace CFD{
using namespace OOGeometry;
using namespace OOMultigrid;
using namespace OOOps;
namespace Interscale{

class UniformAdvector;
class FokkerPlanckSolver;

class FokkerPlanckSolver{
	Grid * g;
	FENEBackwardEulerFactory * factory;
	Interpolator * interpolator;
	Restrictor * restrictor;
	StenciledSmoother * smoother;
	MultigridSolver * solver;

	int n;
	double deltaX, deltaY;
	double deltaT, h, offset, D, H, Q0, lambda;

	CellDoubleArray F;

public:
	Array<double,4> f;
	~FokkerPlanckSolver(){
		delete g;
		delete factory;
		delete interpolator;
		delete restrictor;
		delete solver;
	}
	FokkerPlanckSolver(int n, double deltaX, double deltaT, double h, double offset, double D, double H, double Q0, double lambda){
		this->n = n;
		this->deltaX = deltaX;
		this->deltaY = deltaX;

		this->deltaT = deltaT;
		this->h = h;
		this->offset = offset;
		this->D = D;
		this->H = H;
		this->Q0 = Q0;
		this->lambda = lambda;

		g = new Circle(h,Q0,offset);
		factory = new FENEBackwardEulerFactory(deltaT);
		factory->setD(D);
		factory->setH(H);
		factory->setQmax(Q0);
		interpolator = new OOMultigrid::PiecewiseConstantInterpolator();
		restrictor = new OOMultigrid::VolumeWeightedRestrictor();
		smoother = new OOMultigrid::GSFourPoint();
		solver = new OOMultigrid::MultigridSolver(smoother,interpolator,restrictor);

		F.resize(g->xRange,g->yRange);
		calculateSpringForce(F);

		Range nRange(0,n-1);
		f.resize(nRange,nRange,g->xRange,g->yRange);
		double A = sum(where(g->getCellTypes() != COVERED, g->getVolumes(), 0));
		f = 1.0/A;

		//cout << "f has size " << f.shape() << endl;
	}
	void solveFokkerPlanck(Array<double,3> & U){
		Array<double,2> gradU(shape(2,2));
		TinyVector<int,2> gradUIndex;
		gradUIndex(0) = 1;
		gradUIndex(1) = 1;
		gradU.reindexSelf(gradUIndex);

		int iP, iM, jP, jM;
		for(int i = 0; i < n; i++){
			iP = (i + 1 + n) % n;
			iM = (i - 1 + n) % n;
			for(int j = 0; j < n; j++){
				jP = (j + 1 + n) % n;
				jM = (j - 1 + n) % n;

				gradU(1,1) = (U(iP,j,0) - U(iM,j,0))/(2*deltaX); // Du_11
				gradU(1,2) = (U(i,jP,0) - U(i,jM,0))/(2*deltaY); // Du_12
				gradU(2,1) = (U(iP,j,1) - U(iM,j,1))/(2*deltaX); // Du_21
				gradU(2,2) = (U(i,jP,1) - U(i,jM,1))/(2*deltaY); // Du_22

				factory->setGradU(gradU(1,1),gradU(1,2),gradU(2,1),gradU(2,2));
				CellDoubleArray fij = f(i,j,g->xRange,g->yRange);

				//void solve(CellDoubleArray & u, CellDoubleArray & f, Stencil * stencil, int v1, int v2, int its){

				CellDoubleArray rhs = g->makeCellDoubleArray();
				solver->vCycle(fij,fij,rhs,2,2,g,factory);
//				solver->solve(fij,rhs,stencil,2,2,1);
				//fij = solver->solve(fij,fij,stencil,2,2,1);
				f(i,j,g->xRange,g->yRange) = fij;
			}
		}
	}

	void calculateStress(Array<double,3> & S){
		double S11, S12, S22;

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				CellDoubleArray fij = f(i,j,g->xRange,g->yRange);

				stressAtPoint(fij,S11,S12,S22);
				S(i,j,0) = S11;
				S(i,j,1) = S12;
				S(i,j,2) = S22;
			}
		}
	}

	void stressAtPoint(CellDoubleArray &f, double &S11, double &S12, double &S22){
		S11 = 0;
		S12 = 0;
		S22 = 0;
		for(int i = g->iMin; i <= g->iMax; i++){
			for(int j = g->jMin; j <= g->jMax; j++){
				if(g->isUncovered(i,j)){
					//Coord Q = g->centers(i,j);
					Coord Q = g->cells(i,j)->getCentroid();
					double dQ = g->cells(i,j)->getVolume();
					//double dQ = (g->getVolumes()(i,j))*pow2(g->h);

					S11 += f(i,j)*F(i,j)*pow2(Q(0))*dQ;
					S22 += f(i,j)*F(i,j)*pow2(Q(1))*dQ;
					S12 += f(i,j)*F(i,j)*Q(0)*Q(1)*dQ;
				}
			}
		}
		S11 *= lambda;
		S12 *= lambda;
		S22 *= lambda;
	}
	double magnitude(Coord c){
		return sqrt(pow2(c(0)) + pow2(c(1)));
	}
	void calculateSpringForce(CellDoubleArray &F){
		F = 0;
		for(int i = g->iMin; i <= g->iMax; i++){
			for(int j = g->jMin; j <= g->jMax; j++){
				if(g->isUncovered(i,j)){
					//double Q = magnitude(g->centers(i,j));//magnitude(g->cellCentroids(i,j));
					//Q = min(Q0-h/10.0,Q);
					double Q = magnitude(g->cells(i,j)->getCentroid());
					F(i,j) = H*Q;//H/(1 - pow2(Q/Q0));
					//F(i,j) = H*magnitude(g->centers(i,j));
				}
			}
		}
	}
	Grid * getGrid(){
		return g;
	}
};

class UniformAdvector{

	int n;
	double deltaX, deltaY;
	Grid * grid;
	Array<double,2> Ap;
	Array<double,2> Am;
	Array<double,2> Bp;
	Array<double,2> Bm;
	Array<double,2> F;
	Array<double,2> G;
	Array<double,2> u;
	Array<double,2> v;
	Array<double,2> uM, uP, vM, vP;
public:
	UniformAdvector(int n, double deltaX, double deltaY, Grid * grid);

	/* Figure out wave propagation for each point in physical space
	 * then apply that to all the points in configuration space.
	 *
	 * Uses the corner transport upwinding (CTU) algorithm described in
	 * section 20.5 of Leveque (FVM).
	 */
	void advectFlat(double deltaT, const Array<double,3> & U, Array<double,2> & f);
	double limiter(double theta);
	Array<double,2> limiter(Array<double,2> theta);
	Array<double,4> make4DArray();
	CellDoubleArray copyToSlice(Array<double,4> & f, int i, int j);
	void copyFromSlice(CellDoubleArray & fij, Array<double,4> & f, int i, int j);
	void advectFromFlat(double deltaT, const Array<double,3> & U, Array<double,4> & f);
	void advect(double deltaT, const Array<double,3> & U, Array<double,4> & f);
	/* U is an array of cell-centered, collocated x and y velocities.
	 * This uses averaging to compute the edge velocities at cell
	 * interfaces in the x-direction (array u) and the y direction (array v).
	 */
	void setEdgeVelocities(const Array<double,3> & U, Array<double,2> & u, Array<double,2> & v);

	/*		void riemannSolver(double & u, double & ql, double & qr, double & Am, double & Ap){
			double w = qr - ql;
			double up = max(u,0);
			double um = min(u,0);
			Ap = up*w;
			Am = um*w;
		}
		void transverseRiemannSolver(double & v, double & A, double & BmTrans, double & BpTrans){
			double vp = max(v,0);
			double vm = min(v,0);
			Bm = vm*A;
			Bp = vp*A;
		}*/
};

}
}


#endif /* INTERSCALE_H_ */
