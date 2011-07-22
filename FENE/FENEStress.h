/*
 * FENEStress.h
 *
 *  Created on: Jun 21, 2011
 *      Author: fogelson
 */

#ifndef FENESTRESS_H_
#define FENESTRESS_H_

#include "../Solvers/StenciledAdvectionDiffusionSolver.h"
#include "../Multigrid/Smoothers.h"
#include "MacroscaleObjects.h"
#include <list>
#include <map>
#include <vector>

namespace CFD{
	using namespace Solvers;
	using namespace Geometry;
	using namespace Multigrid;

	typedef TinyVector<int,2> Index;

	struct IndexCompare{
		bool operator()(const Index &I1, const Index &I2){
			if(I1(0) == I2(0))
				return (I1(1) < I2(1));
			else
				return (I1(0) < I2(0));
		}
	};

	struct MeshHierarchy{
		/*
		 * Indices of the upper and lower bounds of this mesh.
		 * These are indices with respect to the finest mesh,
		 * which is where the Stokes equation is solved.
		 */
		int iLower, iUpper, jLower, jUpper;
		int level;

		bool isCoarsest, isFinest;

		MeshHierarchy * parentMesh;

		map<Index,MeshHierarchy*,IndexCompare> childMeshes;

		~MeshHierarchy(){
			map<Index,MeshHierarchy*,IndexCompare>::iterator it = childMeshes.begin();
			while(it != childMeshes.end()){
				delete it->second;
			}
		}
		MeshHierarchy(int iLower, int iUpper, int jLower, int jUpper){
			this->iLower = iLower;
			this->iUpper = iUpper;
			this->jLower = jLower;
			this->jUpper = jUpper;

			isCoarsest = true;
			isFinest = true;

			level = 0;
		}
		void addChild(int iL, int iU, int jL, int jU){
			MeshHierarchy * child = new MeshHierarchy(iL, iU, jL, jU);
			child->parentMesh = this;
			child->isCoarsest = false;
			child->level = level + 1;
			isFinest = false;
			Index I;
			I(0) = iL;
			I(1) = jL;
			if(iLower <= iL && iUpper >= iU && jLower <= jL && jUpper >= jU){
				childMeshes[I] = child;
			}
			else{
				cout << "Tried to add child mesh whose indices exceed those of the parent mesh." << endl;
			}
		}
		void deleteChild(int iL, int iU, int jL, int jU){
			Index I;
			I(0) = iL;
			I(1) = jL;
			childMeshes.erase(I);
			if(childMeshes.size() == 0){
				isFinest = true;
			}
		}
		list<MeshHierarchy*> getFlatList(){
			list<MeshHierarchy*> l;
			map<Index,MeshHierarchy*,IndexCompare>::iterator it = childMeshes.begin();
			while(it != childMeshes.end()){
				if(it->second->isFinest){
					l.push_back(it->second);
				}
				else{
					list<MeshHierarchy*> childList = it->second->getFlatList();
					list<MeshHierarchy*>::iterator childIt = childList.begin();
					while(childIt != childList.end()){
						l.push_back(*childIt);
						childIt++;
					}
				}
				it++;
			}
			return l;
		}
	};


	template <class T>
	class ObjectHierarchy{
		int levels;
		T ** objects;
	public:
		~ObjectHierarchy(){
			for(int l = 0; l < levels; l++){
				delete objects[l];
			}
			delete objects;
		}
		ObjectHierarchy(){}
		ObjectHierarchy(int levels){
			setLevels(levels);
		}
		T* & operator [] (int l){
			return objects[l];
		}
		void setLevels(int levels){
			this->levels = levels;
			objects = new T*[levels];
		}
		int getLevels(){
			return levels;
		}
	};

	class AdaptiveStressSolver{
		int levels;

		Interpolator * interpolator;
		Restrictor * restrictor;
		StenciledSmoother * smoother;
		StenciledMultigridSolver * solver;

		int n;
		double deltaX, deltaY;
		double deltaT, hCoarsest, offset, D, H, Q0;

		vector<Grid*> grids;
		vector<Stencil*> stencils;
		vector<CellDoubleArray*> F;

		Array<CellDoubleArray*,2> f; // Spatial array of pointers to solutions of Fokker-Planck
		Array<int,2> currentLevel; // Current level of solution at f(i,j)

		SymmetricTensorArrayOrder<3> sOrder;
		VelocityArrayOrder<3> uOrder;

		Array<double,3> * S;
		Array<double,3> * U;
		Array<double,4> gradU;

	public:
		~AdaptiveStressSolver(){
			delete interpolator;
			delete restrictor;
			delete smoother;
			delete solver;

			vector<Grid*>::iterator gridsIt;
			gridsIt = grids.begin();
			while(gridsIt != grids.end()){
				delete *gridsIt;
				gridsIt++;
			}
			vector<Stencil*>::iterator stencilsIt;
			stencilsIt = stencils.begin();
			while(stencilsIt != stencils.end()){
				delete *stencilsIt;
				stencilsIt++;
			}
			vector<CellDoubleArray*>::iterator FIt;
			FIt = F.begin();
			while(FIt != F.end()){
				delete *FIt;
				FIt++;
			}

			delete S;
			delete U;

		} // End ~AdaptiveStressSolver

		void setS(double * Sdata){
			S = new Array<double,3>(Sdata,shape(n,n,3),neverDeleteData,sOrder);
		}
		void setU(double * Udata){
			U = new Array<double,3>(Udata,shape(n,n,2),neverDeleteData,uOrder);
		}

		AdaptiveStressSolver(int levels, int n, double deltaX, double deltaT, double hCoarsest, double offset, double D, double H, double Q0){
			this->levels = levels;

			this->n = n;
			this->deltaX = deltaX;
			this->deltaY = deltaX;

			this->deltaT = deltaT;
			this->hCoarsest = hCoarsest;
			this->offset = offset;
			this->D = D;
			this->H = H;
			this->Q0 = Q0;

			grids.resize(levels);
			stencils.resize(levels);
			F.resize(levels);

			if(levels == 1){
				grids[0] = new Circle(hCoarsest,Q0,offset);
				stencils[0] = new FENEStencil(grids[0],deltaT,D,H,Q0);
				F[0] = new CellDoubleArray(grids[0]->xRange,grids[0]->yRange);
			}
			else{
				double h = hCoarsest/pow2(levels-1);
				grids[levels-1] = new Circle(h,Q0,offset);
				stencils[levels-1] = new FENEStencil(grids[levels-1],deltaT,D,H,Q0);
				F[levels-1] = new CellDoubleArray(grids[levels-1]->xRange,grids[levels-1]->yRange);

				for(int l = levels - 2; l >= 0; l--){
					h = hCoarsest/pow2(l);
					grids[l] = grids[l+1]->coarsen();
					stencils[l] = new FENEStencil(grids[l],deltaT,D,H,Q0);
					F[l] = new CellDoubleArray(grids[l]->xRange,grids[l]->yRange);
				}
			}
			calculateSpringForce();

			interpolator = new BilinearInterpolator();
			restrictor = new VolumeWeightedRestrictor();
			smoother = new StenciledFourPointGS();
			solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

			Range nRange(0,n-1);
			f.resize(nRange,nRange);
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					f(i,j) = new CellDoubleArray();
					f(i,j)->resize(grids.back()->xRange,grids.back()->yRange);
					(*f(i,j)) = 1;
				}
			}
			currentLevel.resize(nRange,nRange);
			currentLevel = 0;
		} // End AdaptiveStressSolver constructor

		void refine(int i, int j){
			int c = currentLevel(i,j);
			if(c == levels - 1){
				return;
			}
			else{
				Grid * coarse = grids[c];
				Grid * fine = grids[c+1];
				CellDoubleArray * current = f(i,j);
				CellDoubleArray * refined = new CellDoubleArray(fine->xRange,fine->yRange);
				(*refined) = interpolator->doInterpolate(*current,coarse,fine);
				delete f(i,j);
				f(i,j) = refined;
				currentLevel(i,j) += 1;
			}
		} // End refine

		void coarsen(int i, int j){
			int c = currentLevel(i,j);
			if(c == 0){
				return;
			}
			else{
				Grid * fine = grids[c];
				Grid * coarse = grids[c-1];
				CellDoubleArray * current = f(i,j);
				CellDoubleArray * coarsened = new CellDoubleArray(coarse->xRange,coarse->yRange);
				(*coarsened) = restrictor->doRestrict(*current,fine,coarse);
				delete f(i,j);
				f(i,j) = coarsened;
				currentLevel(i,j) -= 1;
			}
		} // End coarsen

		void calculateSpringForce(){
			vector<Grid*>::iterator gridsIt;
			gridsIt = grids.begin();

			vector<Stencil*>::iterator stencilsIt;
			stencilsIt = stencils.begin();

			vector<CellDoubleArray*>::iterator FIt;
			FIt = F.begin();
			while(gridsIt != grids.end()){
				(*FIt) = 0;
				for(int i = (*gridsIt)->iMin; i <= (*gridsIt)->iMax; i++){
					for(int j = (*gridsIt)->jMin; i <= (*gridsIt)->jMax; j++){
						double Q = magnitude((*gridsIt)->cellCentroids(i,j));
						(*(*FIt))(i,j) = H/(1 - pow2(Q/Q0));
					}
				}

				gridsIt++;
				stencilsIt++;
				FIt++;
			}
			(*FIt) = 0;
		} // End calculateSpringForce


		double magnitude(Coord &c){
			return sqrt(pow2(c(0)) + pow2(c(1)));
		} // End magnitude

		void solvePolymerConfiguration(){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					int c = currentLevel(i,j);
					((FENEStencil*)(stencils[c]))->setGradU(gradU(i,j,Range::all(),Range::all()));
					CellDoubleArray fCurrent = *(f(i,j));
					fCurrent = solver->solve(fCurrent,fCurrent,stencils[c],2,2,1);
				}
			}
		}


/*		void updatePolymersAndCalculateStressTensor(double * Udata, double * Sdata){
#ifdef FeneTiming
			CFD::Timing::callUpdatePolymers++;
			time_t beginCallUpdatePolymers, endCallUpdatePolymers;
			time(&beginCallUpdatePolymers);
#endif
			SymmetricTensorArrayOrder<3> sOrder;
			Array<double,3> S(Sdata,shape(n,n,3),neverDeleteData,sOrder);
			VelocityArrayOrder<3> uOrder;
			Array<double,3> U(Udata,shape(n,n,2),neverDeleteData,uOrder);

			Array<double,2> gradU(shape(2,2));
			TinyVector<int,2> gradUIndex;
			gradUIndex(0) = 1;
			gradUIndex(1) = 1;
			gradU.reindexSelf(gradUIndex);

			int iP, iM, jP, jM;
			for(int i = 0; i < n; i++){
				iP = (i + 1) % n;
				iM = (i - 1) % n;
				for(int j = 0; j < n; j++){
					jP = (j + 1) % n;
					jM = (j - 1) % n;

					gradU(1,1) = (U(iP,j,1) - U(iM,j,1))/(2*deltaX); // Du_11
					gradU(1,2) = (U(i,jP,1) - U(i,jM,1))/(2*deltaY); // Du_12
					gradU(2,1) = (U(iP,j,2) - U(iM,j,2))/(2*deltaX); // Du_21
					gradU(2,2) = (U(i,jP,2) - U(i,jM,2))/(2*deltaY); // Du_22

					stencil->setGradU(gradU);
					CellDoubleArray fij = f(i,j,g->xRange,g->yRange);
#ifdef FeneTiming
					CFD::Timing::callSolver++;
					time_t begin, end;
					time(&begin);
#endif
					fij = solver->solve(fij,fij,stencil,2,2,1);
#ifdef FeneTiming
					time(&end);
					CFD::Timing::callSolverTime += difftime(end,begin);
#endif
					f(i,j,g->xRange,g->yRange) = fij;

					double S11, S12, S22;
					stressAtPoint(fij,S11,S12,S22);
					S(i,j,0) = S11;
					S(i,j,1) = S12;
					S(i,j,2) = S22;
				}
			}
#ifdef FeneTiming
			time(&endCallUpdatePolymers);
			CFD::Timing::callUpdatePolymersTime += difftime(endCallUpdatePolymers,beginCallUpdatePolymers);
#endif
		}*/


	};


	class FokkerPlanckSolver{
		Grid * g;
		FENEStencil * stencil;
		Interpolator * interpolator;
		Restrictor * restrictor;
		StenciledSmoother * smoother;
		StenciledMultigridSolver * solver;

		int n;
		double deltaX, deltaY;
		double deltaT, h, offset, D, H, Q0;

		CellDoubleArray F;

		Array<double,4> f;
	public:
		~FokkerPlanckSolver(){
			delete g;
			delete stencil;
			delete interpolator;
			delete restrictor;
			delete solver;
		}
		FokkerPlanckSolver(int n, double deltaX, double deltaT, double h, double offset, double D, double H, double Q0){
			this->n = n;
			this->deltaX = deltaX;
			this->deltaY = deltaX;

			this->deltaT = deltaT;
			this->h = h;
			this->offset = offset;
			this->D = D;
			this->H = H;
			this->Q0 = Q0;

			g = new Circle(h,Q0,offset);
			stencil = new FENEStencil(g,deltaT,D,H,Q0);
			interpolator = new BilinearInterpolator();
			restrictor = new VolumeWeightedRestrictor();
			smoother = new StenciledFourPointGS();
			solver = new StenciledMultigridSolver(smoother,interpolator,restrictor);

			F.resize(g->xRange,g->yRange);
			calculateSpringForce(F);

			Range nRange(0,n-1);
			f.resize(nRange,nRange,g->xRange,g->yRange);
			f = 1;
		}
		void updatePolymersAndCalculateStressTensor(double * Udata, double * Sdata){
			cout << "n = " << n << endl;
#ifdef FeneTiming
			CFD::Timing::callUpdatePolymers++;
			time_t beginCallUpdatePolymers, endCallUpdatePolymers;
			time(&beginCallUpdatePolymers);
#endif
			SymmetricTensorArrayOrder<3> sOrder;
			//Array<double,3> S;//(shape(n,n,3),sOrder);
			Array<double,3> S(Sdata,shape(n,n,3),neverDeleteData,sOrder);
			//S.setStorage(sOrder);
			//S.resize(n,n,3);
			VelocityArrayOrder<3> uOrder;
			Array<double,3> U(Udata,shape(n,n,2),neverDeleteData,uOrder);
			//Array<double,3> U;//(shape(n,n,2),uOrder);
			//U.setStorage(uOrder);
			//U.resize(n,n,2);

			int a = 0;

			cout << a << endl; a++;

			Array<double,2> gradU(shape(2,2));
			TinyVector<int,2> gradUIndex;
			gradUIndex(0) = 1;
			gradUIndex(1) = 1;
			gradU.reindexSelf(gradUIndex);

			cout << a << endl; a++;

			int iP, iM, jP, jM;
			for(int i = 0; i < n; i++){
				iP = (i + 1) % n;
				iM = (i - 1) % n;
				for(int j = 0; j < n; j++){
					jP = (j + 1) % n;
					jM = (j - 1) % n;

					gradU(1,1) = (U(iP,j,1) - U(iM,j,1))/(2*deltaX); // Du_11
					gradU(1,2) = (U(i,jP,1) - U(i,jM,1))/(2*deltaY); // Du_12
					gradU(2,1) = (U(iP,j,2) - U(iM,j,2))/(2*deltaX); // Du_21
					gradU(2,2) = (U(i,jP,2) - U(i,jM,2))/(2*deltaY); // Du_22

					stencil->setGradU(gradU);
					CellDoubleArray fij = f(i,j,g->xRange,g->yRange);
#ifdef FeneTiming
					CFD::Timing::callSolver++;
					time_t begin, end;
					time(&begin);
#endif
					fij = solver->solve(fij,fij,stencil,2,2,1);
#ifdef FeneTiming
					time(&end);
					CFD::Timing::callSolverTime += difftime(end,begin);
#endif
					f(i,j,g->xRange,g->yRange) = fij;

					double S11, S12, S22;
					stressAtPoint(fij,S11,S12,S22);
					S(i,j,0) = S11;
					S(i,j,1) = S12;
					S(i,j,2) = S22;
				}
			}
			cout << a << endl; a++;
#ifdef FeneTiming
			time(&endCallUpdatePolymers);
			CFD::Timing::callUpdatePolymersTime += difftime(endCallUpdatePolymers,beginCallUpdatePolymers);
#endif
		}
		void stressAtPoint(CellDoubleArray &f, double &S11, double &S12, double &S22){
#ifdef FeneTiming
			CFD::Timing::stressAtPoint++;
			time_t begin, end;
			time(&begin);
#endif
			S11 = 0;
			S12 = 0;
			S22 = 0;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					if(g->isUncovered(i,j)){
						Coord Q = g->cellCentroids(i,j);
						double dQ = (g->volumeFractions(i,j))*pow2(g->h);

						S11 += f(i,j)*F(i,j)*pow2(Q(0))*dQ;
						S22 += f(i,j)*F(i,j)*pow2(Q(1))*dQ;
						S12 += f(i,j)*F(i,j)*Q(0)*Q(1)*dQ;
					}
				}
			}
#ifdef FeneTiming
			time(&end);
			CFD::Timing::stressAtPointTime += difftime(end,begin);
#endif
		}
		double magnitude(Coord c){
			return sqrt(pow2(c(0)) + pow2(c(1)));
		}
		void calculateSpringForce(CellDoubleArray &F){
			F = 0;
			for(int i = g->iMin; i <= g->iMax; i++){
				for(int j = g->jMin; j <= g->jMax; j++){
					if(g->isUncovered(i,j)){
						double Q = magnitude(g->cellCentroids(i,j));
						F(i,j) = H/(1 - pow2(Q/Q0));
					}
				}
			}
		}
	};




/*
	void springForce(Grid * grid, double H, double Q0, CellDoubleArray &F){
		F = 0;
		for(int i = grid->iMin; i <= grid->iMax; i++){
			for(int j = grid->jMin; j <= grid->jMax; j++){
				if(grid->isUncovered(i,j)){
					double Q = magnitude(grid->cellCentroids(i,j));
					F(i,j) = H/(1 - pow2(Q/Q0));
				}
			}
		}
	}*/
}

#endif /* FENESTRESS_H_ */
