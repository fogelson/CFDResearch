/*
 * OperatorFactory.h
 *
 *  Created on: Oct 25, 2011
 *      Author: fogelson
 */

#ifndef OPERATORFACTORY_H_
#define OPERATORFACTORY_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <map>
#include <blitz/tinyvec-et.h>
#include "Operators.h"

using namespace blitz;
using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

class FaceCoefficientContainer;

class FaceCoefficientContainer{
	map<Grid*,FaceDoubleArray*> coefficients;
public:
	~FaceCoefficientContainer();
	void clear();
	FaceDoubleArray & get(Grid * g);
	FaceDoubleArray & operator() (Grid * g);
	bool contains(Grid * g);
};


/*
 * Factory for producing single linear operators,
 * not left and right hand sides to problems
 */

template <class T>
class SingleOperatorFactory{
protected:
	map<Grid*,T*> operators;
	virtual void produce(Grid * g) = 0;
public:
	static void clearMap(map<Grid*,T*> & m){
		//cout << "Starting to clear" << endl;
		typename map<Grid*,T*>::iterator it;
		for(it = m.begin(); it != m.end(); it++){
			//cout << "Checking operator at location " << it->second << endl;
			if(it->second != 0){
				//cout << "Deleting operator at location " << it->second << endl;
				delete it->second;
				//cout << "Zeroing" << endl;
				it->second = 0;
			}
		}
		//cout << "Clearing map" << endl;
		m.clear();
	}
	virtual ~SingleOperatorFactory(){
		clearMap(operators);
	}

	virtual bool contains(Grid * g){
		return operators.count(g) > 0;
	}
	virtual T * get(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return operators[g];
	}
};

/*
 * Factory for producing LHS and RHS operators for a
 * linear problem.
 */
template <class T>
class OperatorFactory{
protected:

public:map<Grid*,T*> lhs, rhs;
protected:
	virtual void produce(Grid * g) = 0;
public:
	static void clearMap(map<Grid*,T*> & m){
		SingleOperatorFactory<T>::clearMap(m);
	}
	virtual ~OperatorFactory(){
		clearMap(lhs);
		clearMap(rhs);
	}
	virtual bool contains(Grid * g){
		return lhs.count(g) > 0 && lhs[g] != 0 && rhs.count(g) > 0 && rhs[g] != 0;
	}
	virtual T * getLHS(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return lhs[g];
	}
	virtual T * getRHS(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return rhs[g];
	}
};

/*template <class T>
class VariableOperatorFactory : public OperatorFactory<T>{
	OperatorFactory<T> * variablePart;
public:
	virtual ~VariableOperatorFactory(){
		delete variablePart;
	}
	virtual bool reinitializeVariablePart() = 0;
	virtual bool contains(Grid * g){
		return lhs.count(g) > 0 && variablePart->contains(g) > 0;
	}
	virtual T * getLHS(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return lhs[g];
	}
	virtual T * getRHS(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return rhs[g];
	}
};*/

class LaplacianFactory;

class LaplacianFactory : public SingleOperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
};

class UpwindFactory;

class UpwindFactory : public SingleOperatorFactory<CellToCellOperator>{
	double aX, aY;
	void produce(Grid * g);
public:
	void setA(double aX, double aY);
};

class FENEUpwindFactory;

class FENEUpwindFactory : public SingleOperatorFactory<CellToCellOperator>{
	double Qmax, H;

	double u11, u12, u21, u22;

	FaceCoefficientContainer springSpeedQ1, springSpeedQ2;

	FaceDoubleArray getSpringSpeedQ1(Grid * g);
	FaceDoubleArray getSpringSpeedQ2(Grid * g);
	FaceDoubleArray getFlowSpeedQ1(Grid * g);
	FaceDoubleArray getFlowSpeedQ2(Grid * g);

	void produce(Grid * g);

public:
	FENEUpwindFactory();
	void setH(double H);
	double getH();

	void setQmax(double Qmax);
	double getQmax();

	void setGradU(double u11, double u12, double u21, double u22);
};

class DiffusionBackwardEulerFactory;

class DiffusionBackwardEulerFactory : public OperatorFactory<CellToCellOperator>{
	double D, deltaT, t;

	void produce(Grid * g);
public:
	void setTime(double t);
	DiffusionBackwardEulerFactory(double D, double deltaT);
	/*DiffusionBackwardEulerFactory(double D, double deltaT){
		this->D = D;
		this->deltaT = deltaT;
	}*/
};

class AdvectionBackwardEulerFactory;

class AdvectionBackwardEulerFactory : public OperatorFactory<CellToCellOperator>{
	double aX, aY, deltaT;
	void produce(Grid * g);
public:
	void setA(double aX, double aY);
	void setDeltaT(double deltaT);
};

class AdvectionDiffusionBackwardEulerFactory;

class AdvectionDiffusionBackwardEulerFactory : public OperatorFactory<CellToCellOperator>{
	double aX, aY, D, deltaT;
	LaplacianFactory laplacian;
	UpwindFactory uw;
	void produce(Grid * g);
public:
	AdvectionDiffusionBackwardEulerFactory(double aX, double aY, double D, double deltaT);
	void setA(double aX, double aY);
	void setD(double D);
	void setDeltaT(double deltaT);
};

class FENEBackwardEulerFactory;

class FENEBackwardEulerFactory : public OperatorFactory<CellToCellOperator>{
	double D, deltaT;
	LaplacianFactory laplacian;
	FENEUpwindFactory uw;
	void produce(Grid * g);
public:
	FENEBackwardEulerFactory(double deltaT);
	void setDeltaT(double deltaT);
	void setD(double D);
	void setQmax(double Qmax);
	void setH(double H);
	void setGradU(double u11, double u12, double u21, double u22);
};

class FENESteadyFactory;

class FENESteadyFactory : public OperatorFactory<CellToCellOperator>{
	double D;
	LaplacianFactory laplacian;
	FENEUpwindFactory uw;
	void produce(Grid * g);
public:
	FENESteadyFactory();
	void setD(double D);
	void setQmax(double Qmax);
	void setH(double H);
	void setGradU(double u11, double u12, double u21, double u22);
};



/*
template <class T>
class OperatorFactory{
protected:
	map<Grid*,T*> operators;
	virtual void produce(Grid * g) = 0;
public:
	virtual ~OperatorFactory(){
		typename map<Grid*,T*>::iterator it;
		for(it = operators.begin(); it != operators.end(); it++){
			delete (*it).second;
		}
	}
	void remove(Grid * g){
		if(!contains(g)){
			return;
		}
		T * op = operators[g];
		operators.erase(g);
		delete op;
	}
	bool contains(Grid * g){
		return operators.count(g) > 0;
	}
	T * get(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return operators[g];
	}
};

template <class T, int levels = 2>
class TimeOperatorFactory{
protected:
	OperatorFactory<T> * op;
public:
	virtual ~TimeOperatorFactory(){}
	virtual T * get(Grid * g, int l) = 0;
};

template <class T>
class BackwardEulerFactory : public TimeOperatorFactory<T>{
public:
	BackwardEulerFactory(OperatorFactory<T> * op);
	T * get(Grid * g, int l);
};

template <class T>
class SplitOperatorFactory{
protected:
	map<Grid*,T*> lhsOperators, rhsOperators;
	virtual void produce(Grid * g) = 0;
public:
	virtual ~SplitOperatorFactory(){
		typename map<Grid*,T*>::iterator it;
		for(it = lhsOperators.begin(); it != lhsOperators.end(); it++){
			delete (*it).second;
		}
		for(it = rhsOperators.begin(); it != rhsOperators.end(); it++){
			delete (*it).second;
		}
	}
	void remove(Grid * g){
		if(!contains(g)){
			return;
		}
		T * lhsOp = lhsOperators[g];
		T * rhsOp = rhsOperators[g];
		lhsOperators.erase(g);
		rhsOperators.erase(g);
		delete lhsOp;
		delete rhsOp;
	}
	bool contains(Grid * g){
		return lhsOperators.count(g) > 0;
	}
	T * getLHS(Grid * g){
		//cout << "Getting LHS" << endl;
		if(!contains(g)){
 			produce(g);
		}
		return lhsOperators[g];
	}
	T * getRHS(Grid * g){
		if(!contains(g)){
			produce(g);
		}
		return rhsOperators[g];
	}
};

class LaplacianFactory;

class CrankNicholsonDiffusionFactory;

class UpwindFactory;

class UpwindFactory : public OperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
	double aX, aY;
	TinyVector<double,2> (* c)(Coord pos);
public:
	UpwindFactory(TinyVector<double,2> (* c)(Coord pos));
};

class CrankNicholsonDiffusionFactory : public SplitOperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
	double deltaT;
public:
	CrankNicholsonDiffusionFactory(double deltaT);
};

class LaplacianFactory : public OperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
};
*/

/*class LaplacianFactory : public OperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
};*/
}
}

#endif /* OPERATORFACTORY_H_ */
