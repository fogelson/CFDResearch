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
			//cout << "Need to produce LHS" << endl;
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

class CrankNicholsonDiffusionFactory : public SplitOperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
	double deltaT;
public:
	CrankNicholsonDiffusionFactory(double deltaT);
};

class LaplacianFactory : public OperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
};


/*class LaplacianFactory : public OperatorFactory<CellToCellOperator>{
	void produce(Grid * g);
};*/
}
}

#endif /* OPERATORFACTORY_H_ */
