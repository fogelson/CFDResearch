/*
 * OperatorFactory.cpp
 *
 *  Created on: Oct 25, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

template <class T>
T & OperatorFactory<T>::get(Grid * g){
	cout << "balls" << endl;
	T out;
	return out;
}

/*template <typename T>
OperatorFactory<T>::~OperatorFactory(){
	clear();
}

void OperatorFactory<CellToCellOperator>::clear(){
	map<Grid*,CellToCellOperator*>::iterator it;
	for(it = created.begin(); it != created.end(); it++){
		delete (*it).second;
	}
	created.clear();
}

void OperatorFactory<CellToCellOperator>::remove(Grid * g){
	if(created.count(g) > 0){
		CellToCellOperator * op = created[g];
		created.erase(g);
		delete op;
	}
}

bool OperatorFactory<CellToCellOperator>::contains(Grid * g){
	return created.count(g) > 0;
}

CellToCellOperator & OperatorFactory<CellToCellOperator>::get(Grid * g){
	if(!contains(g)){
		produce(g);
	}
	return &(created[g]);
}*/


}
}
