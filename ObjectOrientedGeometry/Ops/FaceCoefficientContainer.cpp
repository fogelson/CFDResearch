/*
 * FaceCoefficientContainer.cpp
 *
 *  Created on: Dec 26, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

FaceCoefficientContainer::~FaceCoefficientContainer(){
	clear();
}
void FaceCoefficientContainer::clear(){
	map<Grid*,FaceDoubleArray*>::iterator it;
	for(it = coefficients.begin(); it != coefficients.end(); it++){
		if((*it).second != 0){
			delete (*it).second;
			(*it).second = 0;
		}
	}
	coefficients.clear();
}
FaceDoubleArray & FaceCoefficientContainer::get(Grid * g){
	if(!contains(g) || coefficients[g] == 0){
		coefficients[g] = new FaceDoubleArray(g->makeFaceDoubleArray());
	}
	return (*coefficients[g]);
}
FaceDoubleArray & FaceCoefficientContainer::operator() (Grid * g){
	return get(g);
}
bool FaceCoefficientContainer::contains(Grid * g){
	return coefficients.count(g) != 0;
}


}
}
