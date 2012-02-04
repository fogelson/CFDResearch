/*
 * FENEUpwindFactory.cpp
 *
 *  Created on: Dec 24, 2011
 *      Author: fogelson
 */


#include "Operators.h"
#include "OperatorFactory.h"
#include <iostream>

using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

FaceDoubleArray FENEUpwindFactory::getSpringSpeedQ1(Grid * g){
	//if(!springSpeedQ1.contains(g)){
		FaceDoubleArray out = g->makeFaceDoubleArray();
		FaceDoubleArray q1 = g->getFaceX();
		FaceDoubleArray q2 = g->getFaceY();
		FaceDoubleArray Q = g->makeFaceDoubleArray();
		Q = sqrt(pow2(q1) + pow2(q2));
		out = -H*q1/(1 - pow2(Q)/pow2(Qmax));
		//springSpeedQ1(g) = out;
	//}
	return out;
}

FaceDoubleArray FENEUpwindFactory::getSpringSpeedQ2(Grid * g){
	//if(!springSpeedQ2.contains(g)){
		FaceDoubleArray out = g->makeFaceDoubleArray();
		FaceDoubleArray q1 = g->getFaceX();
		FaceDoubleArray q2 = g->getFaceY();
		FaceDoubleArray Q = g->makeFaceDoubleArray();
		Q = sqrt(pow2(q1) + pow2(q2));
		out = -H*q2/(1 - pow2(Q)/pow2(Qmax));
		//springSpeedQ1(g) = out;
	//}
	return out;
}

FaceDoubleArray FENEUpwindFactory::getFlowSpeedQ1(Grid * g){
	FaceDoubleArray out = g->makeFaceDoubleArray();
	FaceDoubleArray q1 = g->getFaceX();
	FaceDoubleArray q2 = g->getFaceY();
	out = u11*q1 + u12*q2;
	return out;
}

FaceDoubleArray FENEUpwindFactory::getFlowSpeedQ2(Grid * g){
	FaceDoubleArray out = g->makeFaceDoubleArray();
	FaceDoubleArray q1 = g->getFaceX();
	FaceDoubleArray q2 = g->getFaceY();
	out = u21*q1 + u22*q2;
	return out;
}


FENEUpwindFactory::FENEUpwindFactory(){
	H = 1;
	Qmax = 1;
	u11 = 0;
	u21 = 0;
	u12 = 0;
	u22 = 0;
}
/*
if(this->aX == aX && this->aY == aY)
	return;
if(!operators.empty()){
	clearMap(operators);
}
*/
void FENEUpwindFactory::setH(double H){
	if(this->H == H)
		return;
	if(!operators.empty()){
		clearMap(operators);
		springSpeedQ1.clear();
		springSpeedQ2.clear();
	}
	this->H = H;
}
double FENEUpwindFactory::getH(){
	return H;
}

void FENEUpwindFactory::setQmax(double Qmax){
	if(this->Qmax == Qmax)
		return;
	if(!operators.empty()){
		clearMap(operators);
		springSpeedQ1.clear();
		springSpeedQ2.clear();
	}
	this->Qmax = Qmax;
}
double FENEUpwindFactory::getQmax(){
	return Qmax;
}

void FENEUpwindFactory::setGradU(double u11, double u12, double u21, double u22){
	if(this->u11 == u11 && this->u12 == u12 && this->u21 == u21 && this->u22 == u22)
		return;
	if(!operators.empty()){
		clearMap(operators);
	}
	this->u11 = u11;
	this->u12 = u12;
	this->u21 = u21;
	this->u22 = u22;
}

void FENEUpwindFactory::produce(Grid * g){
	if(contains(g)){
		return;
	}
	CellToFaceOperator F;
	F.setGrid(g);
	vector<Face*>::iterator it;

	if(!springSpeedQ1.contains(g)){
		springSpeedQ1(g) = getSpringSpeedQ1(g);
	}
	if(!springSpeedQ2.contains(g)){
		springSpeedQ2(g) = getSpringSpeedQ2(g);
	}

	FaceDoubleArray flowSpeedQ1 = g->makeFaceDoubleArray();
	FaceDoubleArray flowSpeedQ2 = g->makeFaceDoubleArray();

	flowSpeedQ1 = getFlowSpeedQ1(g);
	flowSpeedQ2 = getFlowSpeedQ2(g);

	FaceDoubleArray speedQ1 = g->makeFaceDoubleArray(), speedQ2 = g->makeFaceDoubleArray();

	speedQ1 = springSpeedQ1(g) + flowSpeedQ1;
	speedQ2 = springSpeedQ2(g) + flowSpeedQ2;

	// THIS IS CAUSING CRASH/MEMORY PROBLEMS
	/*********** THIS LINE RIGHT ABOVE HERE. YESSIR. THAT'S THE ONE TO WATCH.
	 * THAT'S THE PROBLEM, RIGHT HERE. DEAL WITH THAT, AND I BETCHA YOU'LL BE GOOD,
	 * I TELL YOU WHHAT.
	 * ********
	 */


	/*
	if(springSpeedQ1Map.count(g) > 0){
		springQ1 = springSpeedQ1Map[g];
	}
	else{
		springSpeedQ1Map[g] = springSpeedQ1(g);
		springQ1 = springSpeedQ1Map[g];
	}

	if(springSpeedQ2Map.count(g) > 0){
		springQ2 = springSpeedQ2Map[g];
	}
	else{
		springSpeedQ2Map[g] = springSpeedQ2(g);
		springQ2 = springSpeedQ2Map[g];
	}

	flowQ1 = flowSpeedQ1(g);
	flowQ2 = flowSpeedQ2(g);*/

	/*FaceDoubleArray aQ1 = g->makeFaceDoubleArray(), aQ2 = g->makeFaceDoubleArray();
	aQ1 = springQ1 + flowQ1;
	aQ2 = springQ2 + flowQ2;*/

	for(it = g->faces.begin(); it != g->faces.end(); it++){
		if((*it)->isInterior()){
			double a = 0;
			if((*it)->isNS()){
				a = speedQ2((*it)->getIndex());
			}
			if((*it)->isEW()){
				a = speedQ1((*it)->getIndex());
			}
			Cell * upwind;
			if(a > 0){
				upwind = (*it)->getFrom();
			}
			else{
				upwind = (*it)->getTo();
			}
			FaceIndex faceIndex = (*it)->getIndex();
			CellIndex cellIndex(upwind->getI(), upwind->getJ());
			F.coefficients[faceIndex][cellIndex] = -a;
		}
	}
	Divergence div(g);

	operators[g] = new CellToCellOperator(div(F));
}


/*void UpwindFactory::setA(double aX, double aY){
	if(this->aX == aX && this->aY == aY)
		return;
	if(!operators.empty()){
		clearMap(operators);
	}
	this->aX = aX;
	this->aY = aY;
}*/

}
}
