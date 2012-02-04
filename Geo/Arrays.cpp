/*
 * Arrays.cpp
 *
 *  Created on: Oct 7, 2011
 *      Author: fogelson
 */

#include "Geometry.h"
#include <iostream>

using namespace std;

namespace CFD{
namespace OOGeometry{

CellDoubleArray Arrays::makeCellDoubleArray(Grid * g){
	CellDoubleArray out;
	out.resize(g->xRange,g->yRange);
	out = 0;
	return out;
}

FaceDoubleArray Arrays::makeFaceDoubleArray(Grid * g){
	FaceDoubleArray out;
	out.resize(g->faces.size());
	out = 0;
	return out;
}

VertexDoubleArray Arrays::makeVertexDoubleArray(Grid * g){
	VertexDoubleArray out;
	out.resize(g->vertices.size());
	out = 0;
	return out;
}

}
}
