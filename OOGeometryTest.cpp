/*
 * OOGeometryTest.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

/*
 * GridElement.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#include "Geo/Geometry.h"
#include <vector>
#include <iostream>

using namespace std;

using namespace CFD::OOGeometry;

int main(){
	double r = 1, h = .2, offset = h/2;
	Circle * g = new Circle(h,r,offset);

	/*	FaceDoubleArray fda = g->makeFaceDoubleArray();
	cout << fda << endl;
	fda.resize(g->faces.size());
	fda = 0;
	cout << fda << endl;*/
	CellDoubleArray cda = g->makeCellDoubleArray();
	cout << cda << endl;

	FaceDoubleArray fda = g->makeFaceDoubleArray();
	cout << fda << endl;

	VertexDoubleArray vda = g->makeVertexDoubleArray();
	cout << vda << endl;

	delete g;

	return 0;
}

