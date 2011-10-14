/*
 * BlitzMatlab.h
 *
 *  Created on: Nov 22, 2010
 *      Author: fogelson
 */

// Helper file for functions for converting to and from MATLAB array pointers
// and the Blitz++ array classes used in the C++ code.

#ifndef BLITZMATLAB_H_
#define BLITZMATLAB_H_

#include "mex.h"
#include "matrix.h"
#include <blitz/array.h>

using namespace blitz;

namespace blitzmatlab{
	mxArray* setMxDouble(double d);

	double getMxDouble(const mxArray *mx);

    mxArray* setMxInt(int i);

	int getMxInt(const mxArray *mx);

	// Create a row-majored, zero-indexed Blitz++ 2D array from MATLAB data
	Array<double,2> getMxArray(const mxArray *mx);

	// Create a pointer to a MATLAB mxArray from row-majored, zero-index Blitz++ 2D array data
	mxArray* setMxArray(Array<double,2> rowMajor);

	mxArray* setMxVector(Array<double,1> rowMajor);

	Array<double,1> getMxVector(const mxArray *mx);
}

#endif /* BLITZMATLAB_H_ */
