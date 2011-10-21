/*
 * BlitzMatlab.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: fogelson
 */

#include "BlitzMatlab.h"

namespace blitzmatlab{
	mxArray* setMxDouble(double d){
        double *mxPointer;

        mxArray *mx = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);

        mxPointer = mxGetPr(mx);
        mxPointer[0] = d;

        return mx;
	}
	double getMxDouble(const mxArray *mx){
		double *mxPointer, *outPointer, out;
		mxPointer = mxGetPr(mx);
		outPointer = &out;
		outPointer[0] = mxPointer[0];
		return out;
	}
    mxArray* setMxInt(int i){
        return setMxDouble((double)i);
    }
	int getMxInt(const mxArray *mx){
		double d = getMxDouble(mx);
		int i = (int)d;
		if(i != d){
			mexErrMsgTxt("Input must be an integer.");
		}
		return (int)getMxDouble(mx);
	}

	// Create a row-majored, zero-indexed Blitz++ 2D array from MATLAB data
	Array<double,2> getMxArray(const mxArray *mx){
		mwSize rows, cols;
		int m, n;
		double *mxPointer, *blitzPointer;

		rows = mxGetM(mx);
		cols = mxGetN(mx);
		m = (int) rows;
		n = (int) cols;

		Array<double,2> columnMajor(m,n,fortranArray);

		mxPointer = mxGetPr(mx);
		blitzPointer = columnMajor.data();

		for(int i = 0; i < columnMajor.size(); i++){
			blitzPointer[i] = mxPointer[i];
		}

		Array<double,2> rowMajor(m,n);
		rowMajor = columnMajor;
		return rowMajor;
	}
	// Create a pointer to a MATLAB mxArray from row-majored, zero-index Blitz++ 2D array data
	mxArray* setMxArray(Array<double,2> rowMajor){
		mwSize rows, cols;
		int m, n;
		double *mxPointer, *blitzPointer;

		m = rowMajor.length(0);
		n = rowMajor.length(1);
		rows = (mwSize) m;
		cols = (mwSize) n;

		Array<double,2> columnMajor(m,n,fortranArray);
		columnMajor = rowMajor;

		mxArray *mx = mxCreateNumericMatrix(rows,cols,mxDOUBLE_CLASS,mxREAL);

		mxPointer = mxGetPr(mx);
		blitzPointer = columnMajor.data();

		for(int i = 0; i < columnMajor.size(); i++){
			mxPointer[i] = blitzPointer[i];
		}
		return mx;
	}

	mxArray* setMxVector(Array<double,1> rowMajor){
		Array<double,2> asArray(rowMajor.size(),1);
		asArray(Range::all(),0) = rowMajor(Range::all());
		return setMxArray(asArray);/*

		mwSize rows, cols;
		int m, n;
		double *mxPointer, *blitzPointer;

		m = rowMajor.length(0);
		n = 1;
		rows = (mwSize) m;
		cols = (mwSize) n;

		Array<double,1> columnMajor(m,fortranArray);

		mxArray *mx = mxCreateNumericMatrix(rows,cols,mxDOUBLE_CLASS,mxREAL);

		mxPointer = mxGetPr(mx);
		blitzPointer = columnMajor.data();

		for(int i = 0; i < columnMajor.size(); i++){
			mxPointer[i] = blitzPointer[i];
		}
		return mx;*/
	}

	Array<double,1> getMxVector(const mxArray *mx){
		Array<double,2> asArray = getMxArray(mx);
		Array<double,1> asVector(asArray.length(0));
		asVector(Range::all()) = asArray(Range::all(),0);
		return asVector;
		/*mwSize rows, cols;
		int m, n;
		double *mxPointer, *blitzPointer;
		rows = mxGetM(mx);
		cols = mxGetN(mx);
		m = (int) rows;
		n = (int) cols;
		Array<double,1> columnMajor(m,fortranArray);

		mxPointer = mxGetPr(mx);
		blitzPointer = columnMajor.data();

		for(int i = 0; i < columnMajor.size(); i++){
			blitzPointer[i] = mxPointer[i];
		}

		Array<double,1> rowMajor(m,n);
		rowMajor = columnMajor;
		return rowMajor;*/
	}

/*	GridScalar restrictToGrid(GridScalar u, CoordTypeArray types){
		Range iRange(u.lbound(0),u.ubound(0)), jRange(u.lbound(1),u.ubound(1));
		GridScalar out(iRange,jRange);
		out = u;
		for(int i = out.lbound(0); i <= out.ubound(0); i++){
			for(int j = out.lbound(1); j <= out.ubound(1); j++){
				if(types(i,j) == EXTERIOR){
					//out(i,j) = mxGetNaN();
				}
			}
		}
	}*/

}
