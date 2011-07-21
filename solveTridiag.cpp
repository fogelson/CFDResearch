/*
 * mxSolveTridiag.cpp
 *
 *  Created on: Feb 25, 2011
 *      Author: fogelson
 */

#include <blitz/array.h>
#include <lapackpp/lapackpp.h>

using namespace la;
using namespace blitz;

// x = mxSolveTridiag(lower,center,upper,b)
int mainNot(){
	for(int j = 0; j < 5; j++){
		cout << j << endl;
		int N = 30;

		Array<double,1> lower(N);
		Array<double,1> center(N);
		Array<double,1> upper(N);
		Array<double,1> b(N);

		lower = 1;
		center = -2;
		upper = 1;
		b = 3;

		LaVectorDouble laLower(N-1);
		LaVectorDouble laUpper(N-1);
		LaVectorDouble laCenter(N);
		LaVectorDouble laB(N);
		laLower(0) = lower(0);
		laCenter(0) = center(0);
		laUpper(0) = upper(1);
		laB(0) = b(0);
		for(int n = 1; n <= N-2; n++){
			laLower(n) = lower(n);
			laCenter(n) = center(n);
			laUpper(n) = upper(n+1);
			laB(n) = b(n);
		}
		laLower(N-1) = lower(N-1);
		laCenter(N-1) = center(N-1);
		laB(N-1) = b(N-1);

		LaTridiagMatDouble L(laCenter, laLower, laUpper);

		LaTridiagFactDouble fact;
		LaGenMatDouble laX(N,1);
		LaTridiagMatFactorize(L,fact);
		LaLinearSolve(fact, laX, laB);

		Array<double,1> x(N);

		for(int n = 0; n < N; n++){
			x(n) = laX(n,0);
		}

	//	cout << laX;
		cout << x << endl;
	/*	cout << L << endl;

		cout << laLower << endl;
		cout << "---" << endl;
		cout << laCenter << endl;
		cout << "---" << endl;
		cout << laUpper << endl;
		cout << "---" << endl;
		cout << laB << endl;*/
	}
	return 0;
}
