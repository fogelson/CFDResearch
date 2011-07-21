#include "mex.h"
#include "matrix.h"
#include "BlitzMatlab.h"
#include "Geometry/Geometry.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace blitzmatlab;
using namespace CFD;
using namespace Geometry;

/* nlhs:	number of outputs to Matlab function
 * plhs:	array of pointers to the memory locations of the outputs
 * nrhs:	number of inputs to Matlab function
 * prhs:	array of pointers to the memory locations of the inputs
 */

/* The following command should be invoked from MATLAB:
 * [B] = mxFileTemplate(c,A)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Desired number of inputs and outputs
	int inputs = 2;
	int outputs = 4;

	/* Check for proper number of arguments. */
	if(nrhs != inputs){
		mexErrMsgTxt("Wrong number of inputs.");
	}
	else if(nlhs > outputs){
		mexErrMsgTxt("Too many outputs.");
	}

	string runname;
	char * input_buf;
	mwSize buflen;

	/* input must be a string */
	if ( mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("Input must be a string.");

	/* input must be a row vector */
	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("Input must be a row vector.");

	/* get the length of the input string */
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

	/* copy the string data from prhs[0] into a C string input_ buf.    */
	input_buf = mxArrayToString(prhs[0]);

	if(input_buf == NULL)
		mexErrMsgTxt("Could not convert input to string.");

	runname = input_buf;

	int n = getMxInt(prhs[1]);

	string configfilename = runname + "_config.cfd";

	ifstream infile;
	infile.open(configfilename.c_str(), ios::in | ios::binary);

	int saveFrequency,totalSaves;
	double h, offset, deltaT, T, D, alpha, H, Q0, u11, u12, u21, u22;

	infile.read((char *)(&h), sizeof(h));
	infile.read((char *)(&deltaT), sizeof(deltaT));
	infile.read((char *)(&T), sizeof(T));
	infile.read((char *)(&saveFrequency), sizeof(saveFrequency));
	infile.read((char *)(&D), sizeof(D));
	infile.read((char *)(&alpha), sizeof(alpha));
	infile.read((char *)(&u11), sizeof(u11));
	infile.read((char *)(&u12), sizeof(u12));
	infile.read((char *)(&u21), sizeof(u21));
	infile.read((char *)(&u22), sizeof(u22));
	infile.read((char *)(&H), sizeof(H));
	infile.read((char *)(&Q0), sizeof(Q0));
	infile.close();

	offset = h/2;

	cout << h << endl;
	cout << deltaT << endl;
	cout << T << endl;
	cout << saveFrequency << endl;
	cout << D << endl;
	cout << alpha << endl;
	cout << u11 << endl;
	cout << u12 << endl;
	cout << u21 << endl;
	cout << u22 << endl;
	cout << H << endl;
	cout << Q0 << endl;


	Grid * grid = new Circle(h,Q0,offset);
	CellDoubleArray u = grid->makeCellDoubleArray();
	double t;

	stringstream ss;
	ss << runname << "_" << n << ".cfddata";
	infile.open(ss.str().c_str(), ios::in | ios::binary);
	infile.read((char *)(&t), sizeof(t));
	infile.read((char *)(u.data()), ((long)u.size())*sizeof(double));
	infile.close();
	u = where(grid->cellTypes == COVERED, mxGetNaN(), u);

	CoordArray cellCentroids = grid->cellCentroids;
	CellDoubleArray xCellCentroids = cellCentroids.extractComponent(double(),0,2);
	CellDoubleArray yCellCentroids = cellCentroids.extractComponent(double(),1,2);

	// Set output pointers to the desired outputs
	plhs[0] = setMxArray(u);
	plhs[1] = setMxDouble(t);
	plhs[2] = setMxArray(xCellCentroids);
	plhs[3] = setMxArray(yCellCentroids);

	delete grid;
}
