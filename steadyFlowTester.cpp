/*
 * steadyFlowTester.cpp
 *
 *  Created on: May 20, 2011
 *      Author: fogelson
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Solvers/StenciledAdvectionDiffusionSolver.h"
#include "Multigrid/Smoothers.h"

using namespace CFD;
using namespace Solvers;
using namespace Geometry;
using namespace Multigrid;
using namespace std;

int main(){
	string runname;

	double h, offset, deltaT, T, D, alpha, H, Q0, u11, u12, u21, u22;

	int saveFrequency;

	cout << "Run name: ";
	getline(cin,runname);
	string configfilename = runname + "_config.cfd";

	cout << "Config data will be stored in " << configfilename << endl << endl;

	cout << "h = ";
	cin >> h;
	offset = h/2;

	cout << "deltaT = ";
	cin >> deltaT;

	cout << "T_final = ";
	cin >> T;

	cout << "Save solution every n timesteps. n = ";
	cin >> saveFrequency;

	cout << "D = ";
	cin >> D;

	cout << "alpha = ";
	cin >> alpha;

	cout << "u11 = ";
	cin >> u11;

	cout << "u12 = ";
	cin >> u12;

	cout << "u21 = ";
	cin >> u21;

	cout << "u22 = ";
	cin >> u22;

	H = 1;
	Q0 = 1;

	ofstream outfile;
	outfile.open(configfilename.c_str(), ios::out | ios::binary | ios::trunc);

	outfile.write((char *)(&h), sizeof(h));
	outfile.write((char *)(&deltaT), sizeof(deltaT));
	outfile.write((char *)(&T), sizeof(T));
	outfile.write((char *)(&saveFrequency), sizeof(saveFrequency));
	outfile.write((char *)(&D), sizeof(D));
	outfile.write((char *)(&alpha), sizeof(D));
	outfile.write((char *)(&u11), sizeof(D));
	outfile.write((char *)(&u12), sizeof(D));
	outfile.write((char *)(&u21), sizeof(D));
	outfile.write((char *)(&u22), sizeof(D));
	outfile.write((char *)(&H), sizeof(D));
	outfile.write((char *)(&Q0), sizeof(D));

	outfile.close();

	Array<double,2> gradU;
	Range trivialRange(1,2);
	gradU.resize(trivialRange,trivialRange);
	gradU(1,1) = u11;
	gradU(1,2) = u12;
	gradU(2,1) = u21;
	gradU(2,2) = u22;

	gradU = alpha*gradU;

	Grid * grid = new Circle(h,Q0,h/2);
	CellDoubleArray u = grid->makeCellDoubleArray();
	CellDoubleArray uNew = grid->makeCellDoubleArray();

	u = where(grid->cellTypes == COVERED, 0, 1);

	StenciledFourPointGS smoother;
	BilinearInterpolator bi;
	VolumeWeightedRestrictor vw;
	StenciledMultigridSolver multigrid(&smoother,&bi,&vw);

	FENEStencil stencil(grid,deltaT,D,H,Q0);
	stencil.setGradU(gradU);

	double t = 0;
	int n = 0;
	while(t <= T){
		if(n % saveFrequency == 0){
			stringstream ss;
			ss << runname << "_" << n << ".cfddata";
			outfile.open(ss.str().c_str(), ios::out | ios::binary | ios::trunc);
			outfile.write((char *)(&t), sizeof(t));
			outfile.write((char *)(u.data()), ((long)u.size())*sizeof(double));
			outfile.close();
		}

		n++;
		t = ((double)n)*deltaT;

		uNew = multigrid.solve(u,u,&stencil,2,2,1);
		u = uNew;
	}

	delete grid;

	return 0;
}
