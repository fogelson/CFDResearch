/*
 * MexTools.h
 *
 *  Created on: Oct 13, 2011
 *      Author: fogelson
 */

#ifndef MEXTOOLS_H_
#define MEXTOOLS_H_

#include <string>
#include "BlitzMatlab.h"
#include "../Geo/Geometry.h"
//#include "/Applications/MATLAB_R2010a.app/extern/include/mex.h"
//#include "/Applications/MATLAB_R2010a.app/extern/include/matrix.h"
#include "mex.h"
#include "matrix.h"

namespace CFD{
using namespace OOGeometry;
using namespace std;
using namespace blitzmatlab;
namespace OOMexTools{

class MexPlotTool;

class MexPlotTool{
	int movieFrame;
public:
	void newFigure();
	void holdOn();
	void holdOff();
	void graphCellDoubleArray(CellDoubleArray & u, Grid * g, string graphCommand);
	void graphFaceDoubleArray(FaceDoubleArray & u, Grid * g);
	void drawGrid(Grid * g);
	void plotFaceDoubleArrayXLine(FaceDoubleArray & u, Grid * g, int i, Direction d);
	void drawNow();
	void graphCellCentroidData(CellDoubleArray & u, Grid * g);
	void colorbar();
	void initializeMovie();
	void captureFrame();
	void playMovie();
	void title(string t);
};

}
}

#endif /* MEXTOOLS_H_ */
