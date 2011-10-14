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
#include "/opt/pkg/mathworks/matlab-2011a/extern/include/mex.h"
#include "/opt/pkg/mathworks/matlab-2011a/extern/include/matrix.h"

namespace CFD{
using namespace OOGeometry;
using namespace std;
using namespace blitzmatlab;
namespace OOMexTools{

class MexPlotTool;

class MexPlotTool{
public:
	void newFigure();
	void holdOn();
	void holdOff();
	void graphCellDoubleArray(CellDoubleArray & u, Grid * g, string graphCommand);
	void graphFaceDoubleArray(FaceDoubleArray & u, Grid * g);
	void drawGrid(Grid * g);
};

}
}

#endif /* MEXTOOLS_H_ */
