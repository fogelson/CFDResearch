/*
 * testParameterFile.cpp
 *
 *  Created on: Jul 22, 2011
 *      Author: fogelson
 */

#include "parameterFile.h"

using namespace CFD;

int mainNot(){
	ParameterFileReader pfr;
	pfr.setFilename("paramTest");
	pfr.readFile();
	return 0;
}
