/*
 * Operators.h
 *
 *  Created on: Oct 12, 2011
 *      Author: fogelson
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <map>
#include <blitz/tinyvec-et.h>
#include "../Geo/Geometry.h"

using namespace blitz;
using namespace std;

namespace CFD{
using namespace OOGeometry;
namespace OOOps{

class CellToFaceOperator;

struct CellToFaceIndex;
struct CellToFaceIndexCompare;

typedef map<CellToFaceIndex,double,CellToFaceIndexCompare> CellToFaceCoefficients;

class Gradient;

class FaceToCellOperator;

class Divergence;

class CellToCellOperator;

struct CellToCellIndex;
struct CellToCellIndexCompare;

typedef map<CellToCellIndex,double,CellToCellIndexCompare> CellToCellCoefficients;

class FaceToFaceOperator;

struct CellToFaceIndex{
	CellToFaceIndex(int i, int j, int faceIndex);
	int i, j, faceIndex;
};
struct CellToFaceIndexCompare{
	bool operator() (const CellToFaceIndex & lhs, const CellToFaceIndex & rhs) const;
};


class CellToFaceOperator{
protected:
	Grid * g;
	CellToFaceCoefficients coefficients;
public:
	FaceDoubleArray apply(CellDoubleArray & u);
	FaceDoubleArray operator() (CellDoubleArray & u);
};

class Gradient : public CellToFaceOperator{
	void interpolateIrregularFace(double alpha, CellToFaceIndex ind1, CellToFaceIndex ind2, CellToFaceIndex ind3, CellToFaceIndex ind4);
public:
	Gradient(Grid * g);
};

class FaceToCellOperator{
protected:
	Grid * g;
	CellToFaceCoefficients coefficients;
public:
	CellDoubleArray apply(FaceDoubleArray & u);
	CellDoubleArray operator() (FaceDoubleArray & u);
};

class Divergence : public FaceToCellOperator{
public:
	Divergence(Grid * g);
};

struct CellToCellIndex{
	int i1, i2, j1, j2;
	CellToCellIndex(int i1, int j1, int i2, int j2);
};

struct CellToCellIndexCompare{
	bool operator() (const CellToCellIndex & lhs, const CellToCellIndex & rhs) const;
};

class CellToCellOperator{
protected:
	Grid * g;
	CellToCellCoefficients coefficients;
public:
	CellDoubleArray apply(CellDoubleArray & u);
	CellDoubleArray operator() (CellDoubleArray & u);
};

}
}

#endif /* OPERATORS_H_ */
