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
class FaceToCellOperator;
class FaceToFaceOperator;
class CellToCellOperator;

struct CellIndex;
struct CellIndexCompare;

typedef int FaceIndex;

typedef map<CellIndex, double, CellIndexCompare> CellCoefficients;
typedef map<FaceIndex, double> FaceCoefficients;

typedef map<CellIndex, FaceCoefficients, CellIndexCompare> NewCellToFaceCoefficients;
typedef map<CellIndex, CellCoefficients, CellIndexCompare> NewCellToCellCoefficients;
typedef map<FaceIndex, CellCoefficients > NewFaceToCellCoefficients;
typedef map<FaceIndex, FaceCoefficients > NewFaceToFaceCoefficients;

struct CellToFaceIndex;
struct CellToFaceIndexCompare;
typedef map<CellToFaceIndex,double,CellToFaceIndexCompare> CellToFaceCoefficients;

class Gradient;
class ConstantAdvectiveFlux;

class Divergence;

struct CellToCellIndex;
struct CellToCellIndexCompare;

typedef map<CellToCellIndex,double,CellToCellIndexCompare> CellToCellCoefficients;

struct CellIndex{
	int i, j;
	CellIndex(int i, int j);
};
struct CellIndexCompare{
	bool operator() (const CellIndex & lhs, const CellIndex & rhs) const;
};

struct CellToFaceIndex{
	CellToFaceIndex(int i, int j, int faceIndex);
	int i, j, faceIndex;
};
struct CellToFaceIndexCompare{
	bool operator() (const CellToFaceIndex & lhs, const CellToFaceIndex & rhs) const;
};


class CellToFaceOperator{
public:
	Grid * g;
	CellToFaceCoefficients coefficients;
	NewCellToFaceCoefficients newCoefficients;
	FaceDoubleArray constantTerm;
	int getCols();
	int getRows();
	int getColIndex(CellToFaceIndex ind);
	int getRowIndex(CellToFaceIndex ind);
public:
	FaceDoubleArray apply(CellDoubleArray & u);
	//FaceDoubleArray operator() (CellDoubleArray & u);
	FaceDoubleArray operator() (CellDoubleArray u);
	Array<double,2> getMatrix();
	Array<double,1> getVector();
};

class ConstantAdvectiveFlux : public CellToFaceOperator{
	double aX, aY;
public:
	ConstantAdvectiveFlux(Grid * g, double aX, double aY);
};

class Gradient : public CellToFaceOperator{
	void interpolateIrregularFace(double alpha, CellToFaceIndex ind1, CellToFaceIndex ind2, CellToFaceIndex ind3, CellToFaceIndex ind4);
	void interpolateIrregularFaceNew(double alpha, CellIndex ind1, CellIndex ind2, CellIndex ind3, CellIndex ind4, FaceIndex faceIndex);
public:
	Gradient(Grid * g);
};

class FaceToCellOperator{
public:
	Grid * g;
	CellToFaceCoefficients coefficients;
	CellDoubleArray constantTerm;
	int getCols();
	int getRows();
	int getColIndex(CellToFaceIndex ind);
	int getRowIndex(CellToFaceIndex ind);
public:
	CellDoubleArray apply(FaceDoubleArray & u);
	//CellDoubleArray operator() (FaceDoubleArray & u);
	CellDoubleArray operator() (FaceDoubleArray u);
	CellToCellOperator operator() (CellToFaceOperator & B);
	Array<double,2> getMatrix();
	Array<double,1> getVector();
};

class Divergence : public FaceToCellOperator{
public:
	Divergence(Grid * g);
};

struct CellToCellIndex{
	int iFrom, jFrom, iTo, jTo;
	CellToCellIndex(int iFrom, int jFrom, int iTo, int jTo);
};

struct CellToCellIndexCompare{
	bool operator() (const CellToCellIndex & lhs, const CellToCellIndex & rhs) const;
};

class CellToCellOperator{
public:
	Grid * g;
	CellToCellCoefficients coefficients;
	CellDoubleArray constantTerm;
public:
	CellToCellOperator & operator= (const CellToCellOperator & rhs);
	CellDoubleArray apply(CellDoubleArray & u);
	CellDoubleArray operator() (CellDoubleArray & u);
};

class Laplacian : public CellToCellOperator{
	double getDirichlet(Coord c);
	double getDirichlet(double x, double y);
	void applyIrregularGradient(double & c1, double & c2, double & c3, double & c4, double a, double v);
public:
	Laplacian(Grid * g);
};

}
}

#endif /* OPERATORS_H_ */
