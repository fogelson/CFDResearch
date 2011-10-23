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

typedef int FaceIndex;

typedef map<CellIndex, double> CellCoefficients;
typedef map<FaceIndex, double> FaceCoefficients;

/*
 * First index is the range ("to"), second index is the domain ("from").
 */
typedef map<FaceIndex, CellCoefficients> CellToFaceCoefficients;
typedef map<CellIndex, CellCoefficients> CellToCellCoefficients;
typedef map<CellIndex, FaceCoefficients> FaceToCellCoefficients;
typedef map<FaceIndex, FaceCoefficients> FaceToFaceCoefficients;

class Gradient;
class ConstantAdvectiveFlux;

class Divergence;

struct CellIndex{
	int i, j;
	CellIndex(int i, int j);
	bool operator== (const CellIndex & rhs) const;
	bool operator!= (const CellIndex & rhs) const;
	bool operator< (const CellIndex & rhs) const;
	bool operator> (const CellIndex & rhs) const;
	bool operator<= (const CellIndex & rhs) const;
	bool operator>= (const CellIndex & rhs) const;
};

class CellToFaceOperator{
public:
	Grid * g;
	CellToFaceCoefficients coefficients;
	FaceDoubleArray constantTerm;
public:
	FaceDoubleArray apply(CellDoubleArray & u);
	//FaceDoubleArray operator() (CellDoubleArray & u);
	FaceDoubleArray operator() (CellDoubleArray u);
	CellToFaceOperator & operator= (const CellToFaceOperator & rhs);
	CellToFaceOperator operator+ (CellToFaceOperator & B);
	CellToFaceOperator operator- (CellToFaceOperator & B);

};
CellToFaceOperator operator* (FaceDoubleArray & a, CellToFaceOperator & B);
CellToFaceOperator operator* (double a, CellToFaceOperator & B);

class ConstantAdvectiveFlux : public CellToFaceOperator{
	double aX, aY;
public:
	ConstantAdvectiveFlux(Grid * g, double aX, double aY);
};

class Gradient : public CellToFaceOperator{
	//void interpolateIrregularFace(double alpha, CellToFaceIndex ind1, CellToFaceIndex ind2, CellToFaceIndex ind3, CellToFaceIndex ind4);
	void interpolateIrregularFace(double alpha, CellIndex ind1, CellIndex ind2, CellIndex ind3, CellIndex ind4, FaceIndex faceIndex);
public:
	Gradient(Grid * g);
};

class FaceToCellOperator{
public:
	Grid * g;
	FaceToCellCoefficients coefficients;
	CellDoubleArray constantTerm;
public:
	CellDoubleArray apply(FaceDoubleArray & u);
	//CellDoubleArray operator() (FaceDoubleArray & u);
	CellDoubleArray operator() (FaceDoubleArray u);
	CellToCellOperator operator() (CellToFaceOperator & B);
};

class Divergence : public FaceToCellOperator{
public:
	Divergence(Grid * g);
};

class CellToCellOperator{
public:
	Grid * g;
	CellToCellCoefficients coefficients;
	CellDoubleArray constantTerm;
public:
	static CellToCellOperator getIdentity(Grid * g);
	CellToCellOperator & operator= (const CellToCellOperator & rhs);
	CellToCellOperator operator+ (CellToCellOperator & B);
	CellToCellOperator operator- (CellToCellOperator & B);
	CellDoubleArray apply(CellDoubleArray & u);
	CellDoubleArray operator() (CellDoubleArray & u);
};

CellToCellOperator operator* (CellDoubleArray & a, CellToCellOperator & B);
CellToCellOperator operator* (double a, CellToCellOperator & B);

/*class Laplacian : public CellToCellOperator{
public:
	Laplacian(Grid * g);
};*/

}
}

#endif /* OPERATORS_H_ */
