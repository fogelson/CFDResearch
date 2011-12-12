/*
 * Geometry.h
 *
 *  Created on: Oct 3, 2011
 *      Author: fogelson
 */

#ifndef OGEOMETRY_H_
#define OGEOMETRY_H_

#include <blitz/array.h>
#include <list>
#include <vector>
#include <blitz/tinyvec-et.h>

using namespace blitz;
using namespace std;

namespace CFD{
namespace OOGeometry{

class GridElement;
class Cell;
class Face;
class Vertex;
class Grid;
class Circle;

typedef Array<double,2> CellDoubleArray;
typedef Array<double,1> FaceDoubleArray;
typedef Array<double,1> VertexDoubleArray;

typedef TinyVector<double,2> Coord;

typedef int Type;
const Type REGULAR = 0;
const Type IRREGULAR = 1;
const Type COVERED = 2;

typedef int Direction;
const Direction N = 0;
const Direction E = 1;
const Direction S = 2;
const Direction W = 3;
const Direction B = 4;
const Direction NE = 4;
const Direction SE = 5;
const Direction NW = 6;
const Direction SW = 7;
const Direction C = 8;

class GridElement{
	int index;
	Type type;
public:
	void setIndex(int i);
	int getIndex();

	Type getType();
	bool isRegular();
	bool isIrregular();
	bool isCovered();
	bool isUncovered();
	void setType(Type type);
};

class Cell : public GridElement{
	int i, j;
	bool overrideVolume;
public:
	bool upToDate;

	double volume;

	vector<Vertex*> vertices;

	Grid * g;
	int numberOfVertices;
	int numberOfFaces;

	TinyVector<Face*,5> faces;

	Coord center, centroid;

	void update();

public:
	Cell();
	double getVolume();
	void setVolume(double volume);
	double getVolumeFraction();
	int getNumberOfVertices();
	int getNumberOfFaces();

	void setI(int i);
	void setJ(int j);
	int getI();
	int getJ();

	vector<Vertex*> getVertices();
	vector<Face*> getFaces();

	Coord getCenter();
	void setCenter(Coord c);

	Coord getCentroid();
	void setCentroid(Coord c);

	Face * operator() (Direction d);
	void setFace(Face * face, Direction d);
	bool hasFace(Direction d);
	Face * getFace(Direction d);

	Face * createBoundary();
};

class Face : public GridElement{
public:
	Grid * g;
	Cell * fromCell, * toCell;
	Vertex * vertexA, * vertexB;
	Cell * cell1, *cell2;
	Coord centroid;
	TinyVector<double,2> normal;
	bool upToDate;
	double area;

	bool isBoundaryFace;

	void update();
public:
	Face();
	Vertex * getVertexA();
	Vertex * getVertexB();
	Vertex * getA();
	Vertex * getB();

	void setFrom(Cell * fromCell);
	void setTo(Cell * toCell);
	Cell * getFrom();
	Cell * getTo();

	void setVertexA(Vertex * vertexA);
	void setVertexB(Vertex * vertexB);
	void setA(Vertex * A);
	void setB(Vertex * B);

	double getArea();
	double getAreaFraction();

	Coord getCentroid();

	void setNormal(TinyVector<double,2> normal);
	TinyVector<double,2> getNormal();

	bool hasA();
	bool hasB();

	bool isBoundary();
	bool isInterior();
	void setIsBoundary(bool isBoundaryFace);
};

class Vertex : public GridElement{
	Coord c;
public:
	double & operator() (int i);
	Coord getCoord();
	void setCoord(Coord c);
};

class Grid{
public:
	double h;

	Array<Cell*,2> cells;
	vector<Vertex*> vertices;
	vector<Face*> faces;

	Grid * coarseGrid;
	Grid * fineGrid;
	bool hasCoarse;
	bool hasFine;

	int iMin, iMax, jMin, jMax;
	Range xRange, yRange, faceRange, vertexRange;

	void setHasCoarse(bool hasCoarse);
public:
	void tile(double xMin, double yMin, double xMax, double yMax);

public:
	~Grid();
	Grid();
	double getH();
	void setH(double h);
	CellDoubleArray makeCellDoubleArray();
	FaceDoubleArray makeFaceDoubleArray();
	VertexDoubleArray makeVertexDoubleArray();
	void resizeCellDoubleArray(CellDoubleArray & c);
	void resizeFaceDoubleArray(FaceDoubleArray & f);
	void resizeVertexDoubleArray(VertexDoubleArray & v);
	Face * createFace();
	Vertex * createVertex();
	Cell * createCell(int i, int j);
	void addFace(Face * f);
	void addVertex(Vertex * v);
	void addCell(Cell * c, int i, int j);
	double getVolume();
	Type getCellType(int i, int j);
	Type getFaceType(int i, int j, Direction d);
	bool isUncovered(int i, int j);
	bool isCovered(int i, int j);
	bool isRegular(int i, int j);
	bool isIrregular(int i, int j);
	bool isFaceUncovered(int i, int j, Direction d);
	bool isFaceCovered(int i, int j, Direction d);
	bool isFaceRegular(int i, int j, Direction d);
	bool isFaceIrregular(int i, int j, Direction d);
	Array<Type,2> getCellTypes();
	Array<Type,1> getFaceTypes();
	Array<Type,1> getVertexTypes();
	CellDoubleArray getCellX();
	CellDoubleArray getCellY();
	CellDoubleArray getCellCentroidX();
	CellDoubleArray getCellCentroidY();
	CellDoubleArray getVolumes();
	FaceDoubleArray getFaceX();
	FaceDoubleArray getFaceY();
	VertexDoubleArray getVertexX();
	VertexDoubleArray getVertexY();

	virtual Grid * getCoarse() = 0;
};

class Circle : public Grid{
	double r, offset;
	void grid();
	bool contains(Cell * c);
	bool contains(Face * f);
	bool contains(Vertex * v);
	bool contains(Coord c);
	bool contains(double x, double y);
public:
	Circle(double h, double r, double offset);
	void setR(double r);
	double getR();
	void setOffset(double offset);
	double getOffset();
	Grid * getCoarse();
};

class GridArray{
	Grid * g;
public:
	void setGrid(Grid * g);
	Grid * getGrid();
};

class Arrays{
public:
	static CellDoubleArray makeCellDoubleArray(Grid * g);
	static FaceDoubleArray makeFaceDoubleArray(Grid * g);
	static VertexDoubleArray makeVertexDoubleArray(Grid * g);
};

/*class CellDoubleArray : public Array<double,2>, public GridArray{

};

class FaceDoubleArray : public Array<double,1>, public GridArray{

};

class VertexDoubleArray : public Array<double,1>, public GridArray{

};*/

}
}

#endif /* GEOMETRY_H_ */
