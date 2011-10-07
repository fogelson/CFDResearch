/*
 * Geometry.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <blitz/array.h>
#include <list>
#include <blitz/tinyvec-et.h>

using namespace blitz;

namespace CFD{
	namespace Geometry{
		class Grid;
		class Rectangle;

		typedef TinyVector<double,2> Coord;
		typedef Array<Coord,2> CoordArray;

		typedef int Direction;
		const Direction N = 0;
		const Direction S = 1;
		const Direction E = 2;
		const Direction W = 3;
		const Direction B = 4;
		const Direction NE = 4;
		const Direction SE = 5;
		const Direction NW = 6;
		const Direction SW = 7;
		const Direction C = 8;

		const int CARDINAL_DIRECTIONS = 8;

		typedef int Type;
		const Type REGULAR = 0;
		const Type IRREGULAR = 1;
		const Type COVERED = 2;

		typedef Array<double,2> CellDoubleArray;
		typedef Array<Type,2> CellTypeArray;

		typedef TinyVector<double,5> FaceDouble;
		typedef Array<FaceDouble,2> FaceDoubleArray;
		typedef TinyVector<Coord,5> FaceCoord;
		typedef Array<FaceCoord,2> FaceCoordArray;
		typedef TinyVector<Type,5> FaceType;
		typedef Array<FaceType,2> FaceTypeArray;

		typedef TinyVector<double,9> Coefficients;
		typedef Array<Coefficients,2> CoefficientArray;

		struct Edge{
			Coord end1, end2;
			double length(){
				return sqrt(pow2(end1(0) - end2(0)) + pow2(end1(1) - end2(1)));
			}
		};


/*		Direction getCornerDirection(Direction a, Direction b){
			if(a == N){
				if(b == E){
					return NE;
				}
				if(b == W){
					return NW;
				}
			}
			if(a == E){
				if(b == N){
					return NE;
				}
				if(b == S){
					return SE;
				}
			}
			if(a == S){
				if(b == E){
					return SE;
				}
				if(b == W){
					return SW;
				}
			}
			if(a == W){
				if(b == N){
					return NE;
				}
				if(b == S){
					return SW;
				}
			}
			cout << "Tried to return corner of two opposing faces." << endl;
			return -1;
		}*/

/*		class Coefficient{
		public:
			int i, j;
			double c;
			Coefficient(int i, int j, double c){
				this->i = i;
				this->j = j;
				this->c = c;
			}
			Coefficient(int i, int j){
				this->i = i;
				this->j = j;
				c = 0;
			}
		};

		double operator + (Coefficient c, double d){
			return c.c + d;
		}
		double operator + (double d, Coefficient c){
			return d + c.c;
		}
		double operator - (Coefficient c, double d){
			return c.c - d;
		}
		double operator - (double d, Coefficient c){
			return d - c.c;
		}
		double operator * (Coefficient c, double d){
			return c.c * d;
		}
		double operator * (double d, Coefficient c){
			return d * c.c;
		}
		double operator / (Coefficient c, double d){
			return c.c / d;
		}
		double operator / (double d, Coefficient c){
			return d / c.c;
		}
		double operator + (Coefficient a, Coefficient b){
			return a.c + b.c;
		}
		double operator - (Coefficient a, Coefficient b){
			return a.c - b.c;
		}
		double operator * (Coefficient a, Coefficient b){
			return a.c * b.c;
		}
		double operator / (Coefficient a, Coefficient b){
			return a.c / b.c;
		}*/

		class Grid{
			Grid * coarseGrid;
			bool hasCoarsened;
		public:
			Grid(){
				hasCoarsened = false;
			}
			virtual ~Grid(){
				if(hasCoarsened){
					delete coarseGrid;
				}
			}
			double h;
			CoordArray centers;
			CoordArray cellCentroids;
			CellDoubleArray volumeFractions;
			CellTypeArray cellTypes;

			CellDoubleArray lengthRatios;

			Array<int,2> numberOfVertices;
			Array<TinyVector<Coord,5>,2> vertices;

			FaceCoordArray centroids;
			FaceDoubleArray areaFractions;
			FaceTypeArray faceTypes;
			Array<TinyVector<double,2>,2> outwardNormals;

			int iMin, jMin, iMax, jMax;
			Range xRange, yRange;

			void resizeAll(){
				centers.resize(xRange,yRange);
				volumeFractions.resize(xRange,yRange);
				cellTypes.resize(xRange,yRange);
				centroids.resize(xRange,yRange);
				areaFractions.resize(xRange,yRange);
				faceTypes.resize(xRange,yRange);
				outwardNormals.resize(xRange,yRange);
				cellCentroids.resize(xRange,yRange);
				numberOfVertices.resize(xRange,yRange);
				vertices.resize(xRange,yRange);
				lengthRatios.resize(xRange,yRange);
			}

			void calculateRatios(){
				for(int i = iMin; i <= iMax; i++){
					for(int j = jMin; j <= jMax; j++){
						if(isRegular(i,j)){
							lengthRatios(i,j) = 1;
						}
						if(isIrregular(i,j)){
							list<double> areas;
							if(isFaceUncovered(i,j,N)){
								areas.push_front(areaFractions(i,j)(N));
							}
							if(isFaceUncovered(i,j,E)){
								areas.push_front(areaFractions(i,j)(E));
							}
							if(isFaceUncovered(i,j,S)){
								areas.push_front(areaFractions(i,j)(S));
							}
							if(isFaceUncovered(i,j,W)){
								areas.push_front(areaFractions(i,j)(W));
							}
							double aMin = 1, aMax = 0;
							list<double>::iterator it;
							for(it = areas.begin(); it != areas.end(); it++){
								if(*it < aMin){
									aMin = *it;
								}
								if(*it > aMax){
									aMax = *it;
								}
							}
							lengthRatios(i,j) = aMax/aMin;
						}
						else{
							lengthRatios(i,j) = -1;
						}
					}
				}
			}

			void countVertices(int i, int j){
				int out;
				if(cellTypes(i,j) == REGULAR){
					out = 4;
				}
				if(cellTypes(i,j) == COVERED){
					out = 0;
				}
				else{
					out = 1;
					if(faceTypes(i,j)(N) != COVERED){
						out++;
					}
					if(faceTypes(i,j)(E) != COVERED){
						out++;
					}
					if(faceTypes(i,j)(S) != COVERED){
						out++;
					}
					if(faceTypes(i,j)(W) != COVERED){
						out++;
					}
				}
				numberOfVertices(i,j) = out;
			}
			void getVertices(int i, int j){
				Coord z;
				z(0) = 0;
				z(1) = 0;
				vertices(i,j) = z;
				int current = 0;
				int n = numberOfVertices(i,j);
				bool addedN = false, addedE = false, addedS = false, addedW = false;

				Coord NWfromN, NEfromN, NEfromE, SEfromE, SEfromS, SWfromS, SWfromW, NWfromW;

				if(faceTypes(i,j)(N) != COVERED){
					double distanceToVertex = (h/2)*areaFractions(i,j)(N);
					NWfromN(0) = centroids(i,j)(N)(0) - distanceToVertex;
					NWfromN(1) = centroids(i,j)(N)(1);

					NEfromN(0) = centroids(i,j)(N)(0) + distanceToVertex;
					NEfromN(1) = centroids(i,j)(N)(1);
					addedN = true;
				}
				if(faceTypes(i,j)(E) != COVERED){
					double distanceToVertex = (h/2)*areaFractions(i,j)(E);
					NEfromE(0) = centroids(i,j)(E)(0);
					NEfromE(1) = centroids(i,j)(E)(1) + distanceToVertex;

					SEfromE(0) = centroids(i,j)(E)(0);
					SEfromE(1) = centroids(i,j)(E)(1) - distanceToVertex;
				}
				if(faceTypes(i,j)(S) != COVERED){
					double distanceToVertex = (h/2)*areaFractions(i,j)(S);
					SEfromS(0) = centroids(i,j)(S)(0) + distanceToVertex;
					SEfromS(1) = centroids(i,j)(S)(1);

					SWfromS(0) = centroids(i,j)(S)(0) - distanceToVertex;
					SWfromS(1) = centroids(i,j)(S)(1);
				}
				if(faceTypes(i,j)(W) != COVERED){
					double distanceToVertex = (h/2)*areaFractions(i,j)(W);
					SWfromW(0) = centroids(i,j)(W)(0);
					SWfromW(1) = centroids(i,j)(W)(1) - distanceToVertex;

					NWfromW(0) = centroids(i,j)(W)(0);
					NWfromW(1) = centroids(i,j)(W)(1) + distanceToVertex;
				}

				if(isFaceUncovered(i,j,N)){
					vertices(i,j)(current) = NWfromN;
					current++;

					vertices(i,j)(current) = NEfromN;
					current++;
				}
				if(isFaceUncovered(i,j,E)){
					if(isFaceCovered(i,j,N) || centroids(i,j)(E)(1) < centers(i,j)(1)){
						vertices(i,j)(current) = NEfromE;
						current++;
					}

					vertices(i,j)(current) = SEfromE;
					current++;
				}
				if(isFaceUncovered(i,j,S)){
					if(isFaceCovered(i,j,E) || centroids(i,j)(S)(0) < centers(i,j)(0)){
						vertices(i,j)(current) = SEfromS;
						current++;
					}

					vertices(i,j)(current) = SWfromS;
					current++;
				}
				if(isFaceUncovered(i,j,W)){
					if(isFaceCovered(i,j,S) || centroids(i,j)(W)(1) > centers(i,j)(1)){
						vertices(i,j)(current) = SWfromW;
						current++;
					}
					if(isFaceCovered(i,j,N) || centroids(i,j)(W)(1) < centers(i,j)(1)){
						vertices(i,j)(current) = NWfromW;
						current++;
					}
				}

				if(current > n){
					cout << "Error. Added too many vertices in cell " << i << ", " << j << endl;
				}
			}

			CellDoubleArray makeCellDoubleArray(){
				CellDoubleArray c;
				c.resize(xRange, yRange);
				c = 0;
				return c;
			}

			Coord calculateCellCentroid(int i, int j){
				if(getCellType(i,j) != IRREGULAR){
					return centers(i,j);
				}
				/* Calculates the centroid of a polygon using
				 * the formula given on Wikipedia:
				 * http://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
				 */
				else{
					double A = 0;
					for(int k = 0; k < numberOfVertices(i,j); k++){
						int kP = (k + 1) % numberOfVertices(i,j);
						A += vertices(i,j)(k)(0)*vertices(i,j)(kP)(1) - vertices(i,j)(kP)(0)*vertices(i,j)(k)(1);
					}
					A = A/2;
					double Cx = 0, Cy = 0;
					for(int k = 0; k < numberOfVertices(i,j); k++){
						int kP = (k + 1) % numberOfVertices(i,j);
						double temp = vertices(i,j)(k)(0)*vertices(i,j)(kP)(1) - vertices(i,j)(kP)(0)*vertices(i,j)(k)(1);
						Cx += (vertices(i,j)(k)(0) + vertices(i,j)(kP)(0))*temp;
						Cy += (vertices(i,j)(k)(1) + vertices(i,j)(kP)(1))*temp;
					}
					Cx = Cx/(6*A);
					Cy = Cy/(6*A);
					Coord out;
					out(0) = Cx;
					out(1) = Cy;
					return out;
				}
			}

			Type getCellType(int i, int j){
				if(i < iMin || j < jMin || i > iMax || j > jMax)
					return COVERED;
				else
					return cellTypes(i,j);
			}
			Type getFaceType(int i, int j, Direction d){
				if(isCovered(i,j))
					return COVERED;
				else
					return faceTypes(i,j)(d);
			}
			bool isFaceUncovered(int i, int j, Direction d){
				return (getFaceType(i,j,d) != COVERED);
			}
			bool isFaceCovered(int i, int j, Direction d){
				return (getFaceType(i,j,d) == COVERED);
			}
			bool isFaceRegular(int i, int j, Direction d){
				return (getFaceType(i,j,d) == REGULAR);
			}
			bool isFaceIrregular(int i, int j, Direction d){
				return (getFaceType(i,j,d) == IRREGULAR);
			}

			bool isUncovered(int i, int j){
				return (getCellType(i,j) != COVERED);
			}

			bool isCovered(int i, int j){
				return (getCellType(i,j) == COVERED);
			}

			bool isIrregular(int i, int j){
				return (getCellType(i,j) == IRREGULAR);
			}

			bool isRegular(int i, int j){
				return (getCellType(i,j) == REGULAR);
			}

			double getCellEdge(int i, int j, Direction d){
				Coord c = centers(i,j);
				double e;
				if(d == N){
					e = c(1) + h/2;
				}
				else if(d == S){
					e = c(1) - h/2;
				}
				else if(d == E){
					e = c(0) + h/2;
				}
				else if(d == W){
					e = c(0) - h/2;
				}
				return e;
			}

			int calculateFacesCovered(int i, int j){
				int facesCovered = 0;
				FaceType ft = faceTypes(i,j);
				if(ft(N) == COVERED){
					facesCovered++;
				}
				if(ft(S) == COVERED){
					facesCovered++;
				}
				if(ft(E) == COVERED){
					facesCovered++;
				}
				if(ft(W) == COVERED){
					facesCovered++;
				}
				if(ft(B) == COVERED){
					facesCovered++;
				}
				return facesCovered;
			}
			virtual double calculateVolumeFraction(int i, int j){
				if(getCellType(i,j) == REGULAR){
					return 1;
				}
				else if(getCellType(i,j) == COVERED){
					return 0;
				}
				else{
					double A = 0;
					for(int k = 0; k < numberOfVertices(i,j); k++){
						int kP = (k + 1) % numberOfVertices(i,j);
						A += vertices(i,j)(k)(0)*vertices(i,j)(kP)(1) - vertices(i,j)(kP)(0)*vertices(i,j)(k)(1);
					}
					A = abs(A/2);
					return A/pow2(h);
				}
			}

			virtual double calculateVolumeFractionHmmm(int i, int j){
				if(cellTypes(i,j) == REGULAR){
					return 1;
				}
				if(cellTypes(i,j) == COVERED){
					return 0;
				}
				if(faceTypes(i,j)(B) == COVERED || faceTypes(i,j)(B) == REGULAR){
					// Something has gone wrong, since the boundary face should
					// always be labeled irregular in an irregular cell
					return -10;
				}
				// Returns the number of covered faces
				int facesCovered = calculateFacesCovered(i,j);
				switch(facesCovered){
				case 0: {// The interior is a square with a triangle cut out of a corner.
						 // So we find the two irregular faces, and the area of the triangle
						 // at their corner, and subtract that from 1.
					double a1, a2;
					bool foundFirstFace = false, foundSecondFace = false;
					if(!foundSecondFace && faceTypes(i,j)(N) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(N);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(N);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(S) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(S);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(S);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(E) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(E);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(E);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(W) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(W);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(W);
							foundSecondFace = true;
						}
					}
					return 1 - (1 - a1)*(1-a2)/2;
					break;
				}
				case 1: { /* This is the most complicated case, because the interior region
							 is a trapezoid. So we need to find all three uncovered faces,
							 determine which two are parallel to one another, and then
							 use the formula for the area of the trapezoid.
						  */
					double a1, a2, aH; // a1 and a2 are the parallel faces, aH is the trapezoid height face.
					if(faceTypes(i,j)(N) != COVERED){
						if(faceTypes(i,j)(S) != COVERED){
							a1 = areaFractions(i,j)(N);
							a2 = areaFractions(i,j)(S);
							if(faceTypes(i,j)(E) != COVERED){ //NSE
								aH = areaFractions(i,j)(E);
							}
							else{
								aH = areaFractions(i,j)(W); //NSW
							}
						}
						else{
							aH = areaFractions(i,j)(N); // EWN
							a1 = areaFractions(i,j)(E);
							a2 = areaFractions(i,j)(W);
						}
					}
					else{
						aH = areaFractions(i,j)(S); // EWS
						a1 = areaFractions(i,j)(E);
						a2 = areaFractions(i,j)(W);
					}
					return aH*(a1 + a2)/2;
					break;
				}
				case 2: {// The interior region is a triangle, so we find the two uncovered faces.
					/*cout << "Discovered a triangle. Here is a lot about it." << endl;
					cout << "\ta_N = " << areaFractions(i,j)(N);
					cout << ", \t\t Type_N = " << faceTypes(i,j)(N) << endl;
					cout << "\ta_S = " << areaFractions(i,j)(S);
					cout << ", \t\t Type_S = " << faceTypes(i,j)(S) << endl;
					cout << "\ta_E = " << areaFractions(i,j)(E);
					cout << ", \t\t Type_E = " << faceTypes(i,j)(E) << endl;
					cout << "\ta_W = " << areaFractions(i,j)(W);
					cout << ", \t\t Type_W = " << faceTypes(i,j)(W) << endl;*/
					//return 1;
					double a1, a2;
					bool foundFirstFace = false, foundSecondFace = false;
					if(!foundSecondFace && faceTypes(i,j)(N) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(N);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(N);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(S) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(S);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(S);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(E) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(E);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(E);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(W) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(W);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(W);
							foundSecondFace = true;
						}
					}
					return a1*a2/2.0;
					break;
				}
				case 3:
					// Something is wrong if this ever gets called.
					cout << "oops" << endl;
					return -1;
					break;
				case 4:
					// Likewise here
					cout << "oops" << endl;
					return -2;
					break;
				default:
					// Likewise here
					cout << "oops" << endl;
					return -3;
					break;
				}
			}
			virtual void respace(double hNew) = 0;
			virtual Grid * copy() const = 0;
			Grid * coarsen(){
				if(hasCoarsened){
					return coarseGrid;
				}
				else{
					coarseGrid = this->copy();
					coarseGrid->respace(2*h);
					hasCoarsened = true;
					return coarseGrid;
				}
			}
			template<typename T>
			void resizeArray(Array<T,2> &u){
				u.resize(xRange,yRange);
			}
		};


		class NewCircle : public Grid{
			double r, offset;

			bool isInterior(double x, double y){
				return (pow2(x) + pow2(y) < pow2(r));
			}
			bool isExterior(double x, double y){
				return (pow2(x) + pow2(y) > pow2(r));
			}
			bool isOnBoundary(double x, double y){
				return (pow2(x) + pow2(y) == pow2(r));
			}
			bool isInterior(Coord c){
				return isInterior(c(0),c(1));
			}
			bool isExterior(Coord c){
				return isExterior(c(0),c(1));
			}
			bool isOnBoundary(Coord c){
				return isOnBoundary(c(0),c(1));
			}

		public:
			void make(){
				double width = 2 * r + offset;
				iMin = 1;
				jMin = 1;
				iMax = ceil(width/h);
				jMax = ceil(width/h);
				iMax = iMax + (iMax % 2);
				jMax = jMax + (jMax % 2);
				xRange.setRange(1,iMax,1);
				yRange.setRange(1,jMax,1);

				centers.resize(xRange,yRange);
				volumeFractions.resize(xRange,yRange);
				cellTypes.resize(xRange,yRange);
				centroids.resize(xRange,yRange);
				areaFractions.resize(xRange,yRange);
				faceTypes.resize(xRange,yRange);
				outwardNormals.resize(xRange,yRange);
				cellCentroids.resize(xRange,yRange);
				numberOfVertices.resize(xRange,yRange);
				vertices.resize(xRange,yRange);

				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						centers(i,j)(0) = h*(i-1) - r - offset + (h/2);
						centers(i,j)(1) = h*(j-1) - r - offset + (h/2);

						Coord NECorner, NWCorner, SECorner, SWCorner;

						NECorner(0) = h*i - r - offset;
						NECorner(1) = h*j - r - offset;

						NWCorner(0) = h*(i-1) - r - offset;
						NWCorner(1) = h*j - r - offset;

						SECorner(0) = h*i - r - offset;
						SECorner(1) = h*(j-1) - r - offset;

						SWCorner(0) = h*(i-1) - r - offset;
						SWCorner(1) = h*(j-1) - r - offset;

						Edge NEdge, SEdge, EEdge, WEdge;

						NEdge.end1 = NECorner;
						NEdge.end2 = NWCorner;

						SEdge.end1 = SECorner;
						SEdge.end2 = SWCorner;

						EEdge.end1 = NECorner;
						EEdge.end2 = SECorner;

						WEdge.end1 = NWCorner;
						WEdge.end2 = SWCorner;


					}
				}
			}
		};

		class Circle : public Grid{
		public:
			double r;
			double x0, y0;
			double offset;
			void makeNope(){
				double width = 2 * r + offset;
				centers.resize(xRange,yRange);
				volumeFractions.resize(xRange,yRange);
				cellTypes.resize(xRange,yRange);
				centroids.resize(xRange,yRange);
				areaFractions.resize(xRange,yRange);
				faceTypes.resize(xRange,yRange);
				outwardNormals.resize(xRange,yRange);
				cellCentroids.resize(xRange,yRange);
				numberOfVertices.resize(xRange,yRange);
				vertices.resize(xRange,yRange);

				for(int i = iMin; i <= iMax; i++){
					for(int j = jMin; j <= jMax; j++){
						// NE
						vertices(i,j)(0)(0) = i*h - offset - r + x0;
						vertices(i,j)(0)(1) = j*h - offset - r + y0;

						// SE
						vertices(i,j)(1)(0) = i*h - offset - r + x0;
						vertices(i,j)(1)(1) = (j-1)*h - offset - r + y0;

						// SW
						vertices(i,j)(2)(0) = (i-1)*h - offset - r + x0;
						vertices(i,j)(2)(1) = (j-1)*h - offset - r + y0;

						// NW
						vertices(i,j)(3)(0) = (i-1)*h - offset - r + x0;
						vertices(i,j)(3)(1) = j*h - offset - r + y0;

						vertices(i,j)(4)(0) = 0;
						vertices(i,j)(4)(1) = 0;

						numberOfVertices(i,j) = 4;

						centers(i,j)(0) = (i-1)*h + (h/2) - offset - r + x0;
						centers(i,j)(1) = (j-1)*h + (h/2) - offset - r + y0;
					}
				}

			}
			void make(){
				double width = 2 * r + offset;
				iMin = 1;
				jMin = 1;
				iMax = ceil(width/h);
				jMax = ceil(width/h);
				iMax = iMax + (iMax % 2);
				jMax = jMax + (jMax % 2);
				xRange.setRange(1,iMax,1);
				yRange.setRange(1,jMax,1);

				resizeAll();

				/*centers.resize(xRange,yRange);
				volumeFractions.resize(xRange,yRange);
				cellTypes.resize(xRange,yRange);
				centroids.resize(xRange,yRange);
				areaFractions.resize(xRange,yRange);
				faceTypes.resize(xRange,yRange);
				outwardNormals.resize(xRange,yRange);
				cellCentroids.resize(xRange,yRange);
				numberOfVertices.resize(xRange,yRange);
				vertices.resize(xRange,yRange);*/

				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						cellTypes(i,j) = REGULAR;
						faceTypes(i,j) = REGULAR;
						faceTypes(i,j)(B) = COVERED;
						volumeFractions(i,j) = 1;
						areaFractions(i,j) = 1;
						areaFractions(i,j)(B) = 0;
						outwardNormals(i,j) = 0;

						centers(i,j)(0) = (i-1)*h + (h/2) - offset - r + x0;
						centers(i,j)(1) = (j-1)*h + (h/2) - offset - r + y0;

						centroids(i,j)(N)(0) = centers(i,j)(0);
						centroids(i,j)(N)(1) = getCellEdge(i,j,N);

						centroids(i,j)(S)(0) = centers(i,j)(0);
						centroids(i,j)(S)(1) = getCellEdge(i,j,S);

						centroids(i,j)(E)(0) = getCellEdge(i,j,E);
						centroids(i,j)(E)(1) = centers(i,j)(1);

						centroids(i,j)(W)(0) = getCellEdge(i,j,W);
						centroids(i,j)(W)(1) = centers(i,j)(1);

						Coord NWCorner, NECorner, SWCorner, SECorner;
						NWCorner(0) = centroids(i,j)(W)(0);
						NWCorner(1) = centroids(i,j)(N)(1);

						NECorner(0) = centroids(i,j)(E)(0);
						NECorner(1) = centroids(i,j)(N)(1);

						SWCorner(0) = centroids(i,j)(W)(0);
						SWCorner(1) = centroids(i,j)(S)(1);

						SECorner(0) = centroids(i,j)(E)(0);
						SECorner(1) = centroids(i,j)(S)(1);

						bool isCovered = true;
						list<Coord> boundaryEdges;

						// North face
						if(contains(NWCorner) && contains(NECorner)){
							isCovered = false;
						}
						else if(contains(NWCorner)){
							isCovered = false;
							faceTypes(i,j)(N) = IRREGULAR;
							Coord b;
							b(1) = centroids(i,j)(N)(1);
							b(0) = sqrt(pow2(r) - pow2(b(1) - y0)) + x0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(N)(0) = (b(0) + NWCorner(0))/2;
							areaFractions(i,j)(N) = (b(0) - NWCorner(0))/h;
						}
						else if(contains(NECorner)){
							isCovered = false;
							faceTypes(i,j)(N) = IRREGULAR;
							Coord b;
							b(1) = centroids(i,j)(N)(1);
							b(0) = -sqrt(pow2(r) - pow2(b(1) - y0)) + x0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(N)(0) = (b(0) + NECorner(0))/2;
							areaFractions(i,j)(N) = (NECorner(0) - b(0))/h;
						}
						else{
							faceTypes(i,j)(N) = COVERED;
							centroids(i,j)(N) = 0;
							areaFractions(i,j)(N) = 0;
						}

						// South face
						if(contains(SWCorner) && contains(SECorner)){
							isCovered = false;
						}
						else if(contains(SWCorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(S) = IRREGULAR;
							Coord b;
							b(1) = centroids(i,j)(S)(1);
							b(0) = sqrt(pow2(r) - pow2(b(1) - y0)) + x0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(S)(0) = (b(0) + SWCorner(0))/2;
							areaFractions(i,j)(S) = (b(0) - SWCorner(0))/h;
						}
						else if(contains(SECorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(S) = IRREGULAR;
							Coord b;
							b(1) = centroids(i,j)(S)(1);
							b(0) = -sqrt(pow2(r) - pow2(b(1) - y0)) + x0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(S)(0) = (b(0) + SECorner(0))/2;
							areaFractions(i,j)(S) = (SECorner(0) - b(0))/h;
						}
						else{
							faceTypes(i,j)(S) = COVERED;
							centroids(i,j)(S) = 0;
							areaFractions(i,j)(S) = 0;
						}

						// East face
						if(contains(NECorner) && contains(SECorner)){
							isCovered = false;
						}
						else if(contains(SECorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(E) = IRREGULAR;
							Coord b;
							b(0) = centroids(i,j)(E)(0);
							b(1) = sqrt(pow2(r) - pow2(b(0) - x0)) + y0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(E)(1) = (b(1) + SECorner(1))/2;
							areaFractions(i,j)(E) = (b(1) - SECorner(1))/h;
						}
						else if(contains(NECorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(E) = IRREGULAR;
							Coord b;
							b(0) = centroids(i,j)(E)(0);
							b(1) = -sqrt(pow2(r) - pow2(b(0) - x0)) + y0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(E)(1) = (b(1) + NECorner(1))/2;
							areaFractions(i,j)(E) = (NECorner(1) - b(1))/h;
						}
						else{
							faceTypes(i,j)(E) = COVERED;
							centroids(i,j)(E) = 0;
							areaFractions(i,j)(E) = 0;
						}

						// West face
						if(contains(NWCorner) && contains(SWCorner)){
							isCovered = false;
						}
						else if(contains(SWCorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(W) = IRREGULAR;
							Coord b;
							b(0) = centroids(i,j)(W)(0);
							b(1) = sqrt(pow2(r) - pow2(b(0) - x0)) + y0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(W)(1) = (b(1) + SWCorner(1))/2;
							areaFractions(i,j)(W) = (b(1) - SWCorner(1))/h;
						}
						else if(contains(NWCorner)){
							isCovered = false;
							cellTypes(i,j) = IRREGULAR;
							faceTypes(i,j)(W) = IRREGULAR;
							Coord b;
							b(0) = centroids(i,j)(W)(0);
							b(1) = -sqrt(pow2(r) - pow2(b(0) - x0)) + y0;
							if(boundaryEdges.size() < 2){
								boundaryEdges.push_front(b);
							}

							centroids(i,j)(W)(1) = (b(1) + NWCorner(1))/2;
							areaFractions(i,j)(W) = (NWCorner(1) - b(1))/h;
						}
						else{
							faceTypes(i,j)(W) = COVERED;
							centroids(i,j)(W) = 0;
							areaFractions(i,j)(W) = 0;
						}

						// Boundary face
						if(!boundaryEdges.empty()){
							faceTypes(i,j)(B) = IRREGULAR;
							Coord b1, b2;
							b1 = boundaryEdges.front();
							b2 = boundaryEdges.back();

							centroids(i,j)(B) = (b1 + b2)/2;
							areaFractions(i,j)(B) = sqrt(pow2(b1(0) - b2(0)) + pow2(b1(1) - b2(1)))/h;
							outwardNormals(i,j) = centroids(i,j)(B) - Coord(x0,y0);
							outwardNormals(i,j) = outwardNormals(i,j)/sqrt(pow2(outwardNormals(i,j)(0)) + pow2(outwardNormals(i,j)(1)));
						}

						if(isCovered){
							cellTypes(i,j) = COVERED;
							faceTypes(i,j) = COVERED;
							volumeFractions(i,j) = 0;
							areaFractions(i,j) = 0;
							outwardNormals(i,j) = 0;
						}
					}
				}
				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						countVertices(i,j);
						getVertices(i,j);
						volumeFractions(i,j) = calculateVolumeFraction(i,j);
						cellCentroids(i,j) = calculateCellCentroid(i,j);
					}
				}
			}
			/*virtual double calculateVolumeFraction(int i, int j){
				if(cellTypes(i,j) == COVERED){
					return 0;
				}
				else if(cellTypes(i,j) == REGULAR){
					return 1;
				}
				else{
					double A;
					if(numberOfVertices(i,j) == 3){
						double x0, x1, y0;
						if(faceTypes(i,j)(N) != COVERED && faceTypes(i,j)(E) != COVERED){
							x0 = abs(getCellEdge(i,j,E));
							y0 = abs(getCellEdge(i,j,N));
							x1 = x0 + h*areaFractions(i,j)(N);
						}
						else if(faceTypes(i,j)(E) != COVERED && faceTypes(i,j)(S) != COVERED){
							x0 = abs(getCellEdge(i,j,E));
							y0 = abs(getCellEdge(i,j,S));
							x1 = x0 + h*areaFractions(i,j)(S);
						}
						else if(faceTypes(i,j)(S) != COVERED && faceTypes(i,j)(W) != COVERED){
							x0 = abs(getCellEdge(i,j,W));
							y0 = abs(getCellEdge(i,j,S));
							x1 = x0 + h*areaFractions(i,j)(S);
						}
						else if(faceTypes(i,j)(W) != COVERED && faceTypes(i,j)(N) != COVERED){
							x0 = abs(getCellEdge(i,j,W));
							y0 = abs(getCellEdge(i,j,N));
							x1 = x0 + h*areaFractions(i,j)(N);
						}
						A = -(pow2(r)/2.0)*(acos(x1/r) - acos(x0/r)) + (x1*sqrt(pow2(r) - pow2(x1)) - x0*sqrt(pow2(r) - pow2(x0))) - y0*(x1 - x0);

					}
					else if(numberOfVertices(i,j) == 4){
						double x0, x1, y0, y1;

					}
					else if(numberOfVertices(i,j) == 5){

					}

					return A;
				}
			}*/
		public:
			virtual double calculateVolumeFractionExactly(int i, int j){
				if(cellTypes(i,j) == REGULAR){
					return 1;
				}
				if(cellTypes(i,j) == COVERED){
					return 0;
				}
				if(faceTypes(i,j)(B) == COVERED || faceTypes(i,j)(B) == REGULAR){
					// Something has gone wrong, since the boundary face should
					// always be labeled irregular in an irregular cell
					return -10;
				}
				// Returns the number of covered faces
				int facesCovered = calculateFacesCovered(i,j);
				switch(facesCovered){
				case 0: {// The interior is a square with a triangle cut out of a corner.
						 // So we find the two irregular faces, and the area of the triangle
						 // at their corner, and subtract that from 1.

					double x0, x1, y0, l;
					bool doIt = false;
					if(isFaceRegular(i,j,W) && isFaceRegular(i,j,S)){
						l = areaFractions(i,j)(N)*h;
						x1 = centroids(i,j)(E)(0);
						x0 = l + x1 - h;
						y0 = abs(centroids(i,j)(S)(1));
						doIt = true;
					}
					else if(isFaceRegular(i,j,W) && isFaceRegular(i,j,N)){
						l = areaFractions(i,j)(S)*h;
						x1 = centroids(i,j)(E)(0);
						x0 = l + x1 - h;
						y0 = abs(centroids(i,j)(N)(1));
						doIt = true;
					}
					else if(isFaceRegular(i,j,E) && isFaceRegular(i,j,S)){
						l = areaFractions(i,j)(W)*h;
						x1 = centroids(i,j)(N)(1);
						x0 = l + x1 - h;
						y0 = abs(centroids(i,j)(E)(0));
						doIt = true;
					}
					else if(isFaceRegular(i,j,E) && isFaceRegular(i,j,N)){
						l = areaFractions(i,j)(W)*h;
						x0 = centroids(i,j)(S)(1);
						x1 = h + x0 - l;
						y0 = abs(centroids(i,j)(E)(0));
						doIt = true;
					}

					doIt = false;

					if(doIt){
						double out = 0;
						out = (0.5*r)*(x1*sqrt(1-pow2(x1/r)) + r*asin(x1/r)
							- x0*sqrt(1-pow2(x0/r))
							- r*asin(x0/r)) - y0*(x1 - x0)
							+ l*h;
						return abs(out/pow2(h));
					}

					double a1, a2;
					bool foundFirstFace = false, foundSecondFace = false;
					if(!foundSecondFace && faceTypes(i,j)(N) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(N);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(N);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(S) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(S);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(S);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(E) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(E);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(E);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(W) == IRREGULAR){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(W);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(W);
							foundSecondFace = true;
						}
					}
					return 1 - (1 - a1)*(1-a2)/2;
					break;
				}
				case 1: {
					double x0, x1, y0;
					if(isFaceRegular(i,j,S)){
						x0 = centroids(i,j)(W)(0);
						x1 = centroids(i,j)(E)(0);
						y0 = abs(centroids(i,j)(S)(1));
					}
					else if(isFaceRegular(i,j,N)){
						x0 = centroids(i,j)(W)(0);
						x1 = centroids(i,j)(E)(0);
						y0 = abs(centroids(i,j)(N)(1));
					}
					else if(isFaceRegular(i,j,E)){
						x0 = centroids(i,j)(S)(1);
						x1 = centroids(i,j)(N)(1);
						y0 = abs(centroids(i,j)(E)(0));
					}
					else if(isFaceRegular(i,j,W)){
						x0 = centroids(i,j)(S)(1);
						x1 = centroids(i,j)(N)(1);
						y0 = abs(centroids(i,j)(W)(0));
					}

					double out = 0;
					out = (0.5*r)*(x1*sqrt(1-pow2(x1/r)) + r*asin(x1/r)
							- x0*sqrt(1-pow2(x0/r)) - r*asin(x0/r))
							- y0*(x1 - x0);
					return abs(out/pow2(h));
					break;

					/* This is the most complicated case, because the interior region
							 is a trapezoid. So we need to find all three uncovered faces,
							 determine which two are parallel to one another, and then
							 use the formula for the area of the trapezoid.
						  */
					double a1, a2, aH; // a1 and a2 are the parallel faces, aH is the trapezoid height face.
					if(faceTypes(i,j)(N) != COVERED){
						if(faceTypes(i,j)(S) != COVERED){
							a1 = areaFractions(i,j)(N);
							a2 = areaFractions(i,j)(S);
							if(faceTypes(i,j)(E) != COVERED){ //NSE
								aH = areaFractions(i,j)(E);
							}
							else{
								aH = areaFractions(i,j)(W); //NSW
							}
						}
						else{
							aH = areaFractions(i,j)(N); // EWN
							a1 = areaFractions(i,j)(E);
							a2 = areaFractions(i,j)(W);
						}
					}
					else{
						aH = areaFractions(i,j)(S); // EWS
						a1 = areaFractions(i,j)(E);
						a2 = areaFractions(i,j)(W);
					}
					return aH*(a1 + a2)/2;
					break;
				}
				case 2: {// The interior region is a triangle, so we find the two uncovered faces.
					/*cout << "Discovered a triangle. Here is a lot about it." << endl;
					cout << "\ta_N = " << areaFractions(i,j)(N);
					cout << ", \t\t Type_N = " << faceTypes(i,j)(N) << endl;
					cout << "\ta_S = " << areaFractions(i,j)(S);
					cout << ", \t\t Type_S = " << faceTypes(i,j)(S) << endl;
					cout << "\ta_E = " << areaFractions(i,j)(E);
					cout << ", \t\t Type_E = " << faceTypes(i,j)(E) << endl;
					cout << "\ta_W = " << areaFractions(i,j)(W);
					cout << ", \t\t Type_W = " << faceTypes(i,j)(W) << endl;*/
					//return 1;

					double a, b, t, s;

					if(isFaceUncovered(i,j,S) && isFaceUncovered(i,j,W)){
						a = areaFractions(i,j)(S)*h;
						b = areaFractions(i,j)(W)*h;
						s = abs(centroids(i,j)(S)(1));
						t = abs(centroids(i,j)(W)(0));
					}
					else if(isFaceUncovered(i,j,S) && isFaceUncovered(i,j,E)){
						a = areaFractions(i,j)(S)*h;
						b = areaFractions(i,j)(E)*h;
						s = abs(centroids(i,j)(S)(1));
						t = abs(centroids(i,j)(E)(0));
					}
					else if(isFaceUncovered(i,j,N) && isFaceUncovered(i,j,W)){
						a = areaFractions(i,j)(N)*h;
						b = areaFractions(i,j)(W)*h;
						s = abs(centroids(i,j)(N)(1));
						t = abs(centroids(i,j)(W)(0));
					}
					else if(isFaceUncovered(i,j,N) && isFaceUncovered(i,j,E)){
						a = areaFractions(i,j)(N)*h;
						b = areaFractions(i,j)(E)*h;
						s = abs(centroids(i,j)(N)(1));
						t = abs(centroids(i,j)(E)(0));
					}

					double out = 0;
					out = (0.5*r)*((t+a)*sqrt(1.0 - pow2((t+a)/r)) + r*asin((t+a)/r) - t*sqrt(1.0 - pow2(t/r)) - r*asin(t/r)) - a*s;

					return out/pow2(h);

					break;

					double a1, a2;
					bool foundFirstFace = false, foundSecondFace = false;
					if(!foundSecondFace && faceTypes(i,j)(N) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(N);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(N);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(S) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(S);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(S);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(E) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(E);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(E);
							foundSecondFace = true;
						}
					}
					if(!foundSecondFace && faceTypes(i,j)(W) != COVERED){
						if(!foundFirstFace){
							a1 = areaFractions(i,j)(W);
							foundFirstFace = true;
						}
						else{
							a2 = areaFractions(i,j)(W);
							foundSecondFace = true;
						}
					}
					return a1*a2/2.0;
					break;
				}
				case 3:
					// Something is wrong if this ever gets called.
					cout << "oops" << endl;
					return -1;
					break;
				case 4:
					// Likewise here
					cout << "oops" << endl;
					return -2;
					break;
				default:
					// Likewise here
					cout << "oops" << endl;
					return -3;
					break;
				}
			}
			bool contains(Coord c){
				return pow2(c(0) - x0) + pow2(c(1) - y0) < pow2(r);
			}
			Circle(double h, double r, double offset){
				this->h = h;
				this->r = r;
				x0 = 0;
				y0 = 0;
				this->offset = offset;
				make();
			}
			Circle(double h, double r, double x0, double y0, double offset){
				this->h = h;
				this->r = r;
				this->x0 = x0;
				this->y0 = y0;
				this->offset = offset;
				make();
			}
			Circle(const Circle & u){
				this->h = u.h;
				this->r = u.r;
				this->x0 = u.x0;
				this->y0 = u.y0;
				this->offset = u.offset;
				make();
				this->volumeFractions = u.volumeFractions;
			}
			Grid * copy() const{
				return new Circle(*this);
			}
			void respace(double hNew){
				this->h = hNew;
				make();
			}
		};

		class UnitSquareBetterBoundaries : public Grid{
			void make(){
				double width = 1 + offset;
				iMin = 1;
				jMin = 1;
				iMax = ceil(width/h);
				jMax = ceil(width/h);
				iMax = iMax + (iMax % 2);
				jMax = jMax + (jMax % 2);
				xRange.setRange(1,iMax,1);
				yRange.setRange(1,jMax,1);

				resizeAll();

				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						cellTypes(i,j) = REGULAR;
						faceTypes(i,j) = REGULAR;
						faceTypes(i,j)(B) = COVERED;
						volumeFractions(i,j) = 1;
						areaFractions(i,j) = 1;
						areaFractions(i,j)(B) = 0;
						outwardNormals(i,j) = 0;

						centers(i,j)(0) = (i-1)*h + (h/2) - offset;
						centers(i,j)(1) = (j-1)*h + (h/2) - offset;

						centroids(i,j)(N)(0) = centers(i,j)(0);
						centroids(i,j)(N)(1) = getCellEdge(i,j,N);

						centroids(i,j)(S)(0) = centers(i,j)(0);
						centroids(i,j)(S)(1) = getCellEdge(i,j,S);

						centroids(i,j)(E)(0) = getCellEdge(i,j,E);
						centroids(i,j)(E)(1) = centers(i,j)(1);

						centroids(i,j)(W)(0) = getCellEdge(i,j,W);
						centroids(i,j)(W)(1) = centers(i,j)(1);

						if(!contains(centers(i,j))
								&&
								((getCellEdge(i,j,N) < 0 || getCellEdge(i,j,N) > 1)
								&& (getCellEdge(i,j,S) < 0 || getCellEdge(i,j,S) > 1))
								||
								((getCellEdge(i,j,E) < 0 || getCellEdge(i,j,E) > 1)
								&& (getCellEdge(i,j,W) < 0 || getCellEdge(i,j,W) > 1))){
							cellTypes(i,j) = COVERED;
							faceTypes(i,j) = COVERED;
							volumeFractions(i,j) = 0;
							areaFractions(i,j) = 0;
							outwardNormals(i,j) = 0;
						}
						else{
							if(abs(centers(i,j)(0)) < h/2){
								if(abs(centers(i,j)(1)) < h/2){ // i = 1, j = 1
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(S) = COVERED;
									areaFractions(i,j)(S) = 0;
									centroids(i,j)(S) = 0;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(E) = IRREGULAR;
									centroids(i,j)(E)(0) = getCellEdge(i,j,E);
									centroids(i,j)(E)(1) = getCellEdge(i,j,N)/2;
									areaFractions(i,j)(E) = getCellEdge(i,j,N)/h;

									faceTypes(i,j)(B) = IRREGULAR;

									// Shouldn't use this centroid. Same with other corners.
									centroids(i,j)(B)(0) = centroids(i,j)(N)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(E)(1);
									areaFractions(i,j)(B) = areaFractions(i,j)(N) + areaFractions(i,j)(E);

									// Shouldn't use this outward normal. Same with other corners.
									outwardNormals(i,j)(0) = -areaFractions(i,j)(E)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = -areaFractions(i,j)(N)/areaFractions(i,j)(B);
								}
								else if(abs(centers(i,j)(1) - 1) < h/2){ // i = 1, j = jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(N) = COVERED;
									areaFractions(i,j)(N) = 0;
									centroids(i,j)(N) = 0;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(E) = IRREGULAR;
									centroids(i,j)(E)(0) = getCellEdge(i,j,E);
									centroids(i,j)(E)(1) = (1 + getCellEdge(i,j,S))/2;
									areaFractions(i,j)(E) = (1 - getCellEdge(i,j,S))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(S)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(E)(1);
									areaFractions(i,j)(B) = areaFractions(i,j)(S) + areaFractions(i,j)(E);

									outwardNormals(i,j)(0) = -areaFractions(i,j)(E)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = areaFractions(i,j)(S)/areaFractions(i,j)(B);
								}
								else{ // i = 1, j != 1, j != jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = 0;
									centroids(i,j)(B)(1) = centers(i,j)(1);
									areaFractions(i,j)(B) = 1;
									outwardNormals(i,j)(0) = -1;
									outwardNormals(i,j)(1) = 0;
								}
							}
							else if(abs(centers(i,j)(0) - 1) < h/2){
								if(abs(centers(i,j)(1)) < h/2){ // i = iMax, j = 1
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(S) = COVERED;
									areaFractions(i,j)(S) = 0;
									centroids(i,j)(S) = 0;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(W) = IRREGULAR;
									centroids(i,j)(W)(0) = getCellEdge(i,j,W);
									centroids(i,j)(W)(1) = getCellEdge(i,j,N)/2;
									areaFractions(i,j)(W) = getCellEdge(i,j,N)/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(N)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(W)(1);
									areaFractions(i,j)(B) = areaFractions(i,j)(N) + areaFractions(i,j)(W);

									outwardNormals(i,j)(0) = areaFractions(i,j)(W)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = -areaFractions(i,j)(N)/areaFractions(i,j)(B);
								}
								else if(abs(centers(i,j)(1) - 1) < h/2){ // i = iMax, j = jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(N) = COVERED;
									areaFractions(i,j)(N) = 0;
									centroids(i,j)(N) = 0;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(W) = IRREGULAR;
									centroids(i,j)(W)(0) = getCellEdge(i,j,W);
									centroids(i,j)(W)(1) = (1 + getCellEdge(i,j,S))/2;
									areaFractions(i,j)(W) = (1 - getCellEdge(i,j,S))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(S)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(W)(1);
									areaFractions(i,j)(B) = areaFractions(i,j)(S) + areaFractions(i,j)(W);

									outwardNormals(i,j)(0) = areaFractions(i,j)(W)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = areaFractions(i,j)(S)/areaFractions(i,j)(B);
								}
								else{ // i = iMax, j != 1, j != jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = 1;
									centroids(i,j)(B)(1) = centers(i,j)(1);
									areaFractions(i,j)(B) = 1;
									outwardNormals(i,j)(0) = 1;
									outwardNormals(i,j)(1) = 0;
								}
							}
							else if(abs(centers(i,j)(1)) < h/2){ // j = 1, i != 1, i != iMax
								cellTypes(i,j) = IRREGULAR;

								faceTypes(i,j)(S) = COVERED;
								areaFractions(i,j)(S) = 0;
								centroids(i,j)(S) = 0;

								faceTypes(i,j)(E) = IRREGULAR;
								centroids(i,j)(E)(0) = getCellEdge(i,j,E);
								centroids(i,j)(E)(1) = getCellEdge(i,j,N)/2;
								areaFractions(i,j)(E) = getCellEdge(i,j,N)/h;

								faceTypes(i,j)(W) = IRREGULAR;
								centroids(i,j)(W)(0) = getCellEdge(i,j,W);
								centroids(i,j)(W)(1) = getCellEdge(i,j,N)/2;
								areaFractions(i,j)(W) = getCellEdge(i,j,N)/h;

								faceTypes(i,j)(B) = IRREGULAR;
								centroids(i,j)(B)(0) = centers(i,j)(0);
								centroids(i,j)(B)(1) = 0;
								areaFractions(i,j)(B) = 1;
								outwardNormals(i,j)(0) = 0;
								outwardNormals(i,j)(1) = -1;
							}
							else if(abs(centers(i,j)(1) - 1) < h/2){ // j = jMax, i != 1, i != iMax
								cellTypes(i,j) = IRREGULAR;

								faceTypes(i,j)(N) = COVERED;
								areaFractions(i,j)(N) = 0;
								centroids(i,j)(N) = 0;

								faceTypes(i,j)(E) = IRREGULAR;
								centroids(i,j)(E)(0) = getCellEdge(i,j,E);
								centroids(i,j)(E)(1) = (1 + getCellEdge(i,j,S))/2;
								areaFractions(i,j)(E) = (1 - getCellEdge(i,j,S))/h;

								faceTypes(i,j)(W) = IRREGULAR;
								centroids(i,j)(W)(0) = getCellEdge(i,j,W);
								centroids(i,j)(W)(1) = (1 + getCellEdge(i,j,S))/2;
								areaFractions(i,j)(W) = (1 - getCellEdge(i,j,S))/h;

								faceTypes(i,j)(B) = IRREGULAR;
								centroids(i,j)(B)(0) = centers(i,j)(0);
								centroids(i,j)(B)(1) = 1;
								areaFractions(i,j)(B) = 1;
								outwardNormals(i,j)(0) = 0;
								outwardNormals(i,j)(1) = 1;
							}
						}
						for(int i = 1; i <= iMax; i++){
							for(int j = 1; j <= jMax; j++){
								volumeFractions(i,j) = calculateVolumeFraction(i,j);
								if(abs(centers(i,j)(0)) < h/2 && abs(centers(i,j)(1)) < h/2){
									volumeFractions(i,j) = areaFractions(i,j)(N)*areaFractions(i,j)(E);
								}
								if(abs(centers(i,j)(0)) < h/2 && abs(1 - centers(i,j)(1)) < h/2){
									volumeFractions(i,j) = areaFractions(i,j)(S) + areaFractions(i,j)(E);
								}
								if(abs(1 - centers(i,j)(0)) < h/2 && abs(centers(i,j)(1)) < h/2){
									volumeFractions(i,j) = areaFractions(i,j)(W) + areaFractions(i,j)(N);
								}
								if(abs(1 - centers(i,j)(0)) < h/2 && abs(1 - centers(i,j)(1)) < h/2){
									volumeFractions(i,j) = areaFractions(i,j)(W) + areaFractions(i,j)(S);
								}
								countVertices(i,j);
								getVertices(i,j);
								cellCentroids(i,j) = calculateCellCentroid(i,j);
							}
						}
					}
				}
			}
		public:
			double offset;
			bool contains(Coord c){
				double x = c(0), y = c(1);
				if(x >= 0 && x <= 1 && y >= 0 && y <= 1)
					return true;
				return false;
			}
			UnitSquareBetterBoundaries(double h, double offset){
				this->h = h;
				this->offset = offset;
				make();
			}
			UnitSquareBetterBoundaries(const UnitSquareBetterBoundaries & u){
				this->h = u.h;
				this->offset = u.offset;
				make();
				this->volumeFractions = u.volumeFractions;
			}
			Grid * copy() const{
				return new UnitSquareBetterBoundaries(*this);
			}
			void respace(double hNew){
				this->h = hNew;
				make();
			}
		};

		class UnitSquare : public Grid{
			void make(){
				double width = 1 + offset;
				iMin = 1;
				jMin = 1;
				iMax = ceil(width/h);
				jMax = ceil(width/h);
				iMax = iMax + (iMax % 2);
				jMax = jMax + (jMax % 2);
				xRange.setRange(1,iMax,1);
				yRange.setRange(1,jMax,1);

				centers.resize(xRange,yRange);
				volumeFractions.resize(xRange,yRange);
				cellTypes.resize(xRange,yRange);
				centroids.resize(xRange,yRange);
				areaFractions.resize(xRange,yRange);
				faceTypes.resize(xRange,yRange);
				outwardNormals.resize(xRange,yRange);
				cellCentroids.resize(xRange,yRange);
				numberOfVertices.resize(xRange,yRange);
				vertices.resize(xRange,yRange);

				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						cellTypes(i,j) = REGULAR;
						faceTypes(i,j) = REGULAR;
						faceTypes(i,j)(B) = COVERED;
						volumeFractions(i,j) = 1;
						areaFractions(i,j) = 1;
						areaFractions(i,j)(B) = 0;
						outwardNormals(i,j) = 0;

						centers(i,j)(0) = (i-1)*h + (h/2) - offset;
						centers(i,j)(1) = (j-1)*h + (h/2) - offset;

						centroids(i,j)(N)(0) = centers(i,j)(0);
						centroids(i,j)(N)(1) = getCellEdge(i,j,N);

						centroids(i,j)(S)(0) = centers(i,j)(0);
						centroids(i,j)(S)(1) = getCellEdge(i,j,S);

						centroids(i,j)(E)(0) = getCellEdge(i,j,E);
						centroids(i,j)(E)(1) = centers(i,j)(1);

						centroids(i,j)(W)(0) = getCellEdge(i,j,W);
						centroids(i,j)(W)(1) = centers(i,j)(1);

						if(!contains(centers(i,j))
								&&
								((getCellEdge(i,j,N) < 0 || getCellEdge(i,j,N) > 1)
								&& (getCellEdge(i,j,S) < 0 || getCellEdge(i,j,S) > 1))
								||
								((getCellEdge(i,j,E) < 0 || getCellEdge(i,j,E) > 1)
								&& (getCellEdge(i,j,W) < 0 || getCellEdge(i,j,W) > 1))){
							cellTypes(i,j) = COVERED;
							faceTypes(i,j) = COVERED;
							volumeFractions(i,j) = 0;
							areaFractions(i,j) = 0;
							outwardNormals(i,j) = 0;
						}
						else{
							if(abs(centers(i,j)(0)) < h/2){
								if(abs(centers(i,j)(1)) < h/2){ // i = 1, j = 1
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(S) = COVERED;
									areaFractions(i,j)(S) = 0;
									centroids(i,j)(S) = 0;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(E) = IRREGULAR;
									centroids(i,j)(E)(0) = getCellEdge(i,j,E);
									centroids(i,j)(E)(1) = getCellEdge(i,j,N)/2;
									areaFractions(i,j)(E) = getCellEdge(i,j,N)/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(N)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(E)(1);
									areaFractions(i,j)(B) = sqrt(pow2(areaFractions(i,j)(N)) + pow2(areaFractions(i,j)(E)));
									outwardNormals(i,j)(0) = -areaFractions(i,j)(E)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = -areaFractions(i,j)(N)/areaFractions(i,j)(B);
								}
								else if(abs(centers(i,j)(1) - 1) < h/2){ // i = 1, j = jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(N) = COVERED;
									areaFractions(i,j)(N) = 0;
									centroids(i,j)(N) = 0;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(E) = IRREGULAR;
									centroids(i,j)(E)(0) = getCellEdge(i,j,E);
									centroids(i,j)(E)(1) = (1 + getCellEdge(i,j,S))/2;
									areaFractions(i,j)(E) = (1 - getCellEdge(i,j,S))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(S)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(E)(1);
									areaFractions(i,j)(B) = sqrt(pow2(areaFractions(i,j)(S)) + pow2(areaFractions(i,j)(E)));
									outwardNormals(i,j)(0) = -areaFractions(i,j)(E)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = areaFractions(i,j)(S)/areaFractions(i,j)(B);
								}
								else{ // i = 1, j != 1, j != jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(W) = COVERED;
									areaFractions(i,j)(W) = 0;
									centroids(i,j)(W) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = getCellEdge(i,j,E)/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = getCellEdge(i,j,E)/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = 0;
									centroids(i,j)(B)(1) = centers(i,j)(1);
									areaFractions(i,j)(B) = 1;
									outwardNormals(i,j)(0) = -1;
									outwardNormals(i,j)(1) = 0;
								}
							}
							else if(abs(centers(i,j)(0) - 1) < h/2){
								if(abs(centers(i,j)(1)) < h/2){ // i = iMax, j = 1
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(S) = COVERED;
									areaFractions(i,j)(S) = 0;
									centroids(i,j)(S) = 0;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(W) = IRREGULAR;
									centroids(i,j)(W)(0) = getCellEdge(i,j,W);
									centroids(i,j)(W)(1) = getCellEdge(i,j,N)/2;
									areaFractions(i,j)(W) = getCellEdge(i,j,N)/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(N)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(W)(1);
									areaFractions(i,j)(B) = sqrt(pow2(areaFractions(i,j)(N)) + pow2(areaFractions(i,j)(W)));
									outwardNormals(i,j)(0) = areaFractions(i,j)(W)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = -areaFractions(i,j)(N)/areaFractions(i,j)(B);
								}
								else if(abs(centers(i,j)(1) - 1) < h/2){ // i = iMax, j = jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(N) = COVERED;
									areaFractions(i,j)(N) = 0;
									centroids(i,j)(N) = 0;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(W) = IRREGULAR;
									centroids(i,j)(W)(0) = getCellEdge(i,j,W);
									centroids(i,j)(W)(1) = (1 + getCellEdge(i,j,S))/2;
									areaFractions(i,j)(W) = (1 - getCellEdge(i,j,S))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = centroids(i,j)(S)(0);
									centroids(i,j)(B)(1) = centroids(i,j)(W)(1);
									areaFractions(i,j)(B) = sqrt(pow2(areaFractions(i,j)(S)) + pow2(areaFractions(i,j)(W)));
									outwardNormals(i,j)(0) = areaFractions(i,j)(W)/areaFractions(i,j)(B);
									outwardNormals(i,j)(1) = areaFractions(i,j)(S)/areaFractions(i,j)(B);
								}
								else{ // i = iMax, j != 1, j != jMax
									cellTypes(i,j) = IRREGULAR;

									faceTypes(i,j)(E) = COVERED;
									areaFractions(i,j)(E) = 0;
									centroids(i,j)(E) = 0;

									faceTypes(i,j)(N) = IRREGULAR;
									centroids(i,j)(N)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(N)(1) = getCellEdge(i,j,N);
									areaFractions(i,j)(N) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(S) = IRREGULAR;
									centroids(i,j)(S)(0) = (1 + getCellEdge(i,j,W))/2;
									centroids(i,j)(S)(1) = getCellEdge(i,j,S);
									areaFractions(i,j)(S) = (1 - getCellEdge(i,j,W))/h;

									faceTypes(i,j)(B) = IRREGULAR;
									centroids(i,j)(B)(0) = 1;
									centroids(i,j)(B)(1) = centers(i,j)(1);
									areaFractions(i,j)(B) = 1;
									outwardNormals(i,j)(0) = 1;
									outwardNormals(i,j)(1) = 0;
								}
							}
							else if(abs(centers(i,j)(1)) < h/2){ // j = 1, i != 1, i != iMax
								cellTypes(i,j) = IRREGULAR;

								faceTypes(i,j)(S) = COVERED;
								areaFractions(i,j)(S) = 0;
								centroids(i,j)(S) = 0;

								faceTypes(i,j)(E) = IRREGULAR;
								centroids(i,j)(E)(0) = getCellEdge(i,j,E);
								centroids(i,j)(E)(1) = getCellEdge(i,j,N)/2;
								areaFractions(i,j)(E) = getCellEdge(i,j,N)/h;

								faceTypes(i,j)(W) = IRREGULAR;
								centroids(i,j)(W)(0) = getCellEdge(i,j,W);
								centroids(i,j)(W)(1) = getCellEdge(i,j,N)/2;
								areaFractions(i,j)(W) = getCellEdge(i,j,N)/h;

								faceTypes(i,j)(B) = IRREGULAR;
								centroids(i,j)(B)(0) = centers(i,j)(0);
								centroids(i,j)(B)(1) = 0;
								areaFractions(i,j)(B) = 1;
								outwardNormals(i,j)(0) = 0;
								outwardNormals(i,j)(1) = -1;
							}
							else if(abs(centers(i,j)(1) - 1) < h/2){ // j = jMax, i != 1, i != iMax
								cellTypes(i,j) = IRREGULAR;

								faceTypes(i,j)(N) = COVERED;
								areaFractions(i,j)(N) = 0;
								centroids(i,j)(N) = 0;

								faceTypes(i,j)(E) = IRREGULAR;
								centroids(i,j)(E)(0) = getCellEdge(i,j,E);
								centroids(i,j)(E)(1) = (1 + getCellEdge(i,j,S))/2;
								areaFractions(i,j)(E) = (1 - getCellEdge(i,j,S))/h;

								faceTypes(i,j)(W) = IRREGULAR;
								centroids(i,j)(W)(0) = getCellEdge(i,j,W);
								centroids(i,j)(W)(1) = (1 + getCellEdge(i,j,S))/2;
								areaFractions(i,j)(W) = (1 - getCellEdge(i,j,S))/h;

								faceTypes(i,j)(B) = IRREGULAR;
								centroids(i,j)(B)(0) = centers(i,j)(0);
								centroids(i,j)(B)(1) = 1;
								areaFractions(i,j)(B) = 1;
								outwardNormals(i,j)(0) = 0;
								outwardNormals(i,j)(1) = 1;
							}
						}
					}
				}
				for(int i = 1; i <= iMax; i++){
					for(int j = 1; j <= jMax; j++){
						volumeFractions(i,j) = calculateVolumeFraction(i,j);
						if(abs(centers(i,j)(0)) < h/2 && abs(centers(i,j)(1)) < h/2){
							volumeFractions(i,j) *= 2.0;
						}
						if(abs(centers(i,j)(0)) < h/2 && abs(1 - centers(i,j)(1)) < h/2){
							volumeFractions(i,j) *= 2.0;
						}
						if(abs(1 - centers(i,j)(0)) < h/2 && abs(centers(i,j)(1)) < h/2){
							volumeFractions(i,j) *= 2.0;
						}
						if(abs(1 - centers(i,j)(0)) < h/2 && abs(1 - centers(i,j)(1)) < h/2){
							volumeFractions(i,j) *= 2.0;
						}
						countVertices(i,j);
						getVertices(i,j);
						cellCentroids(i,j) = calculateCellCentroid(i,j);
					}
				}
			}
		public:
			double offset;
			bool contains(Coord c){
				double x = c(0), y = c(1);
				if(x >= 0 && x <= 1 && y >= 0 && y <= 1)
					return true;
				return false;
			}
			UnitSquare(double h, double offset){
				this->h = h;
				this->offset = offset;
				make();
			}
			UnitSquare(const UnitSquare & u){
				this->h = u.h;
				this->offset = u.offset;
				make();
				this->volumeFractions = u.volumeFractions;
			}
			Grid * copy() const{
				return new UnitSquare(*this);
			}
			void respace(double hNew){
				this->h = hNew;
				make();
			}
		};

		class EBUtilities{
		public:
			/*static list<Coefficient> getGradientCoefficientList(CellDoubleArray &u, int i, int j, Direction dir, Grid *g){
				list<Coefficient> out;
				double h = g->h;
				return out;
			}*/

			static TinyVector<double,9> getGradientCoefficients(CellDoubleArray &u, int i, int j, Direction dir, Grid *g){
				return getGradientCoefficients(i,j,dir,g);
			}
			static Coefficients getGradientCoefficients(int i, int j, Direction dir, Grid * g){
				TinyVector<double,9> out;
				double h = g->h;
				out = 0;
				if(dir != B){
					if(g->faceTypes(i,j)(dir) == REGULAR){
						if(dir == N){
							out(N) = 1/h;
							out(C) = -1/h;
						}
						else if(dir == S){
							out(C) = 1/h;
							out(S) = -1/h;
						}
						else if(dir == E){
							out(E) = 1/h;
							out(C) = -1/h;
						}
						else if(dir == W){
							out(C) = 1/h;
							out(W) = -1/h;
						}
					}
					else if(g->faceTypes(i,j)(dir) == IRREGULAR){
						double a = g->areaFractions(i,j)(dir);
						if(dir == N){
							if(g->isUncovered(i+1,j) && g->isUncovered(i+1,j+1)){
								out(N)  =  ((1 + a)/2)/h;
								out(C)  = -((1 + a)/2)/h;

								out(NE) =  ((1 - a)/2)/h;
								out(E)  = -((1 - a)/2)/h;
							}
							else if(g->isUncovered(i-1,j) && g->isUncovered(i-1,j+1)){
								out(N)  =  ((1 + a)/2)/h;
								out(C)  = -((1 + a)/2)/h;

								out(NW) =  ((1 - a)/2)/h;
								out(W)  = -((1 - a)/2)/h;
							}
							else{
								cout << "Problem" << endl;
							}
						}
						else if(dir == S){
							if(g->isUncovered(i+1,j) && g->isUncovered(i+1,j-1)){
								out(C)  =  ((1 + a)/2)/h;
								out(S)  = -((1 + a)/2)/h;

								out(E)  =  ((1 - a)/2)/h;
								out(SE) = -((1 - a)/2)/h;
							}
							else if(g->isUncovered(i-1,j) && g->isUncovered(i-1,j-1)){
								out(C)  =  ((1 + a)/2)/h;
								out(S)  = -((1 + a)/2)/h;

								out(W)  =  ((1 - a)/2)/h;
								out(SW) = -((1 - a)/2)/h;
							}
							else{
								cout << "Problem" << endl;
							}
						}
						else if(dir == E){
							if(g->isUncovered(i,j+1) && g->isUncovered(i+1,j+1)){
								out(E)  =  ((1 + a)/2)/h;
								out(C)  = -((1 + a)/2)/h;

								out(NE) =  ((1 - a)/2)/h;
								out(N)  = -((1 - a)/2)/h;
							}
							else if(g->isUncovered(i,j-1) && g->isUncovered(i+1,j-1)){
								out(E)  =  ((1 + a)/2)/h;
								out(C)  = -((1 + a)/2)/h;

								out(SE) =  ((1 - a)/2)/h;
								out(S)  = -((1 - a)/2)/h;
							}
							else{
								cout << "Problem" << endl;
							}

						}
						else if(dir == W){
							if(g->isUncovered(i,j+1) && g->isUncovered(i-1,j+1)){
								out(C)  =  ((1 + a)/2)/h;
								out(W)  = -((1 + a)/2)/h;

								out(N)  =  ((1 - a)/2)/h;
								out(NW) = -((1 - a)/2)/h;
							}
							else if(g->isUncovered(i,j-1) && g->isUncovered(i-1,j-1)){
								out(C)  =  ((1 + a)/2)/h;
								out(W)  = -((1 + a)/2)/h;

								out(S)  =  ((1 - a)/2)/h;
								out(SW) = -((1 - a)/2)/h;
							}
							else{
								cout << "Problem" << endl;
							}
						}
					}
				}
				return out;
			}
			static double getFaceGradient(CellDoubleArray &u, int i, int j, Direction dir, Grid *g){
				double h = g->h;
				double grad;
				if(dir != B){
					if(g->faceTypes(i,j)(dir) == REGULAR){
						if(dir == N){
							return (u(i,j+1) - u(i,j))/h;
						}
						else if(dir == S){
							return (u(i,j) - u(i,j-1))/h;
						}
						else if(dir == E){
							return (u(i+1,j) - u(i,j))/h;
						}
						else if(dir == W){
							return (u(i,j) - u(i-1,j))/h;
						}
					}
					else if(g->faceTypes(i,j)(dir) == IRREGULAR){
						double a = g->areaFractions(i,j)(dir);
						if(dir == N){
							if(g->isUncovered(i+1,j) && g->isUncovered(i+1,j+1)){
								grad = ((1 + a)/2)*(u(i,j+1) - u(i,j))/h + ((1 - a)/2)*(u(i+1,j+1) - u(i+1,j))/h;
							}
							else if(g->isUncovered(i-1,j) && g->isUncovered(i-1,j+1)){
								grad = ((1 + a)/2)*(u(i,j+1) - u(i,j))/h + ((1-a)/2)*(u(i-1,j+1) - u(i-1,j))/h;
							}
						}
						else if(dir == S){
							if(g->isUncovered(i+1,j) && g->isUncovered(i+1,j-1)){
								grad = ((1 + a)/2)*(u(i,j) - u(i,j-1))/h + ((1 - a)/2)*(u(i+1,j) - u(i+1,j-1))/h;
							}
							else if(g->isUncovered(i-1,j) && g->isUncovered(i-1,j-1)){
								grad = ((1 + a)/2)*(u(i,j) - u(i,j-1))/h + ((1 - a)/2)*(u(i-1,j) - u(i-1,j-1))/h;
							}
						}
						else if(dir == E){

							if(g->isUncovered(i,j+1) && g->isUncovered(i+1,j+1)){
								grad = ((1 + a)/2)*(u(i+1,j) - u(i,j))/h + ((1 - a)/2)*(u(i+1,j+1) - u(i,j+1))/h;
							}
							else if(g->isUncovered(i,j-1) && g->isUncovered(i+1,j-1)){
								grad = ((1 + a)/2)*(u(i+1,j) - u(i,j))/h + ((1 - a)/2)*(u(i+1,j-1) - u(i,j-1))/h;
							}

						}
						else if(dir == W){
							if(g->isUncovered(i,j+1) && g->isUncovered(i-1,j+1)){
								grad = ((1 + a)/2)*(u(i,j) - u(i-1,j))/h + ((1 - a)/2)*(u(i,j+1) - u(i-1,j+1))/h;
							}
							else if(g->isUncovered(i,j-1) && g->isUncovered(i-1,j-1)){
								grad = ((1 + a)/2)*(u(i,j) - u(i-1,j))/h + ((1 - a)/2)*(u(i,j-1) - u(i-1,j-1))/h;
							}
						}
					}
					else{ // Tried to get gradient through a covered face.
						grad = 0;
					}
				}
				else{ /* Gradient at the boundary is currently unimplemented, since at
					   * present we only need to solve problems with no flux boundary
					   * conditions.
					   */

					cout << "Tried to find boundary gradient, which is currently unimplemented." << endl;

				}
				return grad;
			}
		};

	}
}

#endif /* GEOMETRY_H_ */
