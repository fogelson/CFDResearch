/*
 * Geometry.h
 *
 *  Created on: Jan 29, 2011
 *      Author: bfogelson
 */

/* This file will define classes and methods related to
 * the physical and numerical domains on which we are solving
 * PDEs.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <math.h>
#include <iostream>
#include <list>

using namespace std;
using namespace blitz;

namespace CFD{
	namespace GeometryOld{
		/* The typedef command is used to make an additional
		 * name for a certain type of object. For instance,
		 * the command
		 * 		typedef TinyVector<double,2> Coord;
		 * tells the compiler that every time it sees Coord,
		 * to treat it as a TinyVector with 2 entries, each
		 * of which is a double.
		 *
		 * There are two main reasons for using the typedef
		 * command. First, to come up with a shorter way of
		 * calling a commonly used object. It is a whole lot
		 * easier to type "Coord", than it is to type
		 * "TinyVector<double,2> Coord.
		 *
		 * Second, defining types like this lets you tell the
		 * differences between objects that are actually of a
		 * different type. For instance, I define the types
		 * "Direction", and "CoordType", both of which are
		 * actually just ints. By referring to them in the rest
		 * of my code by the names Direction and CoordType, it
		 * is easy to keep track of what I'm actually doing. I
		 * am much less likely to do something like try to compare
		 * a CoordType (such as INTERIOR or EXTERIOR) to a
		 * Direction (NORTH or SOUTH) than if they were all just
		 * labeled as ints.
		 */

		// Coord (short for Coordinate) is an (x,y) ordered pair.
		typedef TinyVector<double,2> Coord;
		// An array of Coords.
		typedef Array<Coord,2> CoordArray;

		// Spacing is a vector of distances in the N, S, E, W
		// directions. Think of it as the distance from a grid point
		// to its nearest neighbor or the domain boundary in a
		// given direction.
		typedef TinyVector<double,4> Spacing;
		// An array os Spacings.
		typedef Array<Spacing,2> SpacingArray;
		// An int giving the particular direction. Declared below.
		typedef int Direction;
		// NeighborScalar is a vector of values in the N, S, E, W directions.
		// It should be used in functions that return the neighboring
		// values of a grid scalar, be those values at neighboring
		// grid points or at boundary points.
		typedef TinyVector<double,4> NeighborScalar;

		// An int giving a type for a Coord. For instance, could
		// be EXTERIOR, REGULAR, IRREGULAR, or something else.
		// Declared below.
		typedef int CoordType;
		// An array of CoordTypes.
		typedef Array<CoordType,2> CoordTypeArray;

		// An scalar valued array of doubles which should be though
		// of as a scalar valued grid function.
		typedef Array<double,2> GridScalar;

		// An int giving a type for a domain Boundary, e.g. Dirichlet
		// or Neumann. Declared below.
		typedef int BoundaryType;

		/* CoordTypes that might be useful. Currently, the UNVISITED
		 * flag is unused, but may become useful for iterating through
		 * an array in some particular manner (e.g. breadth first
		 * search). BOUNDARY is also unused. INTERIOR is currently only
		 * used when initializing an array of Coords, since those
		 * Coords just eventually be marked as REGULAR or IRREGULAR,
		 * and not just as INTERIOR. EXTERIOR is used to flag anything
		 * outside the domain of interest.
		 */
		const CoordType UNVISITED = -1;
		const CoordType EXTERIOR = 0;
		const CoordType IRREGULAR = 1;
		const CoordType REGULAR = 2;
		const CoordType BOUNDARY = 3;
		const CoordType INTERIOR = 4;

		/* Cardinal directions. These are used by the typedef
		 * Spacing to index elements in that TinyVector. For
		 * instance, if sp is a Spacing, then sp(NORTH) should
		 * be a double that gives the distance from a gridpoint
		 * to its northern neighbor (if it is a regular point)
		 * or to the northern boundary of the domain (if it is
		 * irregular).
		 *
		 * NOTE: The numbers assigned to NORTH, SOUTH, EAST, and
		 * WEST are not entirely arbitrary. Since they are used
		 * as indices to a four element TinyVector, they must be
		 * the integers 0, 1, 2, and 3. Which number corresponds
		 * with which direction is, however, arbitrary.
		 */
		const Direction NORTH = 0;
		const Direction SOUTH = 1;
		const Direction EAST = 2;
		const Direction WEST = 3;
		const Direction NORTHEAST = 4;
		const Direction SOUTHEAST = 5;
		const Direction NORTHWEST = 6;
		const Direction SOUTHWEST = 7;

		const BoundaryType DIRICHLET = 0;
		const BoundaryType NEUMANN = 1;
		const BoundaryType NOFLUX = 2;

		typedef int CellType;
		const CellType REGULARCELL = 0;
		const CellType IRREGULARCELL = 1;
		const CellType EXTERIORCELL = 2;


		// Forward declaration of classes.
		class Boundary;
		class BoundaryCondition;
		class Grid;
		class Circle;
		class GridOperator;
		class Cell;

		/* The BoundaryCondition class gives an object to handle boundaries
		 * of domains. It stores what type of BoundaryType it is, although
		 * it is currently only implemented for Dirichlet boundaries.
		 */

		double zeroBoundaryFunction(double x, double y){
			return 0;
		}

		class Boundary{
		public:
			virtual ~Boundary(){}
			virtual bool contains(Coord c) = 0;
			virtual bool contains(double x, double y){
				Coord c(x,y);
				return contains(c);
			}
			virtual Coord nearestBoundary(Coord c, Direction d) = 0;
		};

		class CircleBoundary : public Boundary{
			Coord c;
			double r;
		public:
			CircleBoundary(Coord c, double r){
				this->c = c;
				this->r = r;
			}
			virtual ~CircleBoundary(){}
			virtual bool contains(Coord p){
				return (pow2(c(0) - p(0)) + pow2(c(1) - p(1)) < pow2(r));
			}
			virtual Coord nearestBoundary(Coord p, Direction d){
				double xC = c(0), yC = c(1), xP = p(0), yP = p(1), xB, yB;
				if(d == NORTH){
					xB = xP;
					yB = sqrt(pow2(r) - pow2(xB - xC)) + yC;
				}
				else if(d == SOUTH){
					xB = xP;
					yB = -sqrt(pow2(r) - pow2(xB-xC)) + yC;
				}
				else if(d == EAST){
					yB = yP;
					xB = sqrt(pow2(r) - pow2(yB - yC)) + xC;
				}
				else if(d == WEST){
					yB = yP;
					xB = -sqrt(pow2(r) - pow2(yB - yC)) + xC;
				}
				Coord b(xB,yB);
				return b;
			}
		};

		class BoundaryCondition{
			BoundaryType type;
			double (*func)(double, double);
		public:
			BoundaryCondition(){
				setType(DIRICHLET);
				setF((*zeroBoundaryFunction));
			}
			BoundaryCondition(BoundaryType type){
				setType(type);
				setF((*zeroBoundaryFunction));
			}
			BoundaryCondition(BoundaryType type, double (*f)(double,double)){
				setType(type);
				setF(f);
			}
			BoundaryCondition(double (*f)(double,double)){
				setType(DIRICHLET);
				setF(f);
			}
			void setType(BoundaryType type){
				this->type = type;
			}
			BoundaryType getType(){
				return type;
			}
			void setF(double (*f)(double,double)){
				func = f;
			}
			double operator ()(double x, double y){
				return (*func)(x,y);
			}
			double operator ()(Coord c){
				return (*func)(c(0),c(1));
			}
		};

		class Grid{
		public:
			virtual ~Grid(){}
			virtual bool contains(Coord c) = 0;
			virtual double getR() = 0;
			virtual double getH() = 0;
			virtual CoordArray getPoints() = 0;
			virtual CoordTypeArray getTypes() = 0;
			virtual SpacingArray getSpacings() = 0;
			virtual GridScalar makeScalar() = 0;
			virtual CoordType getType(int i, int j) = 0;
			virtual Spacing getSpacing(int i, int j) = 0;
			virtual Coord getCoord(int i, int j) = 0;
			virtual double getX(int i, int j) = 0;
			virtual double getY(int i, int j) = 0;
			virtual GridScalar getXes() = 0;
			virtual GridScalar getYes() = 0;
			virtual Coord nearestBoundary(int i, int j, Direction d) = 0;
			virtual double evaluateOnNearestBoundary(int i, int j, Direction d, double (*boundary)(double,double)) = 0;
			virtual int interiorPoints() = 0;
			virtual void setBC(BoundaryCondition bc) = 0;
			virtual BoundaryCondition getBC() = 0;
			virtual double boundaryValue(int i, int j, Direction d) = 0;
			virtual NeighborScalar getNeighbors(GridScalar u, int i, int j) = 0;
			virtual int getM() = 0;
		};
		class ACircle : public Grid{
			Coord center; // Center of the circle
			CoordArray points; // Array of points on the grid
			CoordTypeArray types; // Type of each point on the grid
			SpacingArray spacings; // Spacing from each gridpoint to neighbors
			Range xRange, yRange; // Range of array indices in x and y directions
			double r, h; // Radius and grid spacing
			BoundaryCondition bc;

			// M is the number of entries in the array from the center out to
			// one edge (so the actual arrays are all 2M+1 by 2M+1. N is the
			// number of those points that are interior points (either REGULAR
			// or IRREGULAR).
			int M, N;
			// Returns true if the given Coord is in the interior of the circle.
		public:
			virtual bool contains(Coord p){
				return (pow2(p(0) - center(0)) + pow2(p(1) - center(1)) < pow2(r));
			}
		private:
			// Initializes grid of points, labels them as REGULAR, IRREGULAR, or EXTERIOR.
			void make(){
				M = ceil(r/h); // Number of points needed to get just beyond circle boundary.
				xRange.setRange(-M,M,1);
				yRange.setRange(-M,M,1);
				// We have to make sure all our arrays are the right size, which
				// means using the resize() function.
				points.resize(xRange,yRange);
				types.resize(xRange,yRange);
				spacings.resize(xRange,yRange);
				// Default type for any point is EXTERIOR. If it isn't, we will find out
				// soon enough.
				types = EXTERIOR;
				// Iterate over the grid, and set the (x,y) coordinates of each point
				// on the grid.
				for(int i = -M; i <= M; i++){
					double x = i*h + center(0);
					for(int j = -M; j <= M; j++){
						double y = j*h + center(1);
						Coord p(x,y);
						points(i,j) = p;
						if(contains(p)){
							// If it is inside the circle, it is an interior point.
							types(i,j) = INTERIOR;
						}
					}
				}
				// Count up the interior points.
				N = sum(where(types == INTERIOR, 1, 0));
				/* Iterate over the grid again. At each interior point, check whether
				 * all its neighbors are also interior. If it is, change its type to
				 * REGULAR. If not, change its type to EXTERIOR and set its Spacing
				 * in the appropriate direction to be the distance to the boundary,
				 * not h.
				 */
				for(int i = -M; i <= M; i++){
					for(int j = -M; j <= M; j++){
						if(types(i,j) == INTERIOR){
							double x = points(i,j)(0), y = points(i,j)(1);
							// Spacings in N, S, E, and W directions default to
							// the grid spacing. If the point turns out to be
							// irregular, we will have to compute the distance to
							// the boundary.
							double hN = h, hS = h, hE = h, hW = h;
							// Set type to REGULAR, then check if any neighbors are
							// EXTERIOR. If any of them are, then types(i,j) will be
							// reset to IRREGULAR.
							types(i,j) = REGULAR;
							// Check west
							if(types(i-1,j) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double xB2 = pow2(r) - pow2(y - center(1));
								double xB = -sqrt(xB2) + center(0);
								hW = x - xB;
							}
							// Check east
							if(types(i+1,j) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double xB2 = pow2(r) - pow2(y - center(1));
								double xB = sqrt(xB2) + center(0);
								hE = xB - x;
							}
							// Check south
							if(types(i,j-1) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double yB2 = pow2(r) - pow2(x - center(0));
								double yB = -sqrt(yB2) + center(1);
								hS = y - yB;
							}
							// Check north
							if(types(i,j+1) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double yB2 = pow2(r) - pow2(x - center(0));
								double yB = sqrt(yB2) + center(1);
								hN = yB - y;
							}
							if(types(i+1,j+1) == EXTERIOR){
								types(i,j) = IRREGULAR;
							}
							if(types(i+1,j-1) == EXTERIOR){
								types(i,j) = IRREGULAR;
							}
							if(types(i-1,j+1) == EXTERIOR){
								types(i,j) = IRREGULAR;
							}
							if(types(i-1,j-1) == EXTERIOR){
								types(i,j) = IRREGULAR;
							}
							// Set spacings
							spacings(i,j)(NORTH) = hN;
							spacings(i,j)(SOUTH) = hS;
							spacings(i,j)(EAST) = hE;
							spacings(i,j)(WEST) = hW;
						}
					}
				}
			}
		public:
			virtual ~ACircle(){}
			/* Constructors for the circle. If no center is given,
			 * default to being centered at (0,0).
			 */
			ACircle(double r, double h, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				center = Coord(0,0);
				setBC(bc);
				make();
			}
			ACircle(double r, double h){
				this->r = r;
				this->h = h;
				center = Coord(0,0);
				make();
			}
			ACircle(double r, double h, Coord center, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				this->center = center;
				setBC(bc);
				make();
			}
			ACircle(double r, double h, Coord center){
				this->r = r;
				this->h = h;
				this->center = center;
				make();
			}
			ACircle(double r, double h, double x0, double y0, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				center = Coord(x0,y0);
				setBC(bc);
				make();
			}
			ACircle(double r, double h, double x0, double y0){
				this->r = r;
				this->h = h;
				center = Coord(x0,y0);
				make();
			}
			/* Many of these methods are pretty self-explanatory. The
			 * getPoints(), getTypes(), and getSpacings() methods may
			 * want to be avoided, since they involve copying a whole
			 * array. Rather, use the getPoint(i,j), getType(i,j), and
			 * getSpacing(i,j) methods to only see what the Coord,
			 * CoordType, or Spacing at a particular grid point is, when
			 * needed elsewhere in the code.
			 */
			virtual double getR(){
				return r;
			}
			virtual double getH(){
				return h;
			}
			virtual Coord getCenter(){
				return center;
			}
			virtual CoordArray getPoints(){
				return points.copy();
			}
			virtual CoordTypeArray getTypes(){
				return types.copy();
			}
			virtual SpacingArray getSpacings(){
				return spacings.copy();
			}
			/* Creates a GridScalar object of the appropriate
			 * size and indexing to represent a grid function
			 * on this particular grid, and initializes it to
			 * 0 everywhere.
			 */
			virtual GridScalar makeScalar(){
				GridScalar u(xRange,yRange);
				u = 0;
				return u;
			}
			virtual CoordType getType(int i, int j){
				return types(i,j);
			}
			virtual Spacing getSpacing(int i, int j){
				return spacings(i,j);
			}
			virtual Coord getCoord(int i, int j){
				return points(i,j);
			}
			virtual double getX(int i, int j){
				return points(i,j)(0);
			}
			virtual double getY(int i, int j){
				return points(i,j)(1);
			}
			virtual GridScalar getXes(){
				return points.extractComponent(double(),0,2);
			}
			virtual GridScalar getYes(){
				return points.extractComponent(double(),1,2);
			}
			virtual Coord nearestBoundary(int i, int j, Direction d){
				double x = getX(i,j), y = getY(i,j);
				double xB, yB;
				if(d == NORTH){
					xB = x;
					yB = sqrt(pow2(r) - pow2(x - center(0))) + center(1);
				}
				else if(d == SOUTH){
					xB = x;
					yB = -sqrt(pow2(r) - pow2(x - center(0))) + center(1);
				}
				else if(d == EAST){
					yB = y;
					xB = sqrt(pow2(r) - pow2(y - center(1))) + center(0);
				}
				else if(d == WEST){
					yB = y;
					xB = -sqrt(pow2(r) - pow2(y - center(1))) + center(0);
				}
				return Coord(xB,yB);
			}
			virtual double evaluateOnNearestBoundary(int i, int j, Direction d, double (*boundary)(double,double)){
				Coord b = nearestBoundary(i,j,d);
				return (*boundary)(b(0),b(1));
			}
			virtual int interiorPoints(){
				return N;
			}
			virtual void setBC(BoundaryCondition bc){
				this->bc = bc;
			}
			virtual BoundaryCondition getBC(){
				return bc;
			}
			virtual double boundaryValue(int i, int j, Direction d){
				Coord b = nearestBoundary(i,j,d);
				return bc(b);
			}
			virtual NeighborScalar getNeighbors(GridScalar u, int i, int j){
				double uN, uS, uE, uW;
				if(getType(i,j+1) == EXTERIOR){
					uN = boundaryValue(i,j,NORTH);
				}
				else{
					uN = u(i,j+1);
				}
				if(getType(i,j-1) == EXTERIOR){
					uS = boundaryValue(i,j,SOUTH);
				}
				else{
					uS = u(i,j-1);
				}
				if(getType(i+1,j) == EXTERIOR){
					uE = boundaryValue(i,j,EAST);
				}
				else{
					uE = u(i+1,j);
				}
				if(getType(i-1,j) == EXTERIOR){
					uW = boundaryValue(i,j,WEST);
				}
				else{
					uW = u(i-1,j);
				}
				NeighborScalar out;
				out(NORTH) = uN;
				out(SOUTH) = uS;
				out(EAST) = uE;
				out(WEST) = uW;
				return out;
			}
			virtual int getM(){
				return M;
			}
		};
		/* The circle class stores the center and radius of a circle in
		 * the physical domain, a numerical grid spacing, and all
		 * information about the grid, such as coordinates of gridpoints,
		 * whether those points are interior or exterior, and the like.
		 *
		 * The grid is always centered with a point at the center of the
		 * circle. Grid arrays are always indexed from -M to +M, with
		 * (0,0) indexing the center of the circle.
		 */
		class Circle : public Grid{
			Coord center; // Center of the circle
			CoordArray points; // Array of points on the grid
			CoordTypeArray types; // Type of each point on the grid
			SpacingArray spacings; // Spacing from each gridpoint to neighbors
			Range xRange, yRange; // Range of array indices in x and y directions
			double r, h; // Radius and grid spacing
			BoundaryCondition bc;
			// M is the number of entries in the array from the center out to
			// one edge (so the actual arrays are all 2M+1 by 2M+1. N is the
			// number of those points that are interior points (either REGULAR
			// or IRREGULAR).
			int M, N;
		public:
			// Returns true if the given Coord is in the interior of the circle.
			virtual bool contains(Coord p){
				return (abs(p(0)) < r && abs(p(1)) < r);
//				return (pow2(p(0) - center(0)) + pow2(p(1) - center(1)) < pow2(r));
			}
		private:
			// Initializes grid of points, labels them as REGULAR, IRREGULAR, or EXTERIOR.
			void make(){
				M = ceil(r/h); // Number of points needed to get just beyond circle boundary.
				xRange.setRange(-M,M,1);
				yRange.setRange(-M,M,1);
				// We have to make sure all our arrays are the right size, which
				// means using the resize() function.
				points.resize(xRange,yRange);
				types.resize(xRange,yRange);
				spacings.resize(xRange,yRange);
				// Default type for any point is EXTERIOR. If it isn't, we will find out
				// soon enough.
				types = EXTERIOR;
				// Iterate over the grid, and set the (x,y) coordinates of each point
				// on the grid.
				for(int i = -M; i <= M; i++){
					double x = i*h + center(0);
					for(int j = -M; j <= M; j++){
						double y = j*h + center(1);
						Coord p(x,y);
						points(i,j) = p;
						if(contains(p)){
							// If it is inside the circle, it is an interior point.
							types(i,j) = INTERIOR;
						}
					}
				}
				// Count up the interior points.
				N = sum(where(types == INTERIOR, 1, 0));
				/* Iterate over the grid again. At each interior point, check whether
				 * all its neighbors are also interior. If it is, change its type to
				 * REGULAR. If not, change its type to EXTERIOR and set its Spacing
				 * in the appropriate direction to be the distance to the boundary,
				 * not h.
				 */
				for(int i = -M; i <= M; i++){
					for(int j = -M; j <= M; j++){
						if(types(i,j) == INTERIOR){
							double x = points(i,j)(0), y = points(i,j)(1);
							// Spacings in N, S, E, and W directions default to
							// the grid spacing. If the point turns out to be
							// irregular, we will have to compute the distance to
							// the boundary.
							double hN = h, hS = h, hE = h, hW = h;
							// Set type to REGULAR, then check if any neighbors are
							// EXTERIOR. If any of them are, then types(i,j) will be
							// reset to IRREGULAR.
							types(i,j) = REGULAR;
							// Check west
							if(types(i-1,j) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double xB = -r + center(0);
								hW = x - xB;
							}
							// Check east
							if(types(i+1,j) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double xB = r + center(0);
								hE = xB - x;
							}
							// Check south
							if(types(i,j-1) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double yB = -r + center(1);
								hS = y - yB;
							}
							// Check north
							if(types(i,j+1) == EXTERIOR){
								types(i,j) = IRREGULAR;
								double yB = r + center(1);
								hN = yB - y;
							}
							// Set spacings
							spacings(i,j)(NORTH) = hN;
							spacings(i,j)(SOUTH) = hS;
							spacings(i,j)(EAST) = hE;
							spacings(i,j)(WEST) = hW;
						}
					}
				}
			}
		public:
			virtual ~Circle(){}
			/* Constructors for the circle. If no center is given,
			 * default to being centered at (0,0).
			 */
			Circle(double r, double h, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				center = Coord(0,0);
				setBC(bc);
				make();
			}
			Circle(double r, double h){
				this->r = r;
				this->h = h;
				center = Coord(0,0);
				make();
			}
			Circle(double r, double h, Coord center, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				this->center = center;
				setBC(bc);
				make();
			}
			Circle(double r, double h, Coord center){
				this->r = r;
				this->h = h;
				this->center = center;
				make();
			}
			Circle(double r, double h, double x0, double y0, BoundaryCondition bc){
				this->r = r;
				this->h = h;
				center = Coord(x0,y0);
				setBC(bc);
				make();
			}
			Circle(double r, double h, double x0, double y0){
				this->r = r;
				this->h = h;
				center = Coord(x0,y0);
				make();
			}
			/* Many of these methods are pretty self-explanatory. The
			 * getPoints(), getTypes(), and getSpacings() methods may
			 * want to be avoided, since they involve copying a whole
			 * array. Rather, use the getPoint(i,j), getType(i,j), and
			 * getSpacing(i,j) methods to only see what the Coord,
			 * CoordType, or Spacing at a particular grid point is, when
			 * needed elsewhere in the code.
			 */
			virtual double getR(){
				return r;
			}
			virtual double getH(){
				return h;
			}
			virtual Coord getCenter(){
				return center;
			}
			virtual CoordArray getPoints(){
				return points.copy();
			}
			virtual CoordTypeArray getTypes(){
				return types.copy();
			}
			virtual SpacingArray getSpacings(){
				return spacings.copy();
			}
			/* Creates a GridScalar object of the appropriate
			 * size and indexing to represent a grid function
			 * on this particular grid, and initializes it to
			 * 0 everywhere.
			 */
			virtual GridScalar makeScalar(){
				GridScalar u(xRange,yRange);
				u = 0;
				return u;
			}
			virtual CoordType getType(int i, int j){
				return types(i,j);
			}
			virtual Spacing getSpacing(int i, int j){
				return spacings(i,j);
			}
			virtual Coord getCoord(int i, int j){
				return points(i,j);
			}
			virtual double getX(int i, int j){
				return points(i,j)(0);
			}
			virtual double getY(int i, int j){
				return points(i,j)(1);
			}
			virtual GridScalar getXes(){
				return points.extractComponent(double(),0,2);
			}
			virtual GridScalar getYes(){
				return points.extractComponent(double(),1,2);
			}
			virtual Coord nearestBoundary(int i, int j, Direction d){
				double x = getX(i,j), y = getY(i,j);
				double xB, yB;
				if(d == NORTH){
					xB = x;
					yB = r + center(1);
				}
				else if(d == SOUTH){
					xB = x;
					yB = -r + center(1);
				}
				else if(d == EAST){
					yB = y;
					xB = r + center(0);
				}
				else if(d == WEST){
					yB = y;
					xB = -r + center(0);
				}
				return Coord(xB,yB);
			}
			virtual double evaluateOnNearestBoundary(int i, int j, Direction d, double (*boundary)(double,double)){
				Coord b = nearestBoundary(i,j,d);
				return (*boundary)(b(0),b(1));
			}
			virtual int interiorPoints(){
				return N;
			}
			virtual void setBC(BoundaryCondition bc){
				this->bc = bc;
			}
			virtual BoundaryCondition getBC(){
				return bc;
			}
			virtual double boundaryValue(int i, int j, Direction d){
				Coord b = nearestBoundary(i,j,d);
				return bc(b);
			}
			virtual TinyVector<double,8> extrapolateNeighbors(GridScalar u, int i, int j){
				TinyVector<double,8> out;
				if(getType(i+1,j) == EXTERIOR){
					out(EAST) = 0;
				}
				else{
					out(EAST) = u(i+1,j);
				}
				if(getType(i-1,j) == EXTERIOR){
					out(WEST) = 0;
				}
				else{
					out(WEST) = u(i-1,j);
				}
				if(getType(i,j+1) == EXTERIOR){
					out(NORTH) = 0;
				}
				else{
					out(NORTH) = u(i,j+1);
				}
				if(getType(i,j-1) == EXTERIOR){
					out(SOUTH) = 0;
				}
				else{
					out(SOUTH) = u(i,j-1);
				}
				if(getType(i+1,j+1) == EXTERIOR){
					out(NORTHEAST) = 0;
				}
				else{
					out(NORTHEAST) = u(i+1,j+1);
				}
				if(getType(i+1,j-1) == EXTERIOR){
					out(SOUTHEAST) = 0;
				}
				else{
					out(SOUTHEAST) = u(i+1,j-1);
				}
				if(getType(i-1,j+1) == EXTERIOR){
					out(NORTHWEST) = 0;
				}
				else{
					out(NORTHWEST) = u(i-1,j+1);
				}
				if(getType(i-1,j-1) == EXTERIOR){
					out(SOUTHWEST) = 0;
				}
				else{
					out(SOUTHWEST) = u(i-1,j-1);
				}
				return out;
			}
			virtual TinyVector<double,8> extrapolateNeighborsOld(GridScalar u, int i, int j){
				Spacing sp = getSpacing(i,j);
				NeighborScalar nbr = getNeighbors(u,i,j);

				double hDiag = sqrt(2)*h;
				double xC = getX(i,j), yC = getY(i,j);
				double x0 = center(0), y0 = center(1);
				double hN = sp(NORTH), hS = sp(SOUTH), hE = sp(EAST), hW = sp(WEST), hNE, hSE, hNW, hSW;
				double uC = u(i,j), uN, uS, uE, uW, uNE, uSE, uNW, uSW;
				double bN = nbr(NORTH), bS = nbr(SOUTH), bE = nbr(EAST), bW = nbr(WEST), bNE, bSE, bNW, bSW;

				/*cout << "(x0, y0) = (" << x0 << ", " << y0 << ")." << endl;
				cout << "(xC, yC) = (" << xC << ", " << yC << ")." << endl;
				cout << "r = " << r << endl;*/

				if(getType(i,j+1) == EXTERIOR){
					uN = (h/hN)*(bN-uC) + uC;
				}
				else{
					uN = u(i,j+1);
				}
				if(getType(i,j-1) == EXTERIOR){
					uS = (h/hS)*(bS-uC) + uC;
				}
				else{
					uS = u(i,j-1);
				}
				if(getType(i+1,j) == EXTERIOR){
					uE = (h/hE)*(bE-uC) + uC;
				}
				else{
					uE = u(i+1,j);
				}
				if(getType(i-1,j) == EXTERIOR){
					uW = (h/hW)*(bW-uC) + uC;
				}
				else{
					uW = u(i-1,j);
				}
				if(getType(i+1,j+1) == EXTERIOR){
					/* The boundary points xNE, yNE will lie on the line
					 * y - yC = x - xC.
					 */
					double xNE, yNE;

					/* Either we hit the boundary at a vertical edge, in which
					 * case xNE = x0 + r, or we hit at a horizontal edge, in
					 * which case yNE = y0 + r.
					 *
					 * In the first case, we get yNE = xNE - xC + yC,
					 * = x0 + r - xC + yC, and in the second case we get that
					 * xNE = yNE - yC + xC = y0 + r - yC + xC.
					 *
					 * We just need to find out which edge gets hit first. So
					 * we need to find which possible value of xNE is closer to
					 * xC.
					 */

					xNE = (abs(x0 + r - xC) < abs(y0 + r - yC)) ? (x0 + r) : (y0 + r - yC + xC);
					yNE = xNE - xC + yC;

					hNE = sqrt(pow2(xNE - xC) + pow2(yNE - yC));

					Coord NE(xNE,yNE);
					bNE = bc(NE);

					uNE = (hDiag/hNE)*(bNE-uC) + uC;
					//cout << "(xNE, yNE) = (" << xNE << ", " << yNE << ")." << endl;
				}
				else{
					uNE = u(i+1,j+1);
				}
				if(getType(i+1,j-1) == EXTERIOR){
					/* This is similar to the above, except we are on the line
					 * y - yC = -(x - xC) = xC - x.
					 */
					double xSE, ySE;

					/* So either we get that xSE = x0 + r or xSE = -(y0 - r) + yC + xC.
					 */

					xSE = (abs(x0 + r - xC) < abs(-(y0 - r) + yC)) ? (x0 + r) : (-(y0 - r) + yC + xC);
					ySE = -(xSE - xC) + yC;

					hSE = sqrt(pow2(xSE - xC) + pow2(ySE - yC));

					Coord SE(xSE,ySE);
					bSE = bc(SE);

					uSE = (hDiag/hSE)*(bSE-uC) + uC;
					//cout << "(xSE, ySE) = (" << xSE << ", " << ySE << ")." << endl;
				}
				else{
					uSE = u(i+1,j-1);
				}
				if(getType(i-1,j+1) == EXTERIOR){
					/* This is similar to the above. The line we are on is
					 * y - yC = -(x - xC) = xC - x.
					 */
					double xNW, yNW;

					/* So either we get that xNW = x0 - r or xNW = -(y0 + r) + yC + xC.
					 */

					xNW = (abs(x0 - r - xC) < abs(-(y0 + r) + yC)) ? (x0 - r) : (-(y0 + r) + yC + xC);
					yNW = -(xNW - xC) + yC;

					hNW = sqrt(pow2(xNW - xC) + pow2(yNW - yC));

					Coord NW(xNW,yNW);
					bNW = bc(NW);

					uNW = (hDiag/hNW)*(bNW-uC) + uC;
					//cout << "(xNW, yNW) = (" << xNW << ", " << yNW << ")." << endl;
				}
				else{
					uNW = u(i-1,j+1);
				}
				if(getType(i-1,j-1) == EXTERIOR){
					/* This is similar to the above. The line we are on is
					 * y - yC = x - xC.
					 */
					double xSW, ySW;

					/* So either we get that xSW = x0 - r or xSW = y0 - r - yC + xC.
					 */

					xSW = (abs(x0 - r - xC) < abs(y0 - r - yC)) ? x0 - r : y0 - r - yC + xC;
					ySW = xSW - xC + yC;

					hSW = sqrt(pow2(xSW - xC) + pow2(ySW - yC));

					Coord SW(xSW,ySW);
					bSW = bc(SW);

					uSW = (hDiag/hSW)*(bSW-uC) + uC;
					//cout << "(xSW, ySW) = (" << xSW << ", " << ySW << ")." << endl;
				}
				else{
					uSW = u(i-1,j-1);
				}
				TinyVector<double,8> out;
				out(NORTH) = uN;
				out(SOUTH) = uS;
				out(EAST) = uE;
				out(WEST) = uW;
				out(NORTHEAST) = uNE;
				out(NORTHWEST) = uNW;
				out(SOUTHEAST) = uSE;
				out(SOUTHWEST) = uSW;
				return out;
			}
			virtual NeighborScalar getNeighbors(GridScalar u, int i, int j){
				double uN, uS, uE, uW;
				if(getType(i,j+1) == EXTERIOR){
					uN = boundaryValue(i,j,NORTH);
				}
				else{
					uN = u(i,j+1);
				}
				if(getType(i,j-1) == EXTERIOR){
					uS = boundaryValue(i,j,SOUTH);
				}
				else{
					uS = u(i,j-1);
				}
				if(getType(i+1,j) == EXTERIOR){
					uE = boundaryValue(i,j,EAST);
				}
				else{
					uE = u(i+1,j);
				}
				if(getType(i-1,j) == EXTERIOR){
					uW = boundaryValue(i,j,WEST);
				}
				else{
					uW = u(i-1,j);
				}
				NeighborScalar out;
				out(NORTH) = uN;
				out(SOUTH) = uS;
				out(EAST) = uE;
				out(WEST) = uW;
				return out;
			}
			virtual int getM(){
				return M;
			}
		};

		class Cell{
			CellType type;
		public:
			virtual ~Cell(){}
			virtual CellType getType(){
				return type;
			}
		};

		class RegularCell : public Cell{
		public:
			virtual ~RegularCell(){}
			virtual CellType getType(){
				return REGULARCELL;
			}
		};

		class IrregularCell : public Cell{
			TinyVector<double,2> N; // Inner unit normal to boundary face
			Coord I1, I2; // Intersections of cell faces with the boundary
			double Gamma; // Volume fraction
			Coord mN, mS, mE, mW, mB; // Midpoints of faces
			double A; // Area of the boundary front
			double aN, aS, aE, aW; // Apertures (area fractions) of cell faces
		public:
			virtual ~IrregularCell(){}
			virtual CellType getType(){
				return IRREGULARCELL;
			}
		};

		/* Abstract base class for operators on single grid functions.
		 * The function apply() does exactly what it sounds like.
		 */
		class GridOperator{
		public:
			virtual ~GridOperator(){}
			virtual GridScalar apply(GridScalar u, Circle circ) = 0;
			virtual GridScalar operator() (GridScalar u, Circle circ){
				return apply(u, circ);
			}
		};
	}
}

#endif /* GEOMETRY_H_ */
