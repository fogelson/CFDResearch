/*
 * GeometryTypes.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef GEOMETRYTYPES_H_
#define GEOMETRYTYPES_H_

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <math.h>
#include <iostream>
#include <list>

using namespace std;
using namespace blitz;

namespace CFD{
	namespace Geometry{
		class Cell;
		class Grid;
		class Circle;

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

		// Type for a cell. Might be merged with CoordType and given
		// a more general name.
		typedef int CellType;

		typedef Array<Cell,2> CellArray;
		typedef Array<CellType,2> CellTypeArray;
	}
}

#endif /* GEOMETRYTYPES_H_ */
