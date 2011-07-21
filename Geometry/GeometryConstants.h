/*
 * GeometryConstants.h
 *
 *  Created on: Apr 17, 2011
 *      Author: fogelson
 */

#ifndef GEOMETRYCONSTANTS_H_
#define GEOMETRYCONSTANTS_H_

#include "GeometryTypes.h"

namespace CFD{
	namespace Geometry{
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

		const CellType REGULARCELL = 0;
		const CellType IRREGULARCELL = 1;
		const CellType EXTERIORCELL = 2;
	}
}

#endif /* GEOMETRYCONSTANTS_H_ */
