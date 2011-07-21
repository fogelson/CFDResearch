/*
 * LinearAlgebra.h
 *
 *  Created on: Feb 21, 2011
 *      Author: bfogelson
 */

#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

#include "../Geometry.h"

using namespace CFD;
using namespace Geometry;

namespace CFD{
	namespace LinearAlgebra{
		class LinearOperator{
		public:
			virtual ~LinearOperator(){}
			virtual GridScalar apply(GridScalar in) = 0;
		};
		class Solver{
		public:
			virtual ~Solver(){}
			virtual GridScalar solve(GridScalar in) = 0;
		};

		class IdentitySolver : public Solver{
		public:
			virtual ~IdentitySolver(){}
			virtual GridScalar solve(GridScalar in){
				return in;
			}
		};
	}
}

#endif /* LINEARALGEBRA_H_ */
