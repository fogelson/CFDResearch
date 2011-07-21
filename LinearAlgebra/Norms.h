/*
 * Norms.h
 *
 *  Created on: Feb 21, 2011
 *      Author: bfogelson
 */

#ifndef NORMS_H_
#define NORMS_H_


#endif /* NORMS_H_ */

#include "../Geometry.h"

using namespace CFD;
using namespace Geometry;

namespace CFD{
	namespace LinearAlgebra{
		class Norm{
		public:
			virtual ~Norm(){}
			virtual double apply(GridScalar u, Circle circ) = 0;
			virtual double operator() (GridScalar u, Circle circ){
				return apply(u,circ);
			}
		};

		class MaxNorm : public Norm{
		public:
			virtual ~MaxNorm(){}
			virtual double apply(GridScalar u, Circle circ){
				return max(where(circ.getTypes() != EXTERIOR, abs(u), 0));
			}
		};
	}
}
