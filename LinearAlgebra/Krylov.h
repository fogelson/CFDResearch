/*
 * Krylov.h
 *
 *  Created on: Feb 21, 2011
 *      Author: bfogelson
 */

#ifndef KRYLOV_H_
#define KRYLOV_H_

#include "../Multigrid/GridOperators.h"
#include "../Geometry.h"
#include "LinearAlgebra.h"
#include "Norms.h"

using namespace CFD;
using namespace Geometry;
using namespace Multigrid;

namespace CFD{
	namespace LinearAlgebra{
		class CG{
		public:
			virtual ~CG(){}
			virtual GridScalar solve(GridOperator *A, GridScalar u0, GridScalar b, Circle circ, Norm *norm, double tol, int maxIts){
				GridScalar r = circ.makeScalar();
				GridScalar p = circ.makeScalar();
				GridScalar w = circ.makeScalar();
				GridScalar u = circ.makeScalar();
				double alpha, beta;
				GridScalar rNew = circ.makeScalar();
				GridScalar pNew = circ.makeScalar();
				GridScalar uNew = circ.makeScalar();

				Circle zCirc(circ.getR(),circ.getH());

				u = u0;
				r = b - A->apply(u,circ);
				p = r;
				for(int k = 0; k < maxIts; k++){
					w = A->apply(p,zCirc);
					double rDotR = sum(where(circ.getTypes() != EXTERIOR, r*r, 0));
					double pDotW = sum(where(circ.getTypes() != EXTERIOR, p*w, 0));
					alpha = rDotR/pDotW;
					uNew = u + alpha*p;
					rNew = r - alpha*w;
					if(norm->apply(rNew,circ) < tol){
						cout << "Converged to within desired tolerance in " << k << " iterations. Residual has norm "
								<< norm->apply(rNew,circ) << endl;
						return uNew;
					}
					double rNewDotRNew = sum(where(circ.getTypes() != EXTERIOR, rNew*rNew, 0));
					beta = rNewDotRNew/rDotR;
					pNew = rNew + beta*p;
					r = rNew;
					p = pNew;
					u = uNew;
				}
				cout << "Maximum iterations reached. Residual has norm " << norm->apply(r,circ) << endl;
				return uNew;
			}
		};
		class PCG{
		public:
			virtual ~PCG(){}
			virtual GridScalar solve(Solver *M, GridOperator *A, GridScalar u0, GridScalar b, Circle circ, Norm *norm, double tol, int maxIts){
				GridScalar r = circ.makeScalar();
				GridScalar p = circ.makeScalar();
				GridScalar w = circ.makeScalar();
				GridScalar u = circ.makeScalar();
				GridScalar z = circ.makeScalar();
				double alpha, beta;
				GridScalar rNew = circ.makeScalar();
				GridScalar pNew = circ.makeScalar();
				GridScalar uNew = circ.makeScalar();
				GridScalar zNew = circ.makeScalar();

				Circle zCirc(circ.getR(),circ.getH());

				u = u0;
				r = b - A->apply(u,circ);
				z = M->solve(r);
				p = z;
				for(int k = 0; k < maxIts; k++){
					w = A->apply(p,zCirc);
					double zDotR = sum(where(circ.getTypes() != EXTERIOR, z*r, 0));
					double pDotW = sum(where(circ.getTypes() != EXTERIOR, p*w, 0));
					alpha = zDotR/pDotW;
					uNew = u + alpha*p;
					rNew = r - alpha*w;
					if(norm->apply(rNew,circ) < tol){
						cout << "Converged to within desired tolerance in " << k << " iterations. Residual has norm "
								<< norm->apply(rNew,circ) << endl;
						return uNew;
					}
					zNew = M->solve(rNew);
					double zNewDotRNew = sum(where(circ.getTypes() != EXTERIOR, zNew*rNew, 0));
					beta = zNewDotRNew/zDotR;
					pNew = zNew + beta*p;
					r = rNew;
					p = pNew;
					u = uNew;
					z = zNew;
				}
				cout << "Maximum iterations reached. Residual has norm " << norm->apply(r,circ) << endl;
				return uNew;
			}
		};
	}
}

#endif /* KRYLOV_H_ */
