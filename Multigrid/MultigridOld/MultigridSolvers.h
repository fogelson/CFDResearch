/*
 * MultigridSolvers.h
 *
 *  Created on: Feb 14, 2011
 *      Author: fogelson
 */

#ifndef MULTIGRIDSOLVERS_H_
#define MULTIGRIDSOLVERS_H_

#include "../Geometry.h"
#include "Smoothers.h"
#include "LineSmoothers.h"
#include "IntergridOperators.h"
#include "GridOperators.h"


using namespace CFD;
using namespace Geometry;

namespace CFD{
	namespace Multigrid{


		double maxNorm(GridScalar u, Circle circ){
			return max(where(circ.getTypes() != EXTERIOR, abs(u), 0));
		}

		/* Standard Multigrid Solver for a circular domain.
		 * Stores pointers to GridOperator, Smoother, Interpolator
		 * and Restrictor class objects, and uses those to solve
		 * a PDE on the circle.
		 *
		 * An example constructor would involve first declaring the
		 * operator classes, and passing references to them into the
		 * constructor for Multigrid:
		 *
		 * 		LaplaceOperator Lu;		// GridOperator for difference operator
		 * 		GSRBPoisson gsrb;		// Smoother
		 * 		BilinearInterpolator bi	// Interpolator
		 * 		HalfWeighter hw			// Restrictor
		 *
		 * 		// Declare the solver
		 * 		MultigridSolver mg(&Lu, &gsrb, &bi, &hw);
		 *
		 * 		// Now we can use it to solve problems, by calling:
		 * 		GridScalar u = mg.solve(u0, f, circ, v1, v2, its);
		 *
		 * See comments for the various methods for a discussion of them.
		 */
		class MultigridSolver{
			GridOperator *differenceOperator;
			Smoother *smoother;
			Interpolator *interpolator;
			Restrictor *restrictor;
			/* Runs one Multigrid vCycle with initial value u0, right hand side
			 * f, domain fine, v1 presmooths, and v2 postsmooths.
			 */
			GridScalar vCycle(GridScalar u0, GridScalar f, Circle fine, int v1, int v2){
				//bool printDebug = false;

				double hF = fine.getH();
				double hC = 2.0*hF;

				/* Declare a circle over the fine domain that has zero boundary
				 * condition. This is necessary for transfer operations to
				 * work correctly, since we are interpolating on the error equation
				 * and our error at the boundary should always be zero, since
				 * the boundary condition is prescribed.
				 */
				Circle fineZero(fine.getR(),fine.getH(),fine.getCenter());

				/* Declare the coarse domain. Note that we use a constructor
				 * for the Circle() that doesn't declare the boundary function. This is
				 * because we will be solving the error equation, so it should be zero
				 * on the boundary.
				 */
				Circle coarse(fine.getR(),hC,fine.getCenter());

				/* Initialize GridScalars to store different intermediate values
				 * needed for the vCycle. We could probably make the code more
				 * memory (and possibly speed) efficient by using fewer objects,
				 * but for the time being this is fine.
				 */
				GridScalar u = fine.makeScalar();
				GridScalar Lu = fine.makeScalar();
				GridScalar rF = fine.makeScalar();
				GridScalar rC = coarse.makeScalar();
				GridScalar eC = coarse.makeScalar();
				GridScalar eF = fine.makeScalar();
				GridScalar uCorrected = fine.makeScalar();

				// Presmooth.
				u = smoother->smooth(u0,f,fine,v1);

				// Apply difference operator.
				Lu = differenceOperator->apply(u,fine);

				// Calculate the residual on the fine grid.
				rF = f - Lu;

				// Restrict the residual to the coarse grid.
				rC = restrictor->doRestrict(rF,fineZero,coarse);

				/* If the coarse grid is very coarse, solve directly.
				 * Here we do a "direct" solve by doing just a few
				 * iterations with our smoother.
				 *
				 * If the coarse grid is not very coarse, we recursively
				 * call the vCycle function, but now solving the error
				 * equation on the coarse grid, with initial guess eC = 0.
				 */

				if(coarse.interiorPoints() <= 5){
					// Get the error on the coarse grid.
					//PoissonGSRB gsrb;
					//eC = gsrb.smooth(eC,rC,coarse,300);
					eC = smoother->smooth(eC,rC,coarse,10);
				}
				else{
					// Get the error on the coarse grid.
					eC = vCycle(eC, rC, coarse, v1, v2);
				}

				//PoissonGSRB gsrb;
				//eC = smoother->smooth(eC,rC,coarse,8000);

				// Interpolate the error from the coarse to fine grids.
				eF = interpolator->doInterpolate(eC, fineZero, coarse);

				// Use the error to correct the solution.
				uCorrected = u + eF;

				// Postsmooth.
				u = smoother->smooth(uCorrected,f,fine,v2);
				return u;
			}

		public:
			/* Constructor for MultigridSolver. Requires pointers to instances of the
			 * GridOperator, Smoother, Interpolator, and Restrictor base classes.
			 */
			MultigridSolver(GridOperator *differenceOperator, Smoother *smoother, Interpolator *interpolator, Restrictor *restrictor){
				this->differenceOperator = differenceOperator;
				this->smoother = smoother;
				this->interpolator = interpolator;
				this->restrictor = restrictor;
			}
			/* Solve the PDE given by the differenceOperator and smoother pointers
			 * with initial condition u0, right hand side f, domain circ, v1 presmooths,
			 * and v2 postsmooths by running the number of vCycles given by its.
			 */
			GridScalar solve(GridScalar u0, GridScalar f, Circle circ, int v1, int v2, int its){
				GridScalar u = circ.makeScalar();
				u = u0;
				MaxNorm mn;
				for(int n = 0; n < its; n++){
					u = vCycle(u, f, circ, v1, v2);
				}
				return u;
			}
		};




	}
}

#endif /* MULTIGRIDSOLVERS_H_ */
