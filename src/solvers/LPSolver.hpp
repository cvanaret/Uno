#ifndef LPSOLVER_H
#define LPSOLVER_H

#include <vector>
#include "LocalSolution.hpp"
#include "LP.hpp"

/*! \class LPSolver
* \brief LP solver
*
*/
class LPSolver {
	public:
		/*!
         *  Solve an LP
         * 
         * \param lp: linear program
         * \param d0: initial point
         */
		virtual LocalSolution solve(LP& lp, std::vector<double>& x0) = 0;
		
		virtual void allocate(int n, int m) = 0;
};

#endif // LPSOLVER_H
