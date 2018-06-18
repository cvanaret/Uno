#ifndef QPSOLVER_H
#define QPSOLVER_H

#include <vector>
#include "LocalSolution.hpp"
#include "QP.hpp"

/*! \class QPSolver
* \brief QP solver
*
*/
class QPSolver {
	public:
		/*!
         *  Solve a QP
         * 
         * \param qp: quadratic program
         * \param d0: initial point
         */
        virtual ~QPSolver() {};
		virtual LocalSolution solve(QP& qp, std::vector<double>& x0) = 0;
		
		virtual void allocate(int n, int m) = 0;
};

#endif // QPSOLVER_H
