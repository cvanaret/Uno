#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include "Solver.hpp"
#include <vector>

/*! \class BQPDSolver
* \brief Interface for BQPD
*
*  Interface to the QP/LP solver BQPD
*/
class MA57Solver: public Solver {
	public:
		MA57Solver(int n, int m, int hessian_nnz);
	
		/*!
         *  Solve a QP using BQPD
         */
		Step solve(QP& qp, std::vector<double>& d0);
		
		Step solve(LP& qp, std::vector<double>& d0);
		
	private:
		
};

#endif // MA57SOLVER_H
