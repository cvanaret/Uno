#ifndef QPSOLVER_H
#define QPSOLVER_H

#include <vector>
#include "LPSolver.hpp"
#include "SubproblemSolution.hpp"
#include "Matrix.hpp"

/*! \class QPSolver
* \brief QP solver
*
*/
class QPSolver: public LPSolver {
	public:
		/*!
         *  Solve a QP
         * 
         * \param qp: quadratic program
         * \param d0: initial point
         */
        virtual ~QPSolver() {};
		virtual SubproblemSolution solve_QP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int,double>& linear_objective, std::vector<std::map<int,double> >& constraints_jacobian, CSCMatrix& hessian, std::vector<double>& x0) = 0;
		virtual SubproblemSolution solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int,double>& linear_objective, std::vector<std::map<int,double> >& constraints_jacobian, std::vector<double>& x0) = 0;
		
		//virtual void allocate(int n, int m) = 0;
};

#endif // QPSOLVER_H
