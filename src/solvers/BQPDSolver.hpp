#ifndef BQPDSOLVER_H
#define BQPDSOLVER_H

#include <iostream>
#include <vector>
#include <map>
#include "QPSolver.hpp"
#include "LPSolver.hpp"

/*! \class BQPDSolver
* \brief Interface for BQPD
*
*  Interface to the QP/LP solver BQPD
*/
class BQPDSolver: public QPSolver, public LPSolver {
	public:
		/*!
         *  Constructor
         * 
         * \param n: number of variables
         * \param m: number of constraints
         * \param hessian_nnz: number of nonzero terms in the Hessian
         */
		BQPDSolver(std::vector<int>& hessian_column_start, std::vector<int>& hessian_row_number);
	
		void allocate(int n, int m);
	
		/*!
         *  Solve a QP
         * 
         * \param qp: quadratic program
         * \param d0: initial point
         */
		LocalSolution solve(QP& qp, std::vector<double>& x0);
		
		/*!
         *  Solve an LP
         * 
         * \param lp: linear program
         * \param d0: initial point
         */
		LocalSolution solve(LP& qp, std::vector<double>& x0);

		
	private:
		int kmax_, mlp_, mxwk0_, mxiwk0_;
		std::vector<int> info_;
		std::vector<double> alp_;
		std::vector<int> lp_, ls_;
		std::vector<double> w_, gradient_solution_, residuals_, e_;
		int hessian_nnz_, nhr_, nhi_, mxws_, mxlws_;
		std::vector<double> ws_;
		std::vector<int> lws_;
		int n_, m_, k_, mode_, iprint_, nout_;
		double fmin_, f_solution_;
		int peq_solution_, ifail_;
		std::vector<int> hessian_column_start, hessian_row_number;
		
		/*!
         *  Create a LocalSolution from BQPD's solution
         * 
         * \param d: optimal solution
         */
		LocalSolution generate_solution(std::vector<double>& x);
		
		void build_jacobian(std::vector<double>& full_jacobian, std::vector<int>& full_jacobian_sparsity, std::map<int,double>& jacobian);
};

#endif // BQPDSOLVER_H
