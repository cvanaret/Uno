#ifndef BQPDSOLVER_H
#define BQPDSOLVER_H

#include <iostream>
#include <vector>
#include <map>
#include "QPSolver.hpp"
#include "LPSolver.hpp"

enum BQPDMode {
    COLD_START = 0,
    ACTIVE_SET_EQUALITIES = 1,
    USER_DEFINED = 2,
    ACTIVE_SET_PREVIOUS_CALL = 3
    // goes up to 6
};

/*! \class BQPDSolver
 * \brief Interface for BQPD
 *
 *  Interface to the QP/LP solver BQPD
 */
class BQPDSolver : public QPSolver {
public:
    /*!
     *  Constructor
     * 
     * \param n: number of variables
     * \param m: number of constraints
     */
    BQPDSolver(int number_variables, int number_constraints, int max_number_nonzeros);

    /*!
     *  Solve a QP
     * 
     * \param d0: initial point
     */
    SubproblemSolution solve_QP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, CSCMatrix& hessian, std::vector<double>& x);

    /*!
     *  Solve an LP
     * 
     * \param d0: initial point
     */
    SubproblemSolution solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x);


private:
    int n_, m_;
    int maximum_number_nonzeros;
    std::vector<double> lb, ub; // lower and upper bounds of variables and constraints
    int use_fortran;
    
    std::vector<double> jacobian;
    std::vector<int> jacobian_sparsity;
    int kmax_, mlp_, mxwk0_, mxiwk0_;
    std::vector<int> info_;
    std::vector<double> alp_;
    std::vector<int> lp_, ls_;
    std::vector<double> w_, gradient_solution_, residuals_, e_;
    int nhr_, nhi_, mxws_, mxlws_;
    std::vector<double> ws_;
    std::vector<int> lws_;
    int k_;
    BQPDMode mode_;
    int iprint_, nout_;
    double fmin_, f_solution_;
    int peq_solution_, ifail_;

    /*!
     *  Create a SubproblemSolution from BQPD's solution
     * 
     * \param d: optimal solution
     */
    SubproblemSolution generate_solution(std::vector<double>& x);
    Status int_to_status(int ifail);
    SubproblemSolution solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x, int kmax);
};

#endif // BQPDSOLVER_H
