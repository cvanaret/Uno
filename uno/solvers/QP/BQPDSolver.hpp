#ifndef BQPDSOLVER_H
#define BQPDSOLVER_H

#include <iostream>
#include <vector>
#include <map>
#include "QPSolver.hpp"
#include "LPSolver.hpp"
#include "linear_algebra/CSCSymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"

enum BQPDMode {
   COLD_START = 0, ACTIVE_SET_EQUALITIES = 1, USER_DEFINED = 2, UNCHANGED_ACTIVE_SET = 3, UNCHANGED_ACTIVE_SET_AND_JACOBIAN = 4,
   UNCHANGED_ACTIVE_SET_AND_REDUCED_HESSIAN = 5, UNCHANGED_ACTIVE_SET_AND_JACOBIAN_AND_REDUCED_HESSIAN = 6,
};

/*! \class BQPDSolver
 * \brief Interface for BQPD
 *
 *  Interface to the QP/LP solver BQPD
 */
class BQPDSolver : public QPSolver<CSCSymmetricMatrix> {
public:
   BQPDSolver(size_t number_variables, size_t number_constraints, size_t max_number_nonzeros, bool quadratic_programming);

   Direction solve_QP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector& linear_objective,
         const std::vector<SparseVector>& constraints_jacobian, const CSCSymmetricMatrix& hessian, const std::vector<double>& initial_point) override;

   Direction solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector& linear_objective,
         const std::vector<SparseVector>& constraints_jacobian, const std::vector<double>& initial_point) override;

private:
   size_t n_, m_;
   size_t maximum_number_nonzeros_;
   std::vector<double> lb_, ub_; // lower and upper bounds of variables and constraints
   short use_fortran_;

   std::vector<double> jacobian_;
   std::vector<int> jacobian_sparsity_;
   int kmax_, mlp_, mxwk0_, mxiwk0_;
   std::vector<int> info_;
   std::vector<double> alp_;
   std::vector<int> lp_, ls_;
   std::vector<double> w_, gradient_solution_, residuals_, e_;
   int size_hessian_sparsity_, size_hessian_workspace_, size_hessian_sparsity_workspace_;
   std::vector<double> hessian_;
   std::vector<int> hessian_sparsity_;
   int k_;
   BQPDMode mode_;
   int iprint_, nout_;
   double fmin_, f_solution_;
   int peq_solution_, ifail_;
   std::vector<double> solution;

   /*!
    *  Create a SubproblemSolution from BQPD's solution
    *
    * \param d: optimal solution
    */
   Direction generate_direction();
   Status int_to_status_(int ifail);
   Direction solve_subproblem(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector&
   linear_objective, const std::vector<SparseVector>& constraints_jacobian, const std::vector<double>& x);
};

#endif // BQPDSOLVER_H
