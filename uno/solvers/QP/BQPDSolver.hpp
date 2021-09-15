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

class BQPDSolver: public QPSolver<CSCSymmetricMatrix> {
public:
   BQPDSolver(size_t number_variables, size_t number_constraints, size_t max_number_nonzeros, bool quadratic_programming);

   Direction solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
         linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const std::vector<double>& initial_point) override;

   Direction solve_QP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
         linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const CSCSymmetricMatrix& hessian, const std::vector<double>&
               initial_point) override;

private:
   size_t number_variables, number_constraints;
   size_t maximum_number_nonzeros;
   std::vector<double> lb, ub; // lower and upper bounds of variables and constraints
   const size_t fortran_shift{1};

   std::vector<double> jacobian;
   std::vector<int> jacobian_sparsity;
   int kmax, mlp{1000};
   size_t mxwk0{2000000}, mxiwk0{500000};
   std::array<int, 100> info{};
   std::vector<double> alp;
   std::vector<int> lp, ls;
   std::vector<double> w, gradient_solution, residuals, e;
   size_t size_hessian_sparsity;
   size_t size_hessian_workspace;
   size_t size_hessian_sparsity_workspace;
   std::vector<double> hessian_values;
   std::vector<int> hessian_sparsity;
   int k{0};
   BQPDMode mode{COLD_START};
   int iprint{0}, nout{6};
   double fmin{-1e20};
   int peq_solution{}, ifail{};

   static Status int_to_status(int ifail);
   Direction solve_subproblem(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
   linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const std::vector<double>& initial_point);
   void analyze_constraints(Direction& direction);
};

#endif // BQPDSOLVER_H
