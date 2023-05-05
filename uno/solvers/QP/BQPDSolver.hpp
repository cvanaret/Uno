// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDSOLVER_H
#define UNO_BQPDSOLVER_H

#include <vector>
#include "QPSolver.hpp"
#include "solvers/LP/LPSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "tools/Options.hpp"

// see bqpd.f
enum class BQPDStatus {
   OPTIMAL = 0,
   UNBOUNDED_PROBLEM = 1,
   BOUND_INCONSISTENCY = 2,
   INFEASIBLE = 3,
   INCORRECT_PARAMETER = 4,
   LP_INSUFFICIENT_SPACE = 5,
   HESSIAN_INSUFFICIENT_SPACE = 6,
   SPARSE_INSUFFICIENT_SPACE = 7,
   MAX_RESTARTS_REACHED = 8,
   UNDEFINED = 9
};

enum BQPDMode {
   COLD_START = 0, // cold start
   ACTIVE_SET_EQUALITIES = 1, // hot start
   USER_DEFINED = 2,
   UNCHANGED_ACTIVE_SET = 3,
   UNCHANGED_ACTIVE_SET_AND_JACOBIAN = 4,
   UNCHANGED_ACTIVE_SET_AND_REDUCED_HESSIAN = 5, // warm start
   UNCHANGED_ACTIVE_SET_AND_JACOBIAN_AND_REDUCED_HESSIAN = 6,
};

class BQPDSolver : public QPSolver {
public:
   BQPDSolver(size_t max_number_variables, size_t number_constraints, size_t number_hessian_nonzeros, bool quadratic_programming, const Options& options);

   Direction solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) override;

   Direction solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) override;

private:
   size_t number_hessian_nonzeros;
   std::vector<double> lb, ub; // lower and upper bounds of variables and constraints

   std::vector<double> jacobian;
   std::vector<int> jacobian_sparsity;
   int kmax, mlp{1000};
   size_t mxwk0{2000000}, mxiwk0{500000};
   std::array<int, 100> info{};
   std::vector<double> alp;
   std::vector<int> lp, active_set;
   std::vector<double> w, gradient_solution, residuals, e;
   size_t size_hessian_sparsity;
   size_t size_hessian_workspace;
   size_t size_hessian_sparsity_workspace;
   std::vector<double> hessian_values{};
   std::vector<int> hessian_sparsity{};
   int k{0};
   int iprint{0}, nout{6};
   double fmin{-1e20};
   int peq_solution{}, ifail{};
   const int fortran_shift{1};
   const bool print_subproblem;

   Direction solve_subproblem(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information);
   void analyze_constraints(size_t number_variables, size_t number_constraints, Direction& direction);
   void save_lagrangian_hessian_to_local_format(const SymmetricMatrix<double>& hessian);
   void save_gradients_to_local_format(size_t number_constraints, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian);
   static BQPDMode determine_mode(const WarmstartInformation& warmstart_information);
   static BQPDStatus bqpd_status_from_int(int ifail);
   static SubproblemStatus status_from_bqpd_status(BQPDStatus bqpd_status);
};

#endif // UNO_BQPDSOLVER_H
