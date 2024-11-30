// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDSOLVER_H
#define UNO_BQPDSOLVER_H

#include <array>
#include <vector>
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "solvers/QPSolver.hpp"
#include "solvers/SubproblemStatus.hpp"

namespace uno {
   // forward declaration
   class Options;

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
      COLD_START = 0,
      ACTIVE_SET_EQUALITIES = 1, // cold start
      USER_DEFINED = 2, // hot start
      UNCHANGED_ACTIVE_SET = 3,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN = 4,
      UNCHANGED_ACTIVE_SET_AND_REDUCED_HESSIAN = 5,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN_AND_REDUCED_HESSIAN = 6, // warm start
   };

   enum class BQPDProblemType {LP, QP};

   class BQPDSolver : public QPSolver {
   public:
      BQPDSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, BQPDProblemType problem_type, const Options& options);

      void solve_LP(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

      void solve_QP(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

   private:
      const size_t number_hessian_nonzeros;
      Vector<double> lower_bounds{}, upper_bounds{}; // lower and upper bounds of variables and constraints

      std::vector<double> jacobian{};
      std::vector<int> jacobian_sparsity{};
      int kmax{0}, mlp{1000};
      size_t mxwk0{2000000}, mxiwk0{500000};
      std::array<int, 100> info{};
      std::vector<double> alp{};
      std::vector<int> lp{}, active_set{};
      std::vector<double> w{}, gradient_solution{}, residuals{}, e{};
      size_t size_hessian_sparsity{};
      size_t size_hessian_workspace{};
      size_t size_hessian_sparsity_workspace{};
      std::vector<double> workspace{};
      std::vector<int> workspace_sparsity{};
      int k{0};
      int iprint{0}, nout{6};
      double fmin{-1e20};
      int peq_solution{0}, ifail{0};
      const int fortran_shift{1};

      size_t number_calls{0};
      const bool print_subproblem;

      SparseVector<double> linear_objective; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

      void set_up_subproblem(const WarmstartInformation& warmstart_information, const LagrangeNewtonSubproblem& subproblem);
      void solve_subproblem(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information);
      void categorize_constraints(size_t number_variables, Direction& direction);
      void save_hessian(const LagrangeNewtonSubproblem& subproblem);
      void save_gradients_to_local_format(size_t number_constraints, const SparseVector<double>& linear_objective,
            const RectangularMatrix<double>& constraint_jacobian);
      [[nodiscard]] BQPDMode determine_mode(const WarmstartInformation& warmstart_information) const;
      static BQPDStatus bqpd_status_from_int(int ifail);
      static SubproblemStatus status_from_bqpd_status(BQPDStatus bqpd_status);
   };
} // namespace

#endif // UNO_BQPDSOLVER_H