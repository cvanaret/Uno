// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <vector>
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Statistics;

   struct MA57Workspace {
      int n{};
      int nnz{};
      int lfact{};
      int lifact{};

      std::vector<double> fact{0}; // do not initialize, resize at every iteration
      std::vector<int> ifact{0}; // do not initialize, resize at every iteration
      int lkeep{};
      std::vector<int> keep{};
      std::vector<int> iwork{};
      int lwork{};
      std::vector<double> work{};

      // for ma57id_ (default values of controlling parameters)
      std::array<double, 5> cntl{};
      std::array<int, 20> icntl{};
      std::array<double, 20> rinfo{};
      std::array<int, 40> info{};

      const int nrhs{1}; // number of right hand side being solved
      const int job{1};
      std::vector<double> residuals;

      MA57Workspace() = default;
   };

   class MA57Solver : public DirectSymmetricIndefiniteLinearSolver<size_t, double> {
   public:
      MA57Solver();
      ~MA57Solver() override = default;

      void initialize(const Subproblem& subproblem) override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(const Vector<double>& matrix_values) override;
      void solve_indefinite_system(const Vector<double>& matrix_values, const Vector<double>& rhs, Vector<double>& result) override;
      void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

   private:
      MA57Workspace workspace{};

      // evaluations
      Vector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

      // augmented system
      std::vector<int> row_indices{};
      std::vector<int> column_indices{};
      Vector<double> matrix_values{};
      Vector<double> rhs{};
      Vector<double> solution{};

      static constexpr size_t fortran_shift{1};

      bool use_iterative_refinement{false};
   };
} // namespace

#endif // UNO_MA57SOLVER_H