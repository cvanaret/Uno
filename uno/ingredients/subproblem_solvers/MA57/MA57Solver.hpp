// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <vector>
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/COOEvaluationSpace.hpp"

namespace uno {
   // forward declarations
   class Statistics;
   class Subproblem;

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

   class MA57Solver : public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      MA57Solver();
      ~MA57Solver() override = default;

      void initialize_hessian(const Subproblem& subproblem) override;
      void initialize_augmented_system(const Subproblem& subproblem) override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(const double* matrix_values) override;
      void solve_indefinite_system(const Vector<double>& matrix_values, const Vector<double>& rhs, Vector<double>& result) override;
      void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override;

   private:
      MA57Workspace workspace{};
      COOEvaluationSpace evaluation_space{};

      bool analysis_performed{false};
      bool factorization_performed{false};

      bool use_iterative_refinement{false};
   };
} // namespace

#endif // UNO_MA57SOLVER_H