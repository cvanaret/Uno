// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <cstddef>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "../PivotingToleranceSchedule.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // settings taken from IPOPT
   struct MA57Settings {
      static constexpr double pivoting_threshold = 1e-8;
      static constexpr double max_pivoting_threshold = 1e-4; // https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_ma57_pivtolmax
      static constexpr double pivoting_threshold_increase_exponent = 0.75;
      static constexpr double allocation_safety_factor = 1.05;
      // the workspace grows by at least 20% per enlargement (the factorization continues where it stopped, ICNTL(8) = 1)
      static constexpr size_t max_number_workspace_enlargements = 50;
   };

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

      const int job{1};

      MA57Workspace() = default;
   };

   class MA57Solver : public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      MA57Solver(bool use_scaling = false); // no scaling by default
      ~MA57Solver() override = default;

      void initialize_memory() override;

      void do_symbolic_analysis() override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;
      [[nodiscard]] bool increase_pivoting_tolerance() override;
      [[nodiscard]] std::optional<size_t> number_delayed_pivots() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system() override;
      // whether MA57 scales the matrix itself (MC64)
      [[nodiscard]] bool uses_internal_scaling() const { return this->workspace.icntl[14] != 0; }

   private:
      void factorize(bool is_matrix_positive_definite) override;
      void solve(const double* rhs, double* solution, size_t number_of_rhs) override;

      MA57Workspace workspace{};
      COOLinearSystem linear_system{Indexing::Fortran_indexing};
      PivotingToleranceSchedule pivoting_tolerance_schedule{MA57Settings::max_pivoting_threshold, MA57Settings::pivoting_threshold_increase_exponent};

      bool analysis_performed{false};
      bool factorization_performed{false};

      [[nodiscard]] int& MA57_ICNTL(size_t index);
      [[nodiscard]] double& MA57_CNTL(size_t index);
      [[nodiscard]] int MA57_INFO(size_t index) const;
   };
} // namespace

#endif // UNO_MA57SOLVER_H