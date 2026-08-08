// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // settings taken from IPOPT
   struct MA57Settings {
      static constexpr int mc64_scaling = 0;
      static constexpr double pivoting_threshold = 1e-8;
      static constexpr double allocation_safety_factor = 1.05;
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
      MA57Solver();
      ~MA57Solver() override = default;

      void initialize_memory() override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* solution) override;
      void solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system();

   private:
      MA57Workspace workspace{};
      COOLinearSystem linear_system{Indexing::Fortran_indexing};

      bool analysis_performed{false};
      bool factorization_performed{false};

      [[nodiscard]] int& ICNTL(size_t index);
      [[nodiscard]] double& CNTL(size_t index);
      [[nodiscard]] int INFO(size_t index) const;
   };
} // namespace

#endif // UNO_MA57SOLVER_H