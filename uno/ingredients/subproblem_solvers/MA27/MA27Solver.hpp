// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   struct MA27Workspace {
      int n{};                         // dimension of current factorisation (maximal value here is <= max_dimension)
      int nnz{};                     // number of nonzeros of current factorisation
      std::array<int, 30> icntl{};      // integer array of length 30; integer control values
      std::array<double, 5> cntl{};     // double array of length 5; double control values

      std::vector<int> iw{};           // integer workspace of length liw
      std::vector<int> ikeep{};        // integer array of 3*n; pivot sequence
      std::vector<int> iw1{};          // integer workspace array of length n
      int nsteps{};                    // integer, to be set by ma27
      int iflag{};                     // integer; 0 if pivot order chosen automatically; 1 if the pivot order set by ikeep
      std::array<int, 20> info{};       // integer array of length 20
      double ops{};                    // double, operations count

      std::vector<double> factor{};    // data array of length la;
      int maxfrt{};                    // integer, to be set by ma27
      std::vector<double> w{};         // double workspace
      const size_t number_factorization_attempts{5};
   };

   class MA27Solver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      MA27Solver();
      ~MA27Solver() override = default;

      void initialize_memory() override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* result) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;

   private:
      MA27Workspace workspace{};
      COOLinearSystem linear_system{Indexing::Fortran_indexing};

      bool analysis_performed{false};
      bool factorization_performed{false};

      // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
      void check_factorization_status();
   };
} // namespace

#endif // UNO_MA27SOLVER_H