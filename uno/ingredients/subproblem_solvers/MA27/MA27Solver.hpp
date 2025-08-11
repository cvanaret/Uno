// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/Vector.hpp"

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

      void evaluate_constraint_jacobian(const Subproblem& subproblem) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;

   private:
      MA27Workspace workspace{};

      // evaluations
      Vector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */

      // Jacobian
      size_t number_jacobian_nonzeros{};
      std::vector<size_t> jacobian_row_indices{};
      std::vector<size_t> jacobian_column_indices{};

      // augmented system
      size_t number_hessian_nonzeros{};
      std::vector<int> matrix_row_indices{};          // row index of input
      std::vector<int> matrix_column_indices{};          // col index of input
      Vector<double> matrix_values{};
      Vector<double> rhs{};
      Vector<double> solution{};

      static constexpr size_t fortran_shift{1};

      // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
      void check_factorization_status();
   };
} // namespace

#endif // UNO_MA27SOLVER_H