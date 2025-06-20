// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   class MA27Solver: public DirectSymmetricIndefiniteLinearSolver<size_t, double> {
   public:
      MA27Solver();
      ~MA27Solver() override = default;

      void initialize_memory(size_t number_variables, size_t number_constraints, size_t number_hessian_nonzeros,
         size_t regularization_size) override;

      void do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) override;
      void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
      void solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) override;
      void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Vector<double>& solution,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

   private:
      int n{};                         // dimension of current factorisation (maximal value here is <= max_dimension)
      int nnz{};                     // number of nonzeros of current factorisation
      std::array<int, 30> icntl{};      // integer array of length 30; integer control values
      std::array<double, 5> cntl{};     // double array of length 5; double control values

      // internal matrix representation
      std::vector<int> row_indices{};          // row index of input
      std::vector<int> column_indices{};          // col index of input

      // evaluations
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */
      SymmetricMatrix<size_t, double> augmented_matrix{};
      Vector<double> rhs{};

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


      // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
      void save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix);
      void check_factorization_status();
   };
} // namespace

#endif // UNO_MA27SOLVER_H