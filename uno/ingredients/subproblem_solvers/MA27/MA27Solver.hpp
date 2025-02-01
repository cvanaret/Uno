// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/convexification_strategies/PrimalDualConvexificationStrategy.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   // forward declarations
   class Options;
   template <typename ElementType>
   class Vector;

   class MA27Solver: public DirectSymmetricIndefiniteLinearSolver<size_t, double> {
   public:
      explicit MA27Solver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros,
            const Options& options);
      ~MA27Solver() override = default;

      void do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) override;
      void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
      void solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) override;
      void solve_indefinite_system(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, Vector<double>& result,
            WarmstartInformation& warmstart_information) override;

      [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

   private:
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */
      SymmetricMatrix<size_t, double> hessian;

      size_t dimension;
      size_t number_nonzeros;

      // internal matrix representation
      std::vector<int> row_indices{};          // row index of input
      std::vector<int> column_indices{};          // col index of input
      SymmetricMatrix<size_t, double> augmented_matrix;
      Vector<double> rhs;
      PrimalDualConvexificationStrategy<double> primal_dual_convexification_strategy;

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

      // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
      void set_up_subproblem(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, WarmstartInformation& warmstart_information);
      void assemble_augmented_rhs(LagrangeNewtonSubproblem& subproblem);
      void save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix);
      void check_factorization_status();
   };
} // namespace

#endif // UNO_MA27SOLVER_H