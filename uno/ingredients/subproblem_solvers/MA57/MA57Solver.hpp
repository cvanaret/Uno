// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <vector>
#include "ingredients/subproblem_solvers/DirectEqualityQPSolver.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   struct MA57Factorization {
      int n{};
      int nnz{};
      int lfact{};
      int lifact{};

      MA57Factorization() = default;
   };

   class MA57Solver : public DirectEqualityQPSolver<size_t, double> {
   public:
      MA57Solver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros);
      ~MA57Solver() override = default;

      void do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) override;
      void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
      [[nodiscard]] SubproblemStatus solve_equality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& initial_point, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
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

      const size_t dimension;
      const size_t number_nonzeros;

      // internal matrix representation
      std::vector<int> row_indices;
      std::vector<int> column_indices;
      SymmetricMatrix<size_t, double> augmented_matrix;
      Vector<double> rhs;
      Vector<double> solution;

      // factorization
      MA57Factorization factorization{};
      std::vector<double> fact{0}; // do not initialize, resize at every iteration
      std::vector<int> ifact{0}; // do not initialize, resize at every iteration
      const int lkeep;
      std::vector<int> keep{};
      std::vector<int> iwork{};
      int lwork;
      std::vector<double> work{};

      // for ma57id_ (default values of controlling parameters)
      std::array<double, 5> cntl{};
      std::array<int, 20> icntl{};
      std::array<double, 20> rinfo{};
      std::array<int, 40> info{};

      const int nrhs{1}; // number of right hand side being solved
      const int job{1};
      std::vector<double> residuals;
      const size_t fortran_shift{1};

      bool use_iterative_refinement{false};
      void save_sparsity_pattern_internally(const SymmetricMatrix<size_t, double>& matrix);
      void solve_indefinite_linear_system();
   };
} // namespace

#endif // UNO_MA57SOLVER_H