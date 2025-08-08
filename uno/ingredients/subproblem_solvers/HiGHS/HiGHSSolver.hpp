// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "Highs.h"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public QPSolver {
   public:
      explicit HiGHSSolver(const Options& options);

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override;

      void evaluate_constraint_jacobian(const Subproblem& subproblem) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

   protected:
      HighsModel model;
      Highs highs_solver;
      std::vector<double> constraints{};
      Vector<double> linear_objective{};
      // constraint Jacobian in COO format
      Vector<size_t> jacobian_row_indices{};
      Vector<size_t> jacobian_column_indices{};
      // Lagrangian Hessian in COO format
      Vector<size_t> hessian_row_indices{};
      Vector<size_t> hessian_column_indices{};
      Vector<double> hessian_values{};
      Vector<size_t> permutation_vector{};

      const bool print_subproblem;

      void compute_hessian_sparsity(const Subproblem& subproblem);
      void set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, const WarmstartInformation& warmstart_information);
      void solve_subproblem(const Subproblem& subproblem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H