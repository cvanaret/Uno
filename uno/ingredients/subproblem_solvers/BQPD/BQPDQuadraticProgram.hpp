// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDQUADRATICPROGRAM_H
#define UNO_BQPDQUADRATICPROGRAM_H

#include <functional>
#include <vector>
#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class Iterate;
   class OptimizationProblem;
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   // BQPD-native QuadraticProgram. The objective gradient and the (weak-CSR) constraint Jacobian are
   // packed into a single "gradients" array held by the workspace; the variable and constraint bounds
   // are concatenated in lower_bounds/upper_bounds; the Lagrangian Hessian is exposed to BQPD's gdotx
   // callback through compute_hessian_vector_product(), which is Subproblem-free (it uses either the explicit
   // regularized Hessian stored in the workspace or an operator captured by build()).
   class BQPDQuadraticProgram : public QuadraticProgram, public SolverWorkspace {
   public:
      BQPDQuadraticProgram() = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      // data-driven build: dense objective gradient + COO constraint Jacobian + COO Lagrangian Hessian
      // (one triangle; empty for an LP). Converts COO to BQPD's weak-CSR layout internally.
      void fill(const Vector<double>& linear_objective, const Vector<uno_int>& jacobian_row_indices,
         const Vector<uno_int>& jacobian_column_indices, const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds) override;

      // whether the QP has a Hessian: determines BQPD's kmax (0 for an LP). Valid after a build().
      [[nodiscard]] bool has_curvature() const { return this->use_explicit_hessian || bool(this->hessian_operator); }

      // result <- H * vector, called from the Fortran gdotx callback. Subproblem-free.
      void compute_hessian_vector_product(int dimension, const double* vector, double* result) const;
      
      // data-driven setup: dense objective gradient + COO constraint Jacobian (row = constraint,
      // column = variable) + COO Lagrangian Hessian (one triangle; empty for an LP). Allocates native
      // storage and converts the Jacobian to BQPD's packed weak-CSR layout; stores the Hessian as COO.
      void set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values, const Vector<uno_int>& hessian_row_indices,
         const Vector<uno_int>& hessian_column_indices, const Vector<double>& hessian_values);

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void evaluate_functions(const OptimizationProblem& problem, const Iterate& current_iterate, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information);

      // lower and upper bounds of variables and constraints (concatenated, length n + m)
      std::vector<double> lower_bounds{}, upper_bounds{};
      
      Vector<double> constraints{};
      Vector<double> gradients{};
      Vector<uno_int> gradients_sparsity{};
      // COO constraint Jacobian
      Vector<uno_int> jacobian_row_indices{};
      Vector<uno_int> jacobian_column_indices{};
      Vector<double> jacobian_values{};
      Vector<size_t> permutation_vector{};
      // COO Hessian
      Vector<uno_int> hessian_row_indices{};
      Vector<uno_int> hessian_column_indices{};
      Vector<double> hessian_values{};
      bool hessian_evaluation_required{false};
      mutable Vector<double> hessian_vector_product{};

   private:
      // Hessian representation selected at build() time
      bool use_explicit_hessian{false};
      std::function<void(const double* vector, double* result)> hessian_operator{};
      
      // allocate native storage for the given problem shape
      void allocate_memory(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, bool allocate_explicit_hessian);
      // build BQPD's packed weak-CSR Jacobian sparsity (header + column indices + row-start pointers) and the
      // sorting permutation from the COO arrays already stored in jacobian_row_indices/jacobian_column_indices
      void build_gradients_sparsity_from_jacobian_coo(size_t number_variables, size_t number_constraints);
      // scatter jacobian_values into gradients[number_variables ..] using the sorting permutation
      void scatter_jacobian_values(size_t number_variables);

      void compute_gradients_sparsity(const Subproblem& subproblem);
      void evaluate_jacobian(const OptimizationProblem& problem, const Vector<double>& primals, Evaluations& evaluations);
   };
} // namespace

#endif // UNO_BQPDQUADRATICPROGRAM_H
