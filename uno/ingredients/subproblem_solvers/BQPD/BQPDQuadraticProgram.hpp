// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDQUADRATICPROGRAM_H
#define UNO_BQPDQUADRATICPROGRAM_H

#include <functional>
#include <vector>
#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "BQPDWorkspace.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   // BQPD-native QuadraticProgram. The objective gradient and the (weak-CSR) constraint Jacobian are
   // packed into a single "gradients" array held by the workspace; the variable and constraint bounds
   // are concatenated in lower_bounds/upper_bounds; the Lagrangian Hessian is exposed to BQPD's gdotx
   // callback through compute_hessian_vector_product(), which is Subproblem-free (it uses either the explicit
   // regularized Hessian stored in the workspace or an operator captured by build()).
   class BQPDQuadraticProgram : public QuadraticProgram {
   public:
      BQPDQuadraticProgram() = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      // data-driven build: dense objective gradient + COO constraint Jacobian + COO Lagrangian Hessian
      // (one triangle; empty for an LP). Converts COO to BQPD's weak-CSR layout internally.
      void fill(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) override;

      // whether the QP has a Hessian: determines BQPD's kmax (0 for an LP). Valid after a build().
      [[nodiscard]] bool has_curvature() const { return this->use_explicit_hessian || bool(this->hessian_operator); }

      // result <- H * vector, called from the Fortran gdotx callback. Subproblem-free.
      void compute_hessian_vector_product(int dimension, const double* vector, double* result) const;

      SolverWorkspace& get_workspace();

      BQPDWorkspace workspace{};
      // lower and upper bounds of variables and constraints (concatenated, length n + m)
      std::vector<double> lower_bounds{}, upper_bounds{};

   private:
      // Hessian representation selected at build() time
      bool use_explicit_hessian{false};
      std::function<void(const double* vector, double* result)> hessian_operator{};
   };
} // namespace

#endif // UNO_BQPDQUADRATICPROGRAM_H
