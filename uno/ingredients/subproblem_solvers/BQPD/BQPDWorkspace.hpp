// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDWORKSPACE_H
#define UNO_BQPDWORKSPACE_H

#include <cstddef>
#include "linear_algebra/Vector.hpp"
#include "../SolverWorkspace.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   class Evaluations;
   class Iterate;
   class OptimizationProblem;
   class Subproblem;
   class WarmstartInformation;

   class BQPDWorkspace: public SolverWorkspace {
   public:
      BQPDWorkspace() = default;
      ~BQPDWorkspace() override = default;

      void initialize(const Subproblem& subproblem);

      // data-driven setup: dense objective gradient + COO constraint Jacobian (row = constraint,
      // column = variable) + COO Lagrangian Hessian (one triangle; empty for an LP). Allocates native
      // storage and converts the Jacobian to BQPD's packed weak-CSR layout; stores the Hessian as COO.
      void set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values);

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void evaluate_functions(const OptimizationProblem& problem, const Iterate& current_iterate, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information);

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

   protected:
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

#endif // UNO_BQPDWORKSPACE_H