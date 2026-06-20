// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <stdexcept>
#include "BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   #define BIG 1e30

   void BQPDQuadraticProgram::initialize_memory(const Subproblem& subproblem) {
      this->number_variables = subproblem.number_variables;
      this->number_constraints = subproblem.number_constraints;
      this->number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      if (subproblem.has_curvature() && !subproblem.has_hessian_operator() && !subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The Hessian cannot be evaluated implicitly or explicitly");
      }
      this->workspace.initialize(subproblem);
      this->lower_bounds.resize(this->number_variables + this->number_constraints);
      this->upper_bounds.resize(this->number_variables + this->number_constraints);
   }

   void BQPDQuadraticProgram::build(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // evaluate objective gradient + constraint Jacobian (packed into workspace.gradients)
      this->workspace.evaluate_functions(subproblem.problem, subproblem.current_iterate, current_evaluations,
         warmstart_information);

      // variable bounds
      if (warmstart_information.trust_region_changed) {
         subproblem.set_variables_bounds(this->lower_bounds, this->upper_bounds, trust_region_radius);
      }
      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.new_iterate) {
         auto constraints_lower_bounds = view(this->lower_bounds, this->number_variables,
            this->number_variables + this->number_constraints);
         auto constraints_upper_bounds = view(this->upper_bounds, this->number_variables,
            this->number_variables + this->number_constraints);
         subproblem.set_constraints_bounds(constraints_lower_bounds, constraints_upper_bounds, this->workspace.constraints);
      }
      // replace INFs with large finite values (TODO: is that really useful?)
      for (size_t index: Range(this->number_variables + this->number_constraints)) {
         this->lower_bounds[index] = std::max(-BIG, this->lower_bounds[index]);
         this->upper_bounds[index] = std::min(BIG, this->upper_bounds[index]);
      }

      // choose the Hessian representation for the gdotx callback. This mirrors the logic formerly in
      // the callback: use the explicit matrix if the Hessian must be regularized or no operator exists.
      this->use_explicit_hessian = (!subproblem.is_hessian_positive_definite() && subproblem.performs_primal_regularization()) ||
         !subproblem.has_hessian_operator();
      if (this->use_explicit_hessian) {
         if (!subproblem.has_hessian_matrix()) {
            throw std::runtime_error("The Lagrangian Hessian has no appropriate representation");
         }
         // evaluate + regularize the explicit Hessian once per iterate
         if (this->workspace.hessian_evaluation_required) {
            subproblem.evaluate_lagrangian_hessian(statistics, this->workspace.hessian_values.data());
            subproblem.regularize_lagrangian_hessian(statistics, this->workspace.hessian_values.data());
            this->workspace.hessian_evaluation_required = false;
         }
         this->hessian_operator = nullptr;
      }
      else {
         // matrix-free Hessian-vector product evaluated at the current iterate. The Subproblem outlives
         // this build()/solve() pair (it is owned by the caller for the duration of the solve), so
         // capturing it by pointer is safe; the operator is rebuilt on every build().
         const Subproblem* subproblem_pointer = &subproblem;
         this->hessian_operator = [subproblem_pointer](const double* vector, double* result) {
            subproblem_pointer->compute_hessian_vector_product(subproblem_pointer->current_iterate.primals.data(),
               vector, result);
         };
      }
   }

   void BQPDQuadraticProgram::build(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) {
      // infer the dimensions from the data
      this->number_variables = linear_objective.size();
      this->number_constraints = constraints_lower_bounds.size();

      // allocate native storage and convert the COO Jacobian to BQPD's packed weak-CSR layout; the COO
      // Hessian is stored as-is and consumed by the symmetric matvec in compute_hessian_vector_product()
      this->workspace.set_from_coo(this->number_variables, this->number_constraints, linear_objective,
         jacobian_row_indices, jacobian_column_indices, jacobian_values,
         hessian_row_indices, hessian_column_indices, hessian_values);

      // concatenated variable + constraint bounds, with INFs clamped to large finite values
      this->lower_bounds.resize(this->number_variables + this->number_constraints);
      this->upper_bounds.resize(this->number_variables + this->number_constraints);
      for (size_t variable_index: Range(this->number_variables)) {
         this->lower_bounds[variable_index] = std::max(-BIG, variables_lower_bounds[variable_index]);
         this->upper_bounds[variable_index] = std::min(BIG, variables_upper_bounds[variable_index]);
      }
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->lower_bounds[this->number_variables + constraint_index] = std::max(-BIG, constraints_lower_bounds[constraint_index]);
         this->upper_bounds[this->number_variables + constraint_index] = std::min(BIG, constraints_upper_bounds[constraint_index]);
      }

      // explicit Hessian (no operator) in the data-driven path; an empty Hessian means an LP
      this->use_explicit_hessian = (0 < hessian_values.size());
      this->hessian_operator = nullptr;
   }

   SolverWorkspace& BQPDQuadraticProgram::get_workspace() {
      return this->workspace;
   }

   void BQPDQuadraticProgram::compute_hessian_vector_product(int dimension, const double* vector, double* result) const {
      for (size_t index = 0; index < static_cast<size_t>(dimension); ++index) {
         result[index] = 0.;
      }
      if (this->use_explicit_hessian) {
         // explicit symmetric COO matvec (only the lower/upper triangle is stored)
         const size_t number_hessian_nonzeros = this->workspace.hessian_values.size();
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            const size_t row_index = static_cast<size_t>(this->workspace.hessian_row_indices[nonzero_index]);
            const size_t column_index = static_cast<size_t>(this->workspace.hessian_column_indices[nonzero_index]);
            const double entry = this->workspace.hessian_values[nonzero_index];
            result[row_index] += entry * vector[column_index];
            if (row_index != column_index) {
               result[column_index] += entry * vector[row_index];
            }
         }
      }
      else if (this->hessian_operator) {
         this->hessian_operator(vector, result);
      }
      else {
         throw std::runtime_error("The Lagrangian Hessian has no appropriate representation");
      }
   }
} // namespace
