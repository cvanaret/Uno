// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include "BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   #define BIG 1e30

   void BQPDQuadraticProgram::initialize_memory(const Subproblem& subproblem) {
      if (subproblem.has_curvature() && !subproblem.has_hessian_operator() && !subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The Hessian cannot be evaluated implicitly or explicitly");
      }
      this->number_variables = subproblem.number_variables;
      this->number_constraints = subproblem.number_constraints;
      this->number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      this->lower_bounds.resize(this->number_variables + this->number_constraints);
      this->upper_bounds.resize(this->number_variables + this->number_constraints);
      
      // allocate an explicit Hessian matrix if:
      // - the Hessian is not positive definite and must be regularized, or
      // - the Hessian model only has an explicit representation
      const bool allocate_explicit_hessian = (!subproblem.is_hessian_positive_definite() && subproblem.performs_primal_regularization()) ||
         !subproblem.has_hessian_operator();
      const size_t number_hessian_nonzeros = allocate_explicit_hessian ? subproblem.number_regularized_hessian_nonzeros() : 0;
      this->allocate_memory(subproblem.number_variables, subproblem.number_constraints, subproblem.number_jacobian_nonzeros(),
         number_hessian_nonzeros, allocate_explicit_hessian);

      // save sparsity patterns of objective gradient and constraint Jacobian into the BQPD workspace
      this->compute_gradients_sparsity(subproblem);

      if (allocate_explicit_hessian) {
         subproblem.compute_regularized_hessian_sparsity(this->hessian_row_indices.data(),
            this->hessian_column_indices.data(), Indexing::C_indexing);
      }
      if (subproblem.has_hessian_operator()) {
         this->hessian_vector_product.resize(subproblem.number_variables);
      }
   }
   
   void BQPDQuadraticProgram::allocate_memory(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, bool allocate_explicit_hessian) {
      this->constraints.resize(number_constraints);

      // Jacobian + objective gradient: gradients is a concatenation of the dense objective gradient and the
      // sparse (weak-CSR) constraint Jacobian; gradients_sparsity is BQPD's "la" array
      this->gradients.resize(number_variables + number_jacobian_nonzeros);
      this->gradients_sparsity.resize(1 + number_variables + number_jacobian_nonzeros + 1 + number_constraints + 1);
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      this->jacobian_values.resize(number_jacobian_nonzeros);

      if (allocate_explicit_hessian) {
         this->hessian_row_indices.resize(number_hessian_nonzeros);
         this->hessian_column_indices.resize(number_hessian_nonzeros);
         this->hessian_values.resize(number_hessian_nonzeros);
      }
   }

   void BQPDQuadraticProgram::fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // evaluate objective gradient + constraint Jacobian (packed into workspace.gradients)
      this->evaluate_functions(subproblem.problem, subproblem.current_iterate, current_evaluations,
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
         subproblem.set_constraints_bounds(constraints_lower_bounds, constraints_upper_bounds, this->constraints);
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
         if (this->hessian_evaluation_required) {
            subproblem.evaluate_lagrangian_hessian(statistics, this->hessian_values.data());
            subproblem.regularize_lagrangian_hessian(statistics, this->hessian_values.data());
            this->hessian_evaluation_required = false;
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

   void BQPDQuadraticProgram::fill(const Vector<double>& linear_objective,
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
      this->set_from_coo(this->number_variables, this->number_constraints, linear_objective,
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

   void BQPDQuadraticProgram::compute_hessian_vector_product(int dimension, const double* vector, double* result) const {
      for (size_t index = 0; index < static_cast<size_t>(dimension); ++index) {
         result[index] = 0.;
      }
      if (this->use_explicit_hessian) {
         // explicit symmetric COO matvec (only the lower/upper triangle is stored)
         const size_t number_hessian_nonzeros = this->hessian_values.size();
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            const size_t row_index = static_cast<size_t>(this->hessian_row_indices[nonzero_index]);
            const size_t column_index = static_cast<size_t>(this->hessian_column_indices[nonzero_index]);
            const double entry = this->hessian_values[nonzero_index];
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
   
   void BQPDQuadraticProgram::set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values) {
      const size_t number_jacobian_nonzeros = jacobian_values.size();
      const size_t number_hessian_nonzeros = hessian_values.size();
      const bool has_explicit_hessian = (0 < number_hessian_nonzeros);
      this->allocate_memory(number_variables, number_constraints, number_jacobian_nonzeros, number_hessian_nonzeros,
         has_explicit_hessian);

      // dense objective gradient -> gradients[0 .. number_variables)
      for (size_t variable_index: Range(number_variables)) {
         this->gradients[variable_index] = linear_objective[variable_index];
      }

      // constraint Jacobian: copy the COO sparsity + values, build the weak-CSR sparsity, scatter the values
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         this->jacobian_row_indices[nonzero_index] = jacobian_row_indices[nonzero_index];
         this->jacobian_column_indices[nonzero_index] = jacobian_column_indices[nonzero_index];
         this->jacobian_values[nonzero_index] = jacobian_values[nonzero_index];
      }
      this->build_gradients_sparsity_from_jacobian_coo(number_variables, number_constraints);
      this->scatter_jacobian_values(number_variables);

      // Lagrangian Hessian: stored as COO, consumed directly by the symmetric matvec in
      // BQPDQuadraticProgram::compute_hessian_vector_product()
      for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
         this->hessian_row_indices[nonzero_index] = hessian_row_indices[nonzero_index];
         this->hessian_column_indices[nonzero_index] = hessian_column_indices[nonzero_index];
         this->hessian_values[nonzero_index] = hessian_values[nonzero_index];
      }
      this->hessian_evaluation_required = false;
   }

   double BQPDQuadraticProgram::compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const {
      if (subproblem.has_hessian_operator()) { // linear operator
         // TODO compute the quadratic form directly without temporary result
         // compute Hv
         subproblem.compute_hessian_vector_product(subproblem.current_iterate.primals.data(), vector.data(),
            this->hessian_vector_product.data());
         // compute the dot product <v, Hv>
         return dot(view(vector, 0, subproblem.number_variables), this->hessian_vector_product);
      }
      else if (subproblem.has_hessian_matrix()) { // explicit matrix
         double quadratic_product = 0.;
         const size_t number_hessian_nonzeros = this->hessian_values.size();
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            const size_t row_index = static_cast<size_t>(this->hessian_row_indices[nonzero_index]);
            const size_t column_index = static_cast<size_t>(this->hessian_column_indices[nonzero_index]);
            const double entry = this->hessian_values[nonzero_index];
            if (row_index >= vector.size() || column_index >= vector.size()) {
               throw std::runtime_error("Dimension mismatch");
            }

            const double factor = (row_index != column_index) ? 2. : 1.;
            quadratic_product += factor * entry * vector[row_index] * vector[column_index];
         }
         return quadratic_product;
      }
      else {
         throw std::runtime_error("The Lagrangian Hessian has no appropriate representation");
      }
   }

   void BQPDQuadraticProgram::evaluate_functions(const OptimizationProblem& problem, const Iterate& current_iterate,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // evaluate the functions based on warmstart information
      // gradients is a concatenation of the dense objective gradient and the sparse Jacobian
      if (warmstart_information.new_iterate) {
         this->gradients.fill(0.);
         problem.evaluate_objective_gradient(current_iterate, this->gradients.data(), current_evaluations);
         problem.evaluate_constraints(current_iterate, this->constraints.data(), current_evaluations);
         this->evaluate_jacobian(problem, current_iterate.primals, current_evaluations);
         this->hessian_evaluation_required = true;
      }
   }

   // objective gradient + row-major constraint Jacobian
   void BQPDQuadraticProgram::compute_gradients_sparsity(const Subproblem& subproblem) {
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      // get the Jacobian sparsity in COO format (row = constraint, column = variable)
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      subproblem.compute_jacobian_sparsity(this->jacobian_row_indices.data(), this->jacobian_column_indices.data(),
         0, 0, Indexing::C_indexing, MatrixOrder::ROW_MAJOR);
      // convert COO -> BQPD's packed weak-CSR layout
      this->build_gradients_sparsity_from_jacobian_coo(subproblem.number_variables, subproblem.number_constraints);
   }

   // build BQPD's "la" sparsity array (gradients_sparsity) and the sorting permutation from the COO Jacobian
   // already present in jacobian_row_indices/jacobian_column_indices. Shared by the Subproblem path and the
   // data-driven path so that both produce exactly the same native layout.
   void BQPDQuadraticProgram::build_gradients_sparsity_from_jacobian_coo(size_t number_variables, size_t number_constraints) {
      const size_t number_jacobian_nonzeros = this->jacobian_row_indices.size();

      // header
      const size_t position_of_row_starts = 1 + number_variables + number_jacobian_nonzeros;
      this->gradients_sparsity[0] = static_cast<int>(position_of_row_starts);

      // dense objective gradient
      for (size_t variable_index: Range(number_variables)) {
         this->gradients_sparsity[1 + variable_index] = static_cast<int>(variable_index + Indexing::Fortran_indexing);
      }

      this->gradients_sparsity[position_of_row_starts] = 1; // always starts at 1
      this->gradients_sparsity[position_of_row_starts + 1] = static_cast<int>(1 + number_variables); // dense objective gradient

      // BQPD (sparse) requires a (weak) CSR Jacobian: the entries should be in increasing constraint indices.
      // Since the COO format does not require this, we need to convert from COO to CSR by permutating the entries. To
      // this end, we compute a permutation vector once and for all that contains the correct ordering of terms.
      // The permutation vector is initially filled with [0, 1, ..., number_jacobian_nonzeros-1]
      this->permutation_vector.resize(number_jacobian_nonzeros);
      std::iota(this->permutation_vector.begin(), this->permutation_vector.end(), 0);
      // sort the permutation vector such that the row indices (constraints) of the Jacobian sparsity are in increasing order
      // see https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
      std::sort(this->permutation_vector.begin(), this->permutation_vector.end(),
          [&](const size_t& i, const size_t& j) {
             return (this->jacobian_row_indices[i] < this->jacobian_row_indices[j]);
          }
      );

      // copy the COO format into BQPD's CSR format
      int current_constraint = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t permuted_nonzero_index = this->permutation_vector[jacobian_nonzero_index];
         // variable index
         const uno_int variable_index = this->jacobian_column_indices[permuted_nonzero_index];
         this->gradients_sparsity[1 + number_variables + jacobian_nonzero_index] = variable_index +
            Indexing::Fortran_indexing;

         // constraint index
         const uno_int constraint_index = this->jacobian_row_indices[permuted_nonzero_index];
         if (current_constraint > constraint_index) {
            throw std::runtime_error("Dimension mismatch");
         }
         while (current_constraint < constraint_index) {
            ++current_constraint;
            this->gradients_sparsity[1 + number_variables + number_jacobian_nonzeros + 1 +
               static_cast<size_t>(current_constraint)] = static_cast<int>(number_variables +
                  jacobian_nonzero_index) + Indexing::Fortran_indexing;
         }
      }
      // since there cannot be empty rows, we don't need to loop over empty rows like we do for the HiGHS Hessian
      this->gradients_sparsity[1 + number_variables + number_jacobian_nonzeros + 1 + number_constraints] =
         static_cast<int>(number_variables + number_jacobian_nonzeros) + Indexing::Fortran_indexing;
   }

   void BQPDQuadraticProgram::scatter_jacobian_values(size_t number_variables) {
      // copy the Jacobian values with permutation into &this->gradients[number_variables]
      const size_t number_jacobian_nonzeros = this->jacobian_values.size();
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t permuted_nonzero_index = this->permutation_vector[nonzero_index];
         this->gradients[number_variables + nonzero_index] = this->jacobian_values[permuted_nonzero_index];
      }
   }

   void BQPDQuadraticProgram::evaluate_jacobian(const OptimizationProblem& problem, const Vector<double>& primals,
         Evaluations& evaluations) {
      problem.evaluate_jacobian(primals, this->jacobian_values.data(), evaluations);
      this->scatter_jacobian_values(problem.number_variables);
   }
} // namespace
