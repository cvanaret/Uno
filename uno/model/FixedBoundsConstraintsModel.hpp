// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H
#define UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H

#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // forward declaration
   class Options;

   // move the fixed variables to the set of general constraints (required in barrier methods)
   // for instance, the constraint x[2] == 1 is not interpreted as 1 <= x[2] <= 1, but as the
   // linear constraint 0 x[1] + 1 x[2] + ... + 0 x[n] == 1
   class FixedBoundsConstraintsModel: public Model {
   public:
      explicit FixedBoundsConstraintsModel(const Model& original_model);

      // availability of linear operators
      [[nodiscard]] bool has_jacobian_operator() const override;
      [[nodiscard]] bool has_jacobian_transposed_operator() const override;
      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      // sparsity patterns of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const override;
      void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const override;

      // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
      void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& multipliers, double* result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;

      void postprocess_solution(Iterate& iterate) const override;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   private:
      const Model& model;
      Vector<size_t> fixed_variables{};
      Concatenation<const Collection<size_t>&, ForwardRange> equality_constraints;
      Concatenation<const Collection<size_t>&, ForwardRange> linear_constraints;
   };
} // namespace

#endif // UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H