// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REDUCEDMODEL_H
#define UNO_REDUCEDMODEL_H

#include <vector>
#include "Model.hpp"
#include "ModelReduction.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace uno {
   // model posed in the space of the non-fixed variables, in which the bound constraints are variable bounds.
   // postprocess_solution maps the primal-dual solution back to the original space
   class ReducedModel: public Model {
   public:
      ReducedModel(const Model& original_model, ModelReduction&& reduction);
      // the constraint collections refer to members of the object
      ReducedModel(const ReducedModel&) = delete;
      ReducedModel& operator=(const ReducedModel&) = delete;

      [[nodiscard]] ProblemType get_problem_type() const override;

      // availability of linear operators
      [[nodiscard]] bool has_jacobian_operator() const override;
      [[nodiscard]] bool has_jacobian_transposed_operator() const override;
      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override;
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      // sparsity patterns of Jacobian and Hessian
      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const override;
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const override;
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const override;

      // linear operators
      void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const override;

      [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override;
      [[nodiscard]] const std::vector<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const std::vector<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_nonlinear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;
      void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const override;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

      [[nodiscard]] size_t number_model_objective_evaluations() const override;
      [[nodiscard]] size_t number_model_constraints_evaluations() const override;
      [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override;
      [[nodiscard]] size_t number_model_jacobian_evaluations() const override;
      [[nodiscard]] size_t number_model_hessian_evaluations() const override;
      void reset_number_evaluations() const override;

   private:
      const Model& model;
      const ModelReduction reduction;
      std::vector<double> variables_lower_bounds{}, variables_upper_bounds{};
      std::vector<double> constraints_lower_bounds{}, constraints_upper_bounds{};
      std::vector<size_t> equality_constraints{}, inequality_constraints{};
      std::vector<size_t> linear_constraints{}, nonlinear_constraints{};
      const CollectionAdapter<std::vector<size_t>> linear_constraints_collection{this->linear_constraints};
      const CollectionAdapter<std::vector<size_t>> nonlinear_constraints_collection{this->nonlinear_constraints};
      // nonzeros of the original Jacobian and Hessian that are kept, and their reduced (0-based) indices
      std::vector<size_t> jacobian_nonzeros{}, hessian_nonzeros{};
      std::vector<uno_int> jacobian_row_indices{}, jacobian_column_indices{};
      std::vector<uno_int> hessian_row_indices{}, hessian_column_indices{};
      // buffers in the original space. The fixed components of the primals and the multipliers of the bound
      // constraints are set once and for all
      mutable Vector<double> original_primals, original_multipliers, original_constraints, original_gradient;
      mutable Vector<double> original_vector, original_result;
      mutable Vector<double> original_jacobian_values, original_hessian_values;

      template <typename Array>
      const Vector<double>& expand_primals(const Array& x) const;
      void recover_multipliers(const Iterate& reduced_iterate, Iterate& iterate, Evaluations& evaluations) const;
   };
} // namespace

#endif // UNO_REDUCEDMODEL_H
