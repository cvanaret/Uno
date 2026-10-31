// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALEDMODEL_H
#define UNO_SCALEDMODEL_H

#include "Model.hpp"
#include "optimization/Scaling.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class Options;

   class ScaledModel: public Model {
   public:
      ScaledModel(const Model& original_model, const Vector<double>& initial_primals, const Options& options);
      // identity scaling, until the constraints are rescaled
      ScaledModel(const Model& original_model, const Options& options);

      // recompute the scaling of the constraints from the Jacobian of the original model at a given point
      void rescale_constraints(const Vector<double>& original_jacobian_values);
      [[nodiscard]] const std::vector<double>& get_constraint_scaling() const;

      [[nodiscard]] ProblemType get_problem_type() const override;

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
      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const override;
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const override;
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const override;

      // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
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

      [[nodiscard]] double unscaled_objective(double objective_value) const override;
      [[nodiscard]] double unscaled_stationarity(double scaled_stationarity) const override;
      [[nodiscard]] double unscaled_complementarity(double scaled_complementarity) const override;
      [[nodiscard]] double unscaled_primal_infeasibility(double scaled_primal_infeasibility,
         const Vector<double>& scaled_constraints, Norm residual_norm) const override;

   private:
      const Model& model;
      Scaling scaling;
      mutable Vector<double> scaled_multipliers{};
      std::vector<double> constraints_lower_bounds{};
      std::vector<double> constraints_upper_bounds{};
      mutable Vector<double> unscaled_constraints;

      void scale_constraint_bounds();
   };
} // namespace

#endif // UNO_SCALEDMODEL_H