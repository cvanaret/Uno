// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NLMODEL_H
#define UNO_NLMODEL_H

#include <vector>
#include "model/Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/ProblemType.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "tools/NumberModelEvaluations.hpp"
#include "cppasl/nl_model.hpp"

namespace uno {
   // forward declaration
   class Result;

   // Model read from an .nl file by cppasl
   class NLModel: public Model {
   public:
      explicit NLModel(const std::string& file_name);

      // same convention as AMPLModel: ∇²L(x, y) = σ ∇²f(x) - Σ_j y_j ∇²c_j(x)
      static constexpr double lagrangian_sign_convention{-1.};

      [[nodiscard]] ProblemType get_problem_type() const override;

      [[nodiscard]] bool has_jacobian_operator() const override;
      [[nodiscard]] bool has_jacobian_transposed_operator() const override;
      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override;
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const override;
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const override;
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override;

      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const override;

      void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const override;

      [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_nonlinear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;
      void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const override;

      void write_solution_to_file(Result& result) const;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

      [[nodiscard]] size_t number_model_objective_evaluations() const override;
      [[nodiscard]] size_t number_model_constraints_evaluations() const override;
      [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override;
      [[nodiscard]] size_t number_model_jacobian_evaluations() const override;
      [[nodiscard]] size_t number_model_hessian_evaluations() const override;
      void reset_number_evaluations() const override;

   private:
      NLModel(const std::string& file_name, cppasl::NlModel&& nl_model);

      const cppasl::NlModel nl_model;
      // single workspace (not thread-safe): one per thread would be needed for parallel evaluations
      mutable cppasl::EvaluationWorkspace workspace;
      const bool has_objective;

      Vector<uno_int> jacobian_row_indices;
      Vector<uno_int> jacobian_column_indices;
      mutable std::vector<double> jacobian_values; // for the Jacobian operators
      mutable std::vector<double> cppasl_multipliers; // cppasl uses σ ∇²f + Σ_j y_j ∇²c_j: we pass -y

      std::vector<size_t> linear_constraints{};
      CollectionAdapter<std::vector<size_t>&> linear_constraints_collection;
      std::vector<size_t> nonlinear_constraints{};
      CollectionAdapter<std::vector<size_t>&> nonlinear_constraints_collection;
      std::vector<size_t> equality_constraints{};
      CollectionAdapter<std::vector<size_t>&> equality_constraints_collection;
      std::vector<size_t> inequality_constraints{};
      CollectionAdapter<std::vector<size_t>&> inequality_constraints_collection;
      Vector<size_t> fixed_variables;
      ProblemType problem_type{ProblemType::NONLINEAR};

      mutable NumberModelEvaluations number_model_evaluations{};

      [[nodiscard]] ProblemType determine_problem_type() const;
      const double* negated_multipliers(const Vector<double>& multipliers) const;
   };
} // namespace

#endif // UNO_NLMODEL_H
