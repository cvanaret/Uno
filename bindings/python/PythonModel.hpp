// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONMODEL_H
#define UNO_PYTHONMODEL_H

#include <vector>
#include "model/Model.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace uno {
   // forward reference
   class Options;
   template <typename ElementType>
   class SparseVector;
   
   class PythonModel: public Model {
   public:
      PythonModel(const std::string& file_name, size_t number_variables, size_t number_constraints,
         double objective_sign, double (*objective)(const Vector<double>&), const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds);
      ~PythonModel() override = default;

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override;
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override;
      void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
      void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const override;
      void compute_hessian_vector_product(const Vector<double>& vector, double objective_multiplier, const Vector<double>& multipliers,
            Vector<double>& result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override;
      [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;
      void postprocess_solution(Iterate& iterate, IterateStatus iterate_status) const override;

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   protected:
      double (*objective)(const Vector<double>&);

      const std::vector<double>& variables_lower_bounds;
      const std::vector<double>& variables_upper_bounds;
      const std::vector<double>& constraints_lower_bounds;
      const std::vector<double>& constraints_upper_bounds;
   };
} // namespace

#endif // UNO_PYTHONMODEL_H