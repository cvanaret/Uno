// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONMODEL_H
#define UNO_PYTHONMODEL_H

#include <vector>
#include "PythonTypes.hpp"
#include "UnoSolverWrapper.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "model/Model.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace uno {
   // forward reference
   class Options;

   class PythonModel: public Model {
   public:
      PythonModel(const std::string& file_name, size_t number_variables, size_t number_constraints,
         double objective_sign, const objective_function_type& objective, const constraint_functions_type& constraints,
         const objective_gradient_type& evaluate_objective_gradient, const jacobian_type& evaluate_jacobian,
         const lagrangian_hessian_type& evaluate_lagrangian_hessian, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const std::vector<double>& primal_initial_point,
         const std::vector<double>& dual_initial_point);
      ~PythonModel() override = default;

      // Hessian representation
      [[nodiscard]] bool has_implicit_hessian_representation() const override;
      [[nodiscard]] bool has_explicit_hessian_representation() const override;

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      // sparsity patterns of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_format) const override;
      void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const override;
      // here we use pointers, since the vector and the result may be provided by a low-level subproblem solver
      void compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const override;

      // purely functions
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
      void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const override;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   protected:
      // functions
      const objective_function_type& objective;
      const constraint_functions_type& constraints;
      const objective_gradient_type& objective_gradient;
      const jacobian_type& jacobian;
      const lagrangian_hessian_type& hessian;

      const size_t jacobian_nnz;
      const size_t hessian_nnz;

      const std::vector<double>& variables_lower_bounds;
      const std::vector<double>& variables_upper_bounds;
      const std::vector<double>& constraints_lower_bounds;
      const std::vector<double>& constraints_upper_bounds;

      const std::vector<double>& primal_initial_point;
      const std::vector<double>& dual_initial_point;

      std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */

      // lists of variables and constraints + corresponding collection objects
      ForwardRange linear_constraints{0};
      std::vector<size_t> equality_constraints{};
      CollectionAdapter<std::vector<size_t>&> equality_constraints_collection;
      std::vector<size_t> inequality_constraints{};
      CollectionAdapter<std::vector<size_t>&> inequality_constraints_collection;
      SparseVector<size_t> slacks{};
      Vector<size_t> fixed_variables;
   };
} // namespace

#endif // UNO_PYTHONMODEL_H