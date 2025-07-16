// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONMODEL_H
#define UNO_PYTHONMODEL_H

#include <vector>
#include "PythonTypes.hpp"
#include "UnoSolverWrapper.hpp"
#include "model/Model.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace uno {
   // forward reference
   class Options;

   class PythonModel: public Model {
   public:
      PythonModel(const std::string& file_name, size_t number_variables, size_t number_constraints,
         double objective_sign, const objective_function_type& evaluate_objective, const constraint_functions_type& evaluate_constraints,
         const objective_gradient_type& evaluate_objective_gradient, const jacobian_type& evaluate_jacobian,
         const lagrangian_hessian_type& evaluate_lagrangian_hessian, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const std::vector<double>& primal_initial_point,
         const std::vector<double>& dual_initial_point);
      ~PythonModel() override = default;

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override;
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override;
      void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
      void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const override;
      void compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const override;

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
      // functions
      const objective_function_type& objective;
      const constraint_functions_type& constraints;
      const objective_gradient_type& objective_gradient;
      const jacobian_type& jacobian;
      const lagrangian_hessian_type& hessian;

      const std::vector<double>& variables_lower_bounds;
      const std::vector<double>& variables_upper_bounds;
      const std::vector<double>& constraints_lower_bounds;
      const std::vector<double>& constraints_upper_bounds;

      const std::vector<double>& primal_initial_point;
      const std::vector<double>& dual_initial_point;

      std::vector<BoundType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
      std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
      std::vector<BoundType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
    * UNBOUNDED) */

      // lists of variables and constraints + corresponding collection objects
      std::vector<size_t> linear_constraints{};
      CollectionAdapter<std::vector<size_t>&> linear_constraints_collection;
      std::vector<size_t> equality_constraints{};
      CollectionAdapter<std::vector<size_t>&> equality_constraints_collection;
      std::vector<size_t> inequality_constraints{};
      CollectionAdapter<std::vector<size_t>&> inequality_constraints_collection;
      SparseVector<size_t> slacks{};
      std::vector<size_t> lower_bounded_variables;
      CollectionAdapter<std::vector<size_t>&> lower_bounded_variables_collection;
      std::vector<size_t> upper_bounded_variables;
      CollectionAdapter<std::vector<size_t>&> upper_bounded_variables_collection;
      std::vector<size_t> single_lower_bounded_variables{}; // indices of the single lower-bounded variables
      CollectionAdapter<std::vector<size_t>&> single_lower_bounded_variables_collection;
      std::vector<size_t> single_upper_bounded_variables{}; // indices of the single upper-bounded variables
      CollectionAdapter<std::vector<size_t>&> single_upper_bounded_variables_collection;
      Vector<size_t> fixed_variables;

   protected:
      void categorize_bounded_variables();
   };
} // namespace

#endif // UNO_PYTHONMODEL_H