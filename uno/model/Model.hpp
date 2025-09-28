// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODEL_H
#define UNO_MODEL_H

#include <string>
#include <vector>
#include "linear_algebra/MatrixOrder.hpp"
#include "linear_algebra/Norm.hpp"
#include "optimization/SolutionStatus.hpp"
#include "symbolic/VectorExpression.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   template <typename ElementType>
   class SparseVector;
   template <typename ElementType>
   class Vector;

   enum FunctionType {LINEAR, NONLINEAR};

   // forward declaration
   class Iterate;

   class Model {
   public:
      Model(std::string name, size_t number_variables, size_t number_constraints, double optimization_sense);
      virtual ~Model() = default;

      const std::string name;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */
      const double optimization_sense; /*!< 1: minimization, -1: maximization */

      // availability of linear operators
      [[nodiscard]] virtual bool has_jacobian_operator() const = 0;
      [[nodiscard]] virtual bool has_jacobian_transposed_operator() const = 0;
      [[nodiscard]] virtual bool has_hessian_operator() const = 0;
      [[nodiscard]] virtual bool has_hessian_matrix() const = 0;

      // function evaluations
      [[nodiscard]] virtual double evaluate_objective(const Vector<double>& x) const = 0;
      virtual void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const = 0;

      // dense objective gradient
      virtual void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const = 0;

      // sparsity patterns of Jacobian and Hessian
      virtual void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const = 0;
      virtual void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const = 0;

      // numerical evaluations of Jacobian and Hessian
      virtual void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const = 0;
      virtual void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const = 0;

      // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
      // here we use pointers, since the vector and the result may be provided by a low-level subproblem solver
      virtual void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const = 0;
      virtual void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const = 0;
      virtual void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& multipliers, double* result) const = 0;

      // purely virtual functions
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual const SparseVector<size_t>& get_slacks() const = 0;
      [[nodiscard]] virtual const Vector<size_t>& get_fixed_variables() const = 0;

      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_equality_constraints() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_linear_constraints() const = 0;

      virtual void initial_primal_point(Vector<double>& x) const = 0;
      virtual void initial_dual_point(Vector<double>& multipliers) const = 0;
      virtual void postprocess_solution(Iterate& iterate) const = 0;

      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const = 0;
      [[nodiscard]] virtual size_t number_hessian_nonzeros() const = 0;

      // auxiliary functions
      void project_onto_variable_bounds(Vector<double>& x) const;
      [[nodiscard]] bool is_constrained() const;

      // constraint violation
      [[nodiscard]] virtual double constraint_violation(double constraint_value, size_t constraint_index) const;
      template <typename Array>
      double constraint_violation(const Array& constraints, Norm residual_norm) const;

      void find_fixed_variables(Vector<size_t>& fixed_variables) const;
      void partition_constraints(std::vector<size_t>& equality_constraints, std::vector<size_t>& inequality_constraints) const;
   };

   // compute ||c||
   template <typename Array>
   double Model::constraint_violation(const Array& constraints, Norm residual_norm) const {
      const Range constraints_range = Range(constraints.size());
      const VectorExpression constraint_violation{constraints_range, [&](size_t constraint_index) {
         return this->constraint_violation(constraints[constraint_index], constraint_index);
      }};
      return norm(residual_norm, constraint_violation);
   }
} // namespace

#endif // UNO_MODEL_H