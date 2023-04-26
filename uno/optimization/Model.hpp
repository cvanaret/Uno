// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODEL_H
#define UNO_MODEL_H

#include <string>
#include <vector>
#include <map>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "TerminationStatus.hpp"
#include "EvaluationErrors.hpp"

struct Interval {
   double lb;
   double ub;
};

enum FunctionType {
   LINEAR = 0, /*!< Linear function */
   QUADRATIC, /*!< Quadratic function */
   NONLINEAR /*!< Nonlinear function */
};

enum BoundType { EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED };

// forward declaration
class Iterate;

/*! \class Problem
 * \brief Optimization problem
 *
 *  Description of an optimization problem
 */
class Model {
public:
   Model(std::string name, size_t number_variables, size_t number_constraints, FunctionType problem_type);
   virtual ~Model() = default;

   std::string name;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */
   FunctionType problem_type;

   // objective
   double objective_sign{1.}; /*!< Sign of the objective function (1: minimization, -1: maximization) */

   // data structures to access certain types of variables and constraints
   std::vector<size_t> equality_constraints{};
   std::vector<size_t> inequality_constraints{};
   SparseVector<size_t> slacks;
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables
   std::vector<size_t> single_lower_bounded_variables{}; // indices of the single lower-bounded variables
   std::vector<size_t> single_upper_bounded_variables{}; // indices of the single upper-bounded variables

   // Hessian
   const bool fixed_hessian_sparsity{true};

   // purely virtual functions
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;

   [[nodiscard]] virtual BoundType get_variable_bound_type(size_t i) const = 0;
   [[nodiscard]] virtual FunctionType get_constraint_type(size_t j) const = 0;
   [[nodiscard]] virtual BoundType get_constraint_bound_type(size_t j) const = 0;

   [[nodiscard]] virtual size_t get_number_objective_gradient_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_number_jacobian_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_number_hessian_nonzeros() const = 0;

   [[nodiscard]] virtual double evaluate_objective(const std::vector<double>& x) const = 0;
   virtual void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const = 0;
   virtual void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const = 0;
   virtual void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const = 0;

   virtual void get_initial_primal_point(std::vector<double>& x) const = 0;
   virtual void get_initial_dual_point(std::vector<double>& multipliers) const = 0;
   virtual void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const = 0;

   // constraints
   [[nodiscard]] virtual const std::vector<size_t>& get_linear_constraints() const = 0;

   // auxiliary functions
   static void determine_bounds_types(std::vector<Interval>& variables_bounds, std::vector<BoundType>& status);
   void project_primals_onto_bounds(std::vector<double>& x) const;
   [[nodiscard]] bool is_constrained() const;
   // constraint violation
   [[nodiscard]] double compute_constraint_lower_bound_violation(double constraint_value, size_t j) const;
   [[nodiscard]] double compute_constraint_upper_bound_violation(double constraint_value, size_t j) const;
   [[nodiscard]] virtual double compute_constraint_violation(double constraint_value, size_t j) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const;

protected:
   size_t number_objective_gradient_nonzeros{0}; /*!< Number of nonzero elements in the objective gradient */
   size_t number_jacobian_nonzeros{0}; /*!< Number of nonzero elements in the constraint Jacobian */
   size_t number_hessian_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */
   void determine_constraints();
};

std::string type_to_string(FunctionType function_type);

#endif // UNO_MODEL_H
