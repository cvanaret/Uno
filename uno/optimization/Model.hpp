// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_MODEL_H
#define UNO_MODEL_H

#include <string>
#include <vector>
#include <map>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

struct Interval {
   double lb;
   double ub;
};

enum FunctionType {
   LINEAR = 0, /*!< Linear function */
   QUADRATIC, /*!< Quadratic function */
   NONLINEAR /*!< Nonlinear function */
};

enum ConstraintType { EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED };

enum ConstraintFeasibility { FEASIBLE, INFEASIBLE_LOWER, INFEASIBLE_UPPER };

struct NumericalError : public std::exception {
   [[nodiscard]] const char* what() const throw() override = 0;
};

struct HessianNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() override {
      return "A numerical error was encountered while evaluating the Hessian";
   }
};

struct GradientNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() override {
      return "A numerical error was encountered while evaluating a gradient";
   }
};

struct FunctionNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() override {
      return "A numerical error was encountered while evaluating a function";
   }
};

/*! \class Problem
 * \brief Optimization problem
 *
 *  Description of an optimization problem
 */
class Model {
public:
   Model(std::string name, size_t number_variables, size_t number_constraints, FunctionType problem_type);
   virtual ~Model() = default;

   static std::map<FunctionType, std::string> type_to_string;

   std::string name;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */
   FunctionType problem_type;

   // objective
   double objective_sign{1.}; /*!< Sign of the objective function (1: minimization, -1: maximization) */
   FunctionType objective_type{NONLINEAR}; /*!< Type of the objective (LINEAR, QUADRATIC, NONLINEAR) */

   SparseVector<size_t> equality_constraints; /*!< inequality constraints */
   SparseVector<size_t> inequality_constraints; /*!< inequality constraints */
   SparseVector<size_t> linear_constraints;
   // lists of bounded variables
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   // Hessian
   const bool fixed_hessian_sparsity{false};

   // purely virtual functions
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;

   [[nodiscard]] virtual ConstraintType get_variable_status(size_t i) const = 0;
   [[nodiscard]] virtual FunctionType get_constraint_type(size_t j) const = 0;
   [[nodiscard]] virtual ConstraintType get_constraint_status(size_t j) const = 0;
   [[nodiscard]] virtual size_t get_maximum_number_hessian_nonzeros() const = 0;

   [[nodiscard]] virtual double evaluate_objective(const std::vector<double>& x) const = 0;
   virtual void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const = 0;
   virtual void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const = 0;
   virtual void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix& hessian) const = 0;

   virtual void get_initial_primal_point(std::vector<double>& x) const = 0;
   virtual void get_initial_dual_point(std::vector<double>& multipliers) const = 0;

   // auxiliary functions
   static void determine_bounds_types(std::vector<Interval>& variables_bounds, std::vector<ConstraintType>& status);
   void project_point_onto_bounds(std::vector<double>& x) const;
   [[nodiscard]] bool is_constrained() const;
   // constraint violation
   [[nodiscard]] double compute_constraint_lower_bound_violation(double constraint, size_t j) const;
   [[nodiscard]] double compute_constraint_upper_bound_violation(double constraint, size_t j) const;
   [[nodiscard]] virtual double compute_constraint_violation(double constraint, size_t j) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set,
         Norm residual_norm) const;
   [[nodiscard]] double compute_complementarity_error(const std::vector<double>& x, const std::vector<double>& constraints,
         const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
         const std::vector<double>& upper_bounds_multipliers) const;

protected:
   size_t hessian_maximum_number_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */
   void determine_constraints();
};

#endif // UNO_MODEL_H
