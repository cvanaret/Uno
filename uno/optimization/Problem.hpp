#ifndef PROBLEM_H
#define PROBLEM_H

#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "Constraint.hpp"
#include "linear_algebra/CSCSymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

enum FunctionType {
   LINEAR = 0, /*!< Linear function */
   QUADRATIC, /*!< Quadratic function */
   NONLINEAR /*!< Nonlinear function */
};

struct NumericalError : public std::exception {
   virtual const char* what() const throw() = 0;
};

struct HessianNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() {
      return "A numerical error was encountered while evaluating the Hessian";
   }
};

struct GradientNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() {
      return "A numerical error was encountered while evaluating a gradient";
   }
};

struct FunctionNumericalError : NumericalError {
   [[nodiscard]] const char* what() const throw() {
      return "A numerical error was encountered while evaluating a function";
   }
};

/*! \class Problem
 * \brief Optimization problem
 *
 *  Description of an optimization problem
 */
class Problem {
public:
   Problem(std::string name, int number_variables, int number_constraints, FunctionType problem_type);
   virtual ~Problem() = default;

   static std::map<FunctionType, std::string> type_to_string;

   std::string name;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */
   FunctionType type;

   /* objective */
   double objective_sign{1.}; /*!< Sign of the evaluate_objective function (1: minimization, -1: maximization) */
   FunctionType objective_type{NONLINEAR}; /*!< Type of the evaluate_objective (LINEAR, QUADRATIC, NONLINEAR) */

   /* variables */
   std::vector<Range> variables_bounds;
   std::vector<ConstraintType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */

   /* constraints */
   std::vector<Range> constraint_bounds;
   std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
   std::vector<ConstraintType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
 * UNBOUNDED) */
   std::map<size_t, size_t> equality_constraints{}; /*!< inequality constraints */
   std::map<size_t, size_t> inequality_constraints{}; /*!< inequality constraints */
   std::map<size_t, size_t> linear_constraints{};

   /* Hessian */
   int hessian_maximum_number_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */
   const bool fixed_hessian_sparsity{false};

   // purely virtual functions
   [[nodiscard]] virtual double evaluate_objective(const std::vector<double>& x) const = 0;
   virtual void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const = 0;
   [[nodiscard]] virtual double evaluate_constraint(int j, const std::vector<double>& x) const = 0;
   virtual void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const = 0;
   virtual void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         CSCSymmetricMatrix& hessian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         COOSymmetricMatrix& hessian) const = 0;

   virtual void set_initial_primal_point(std::vector<double>& x) = 0;
   virtual void set_initial_dual_point(std::vector<double>& multipliers) = 0;

   // auxiliary functions
   [[nodiscard]] std::vector<double> evaluate_constraints(const std::vector<double>& x) const;
   void determine_bounds_types(std::vector<Range>& variables_bounds, std::vector<ConstraintType>& status);
   void project_point_in_bounds(std::vector<double>& x) const;
   [[nodiscard]] double compute_constraint_violation(double constraint, size_t j) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set, Norm
      residual_norm) const;
   [[nodiscard]] bool is_constrained() const;

protected:
   void determine_constraints();
};

//template<size_t N, size_t M, size_t MAX_HESSIAN_NNZ>
//class CppProblem : public Problem {
//public:
//   CppProblem(std::string name, double (* objective)(const std::vector<double>& x), void (* objective_gradient)(const std::vector<double>& x,
//         std::vector<double>& gradient)): Problem(name, N, M, NONLINEAR), objective(objective), objective_gradient(objective_gradient),
//         current_variable(0), current_constraint(0) {
//   }
//
//   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override {
//      return this->objective(x);
//   }
//
//   void evaluate_objective_gradient(const std::vector<double>& /*x*/, SparseVector<double>& /*gradient*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   [[nodiscard]] double evaluate_constraint(int /*j*/, const std::vector<double>& /*x*/) const override {
//      return 0.;
//   }
//
//   void evaluate_constraints(const std::vector<double>& /*x*/, std::vector<double>& /*constraints*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   void constraint_gradient(const std::vector<double>& /*x*/, int /*j*/, SparseVector<double>& /*gradient*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   void constraint_jacobian(const std::vector<double>& /*x*/, std::vector<SparseVector<double>>& /*constraint_jacobian*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   void lagrangian_hessian(const std::vector<double>& /*x*/, double /*objective_multiplier*/, const std::vector<double>& /*multipliers*/,
//         CSCSymmetricMatrix& /*hessian*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   void lagrangian_hessian(const std::vector<double>& /*x*/, double /*objective_multiplier*/, const std::vector<double>& /*multipliers*/,
//         COOSymmetricMatrix& /*hessian*/) const override {
//      assert(false && "not yet implemented");
//   }
//
//   // initial point
//   void set_initial_primal_point(std::vector<double>& /*x*/) override {
//      assert(false && "not yet implemented");
//   }
//
//   void set_initial_dual_point(std::vector<double>& /*multipliers*/) override {
//      assert(false && "not yet implemented");
//   }
//
//private:
//   double (* objective)(const std::vector<double>& x);
//   void (* objective_gradient)(const std::vector<double>& x, std::vector<double>& gradient);
//   std::vector<double (*)(const std::vector<double>& x)> constraints;
//   int current_variable;
//   int current_constraint;
//};

// https://stackoverflow.com/questions/53252321/how-to-write-a-function-that-can-take-in-an-array-or-a-vector
//template<class Container>
//void f_container(const Container& container) {
//   std::for_each(std::begin(container), std::end(container), [](double xi) {
//      std::cout << xi << "\n";
//   });
//}

#endif // PROBLEM_H
