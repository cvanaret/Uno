#ifndef UNO_MODEL_H
#define UNO_MODEL_H

#include <string>
#include <vector>
#include <map>
#include "linear_algebra/CSCSymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

struct Range {
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

bool is_finite(double value);

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
   static void determine_bounds_types(std::vector<Range>& variables_bounds, std::vector<ConstraintType>& status);
   void project_point_in_bounds(std::vector<double>& x) const;
   [[nodiscard]] bool is_constrained() const;

protected:
   size_t hessian_maximum_number_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */
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

#endif // UNO_MODEL_H
