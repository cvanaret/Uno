#ifndef PROBLEM_H
#define PROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "Constraint.hpp"
#include "Matrix.hpp"
#include "SparseGradient.hpp"

enum FunctionType {
    LINEAR = 0, /*!< Linear function */
    QUADRATIC, /*!< Quadratic function */
    NONLINEAR /*!< Nonlinear function */
};

struct NumericalError : public std::exception {
    virtual const char* what() const throw () = 0;
};

struct HessianNumericalError : NumericalError {

    const char* what() const throw () {
        return "A numerical error was encountered while evaluating the Hessian";
    }
};

struct GradientNumericalError : NumericalError {

    const char* what() const throw () {
        return "A numerical error was encountered while evaluating a gradient";
    }
};

struct FunctionNumericalError : NumericalError {

    const char* what() const throw () {
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
    Problem(std::string& name, int number_variables, int number_constraints, FunctionType problem_type);
    virtual ~Problem() = default;

    static std::map<FunctionType, std::string> type_to_string;

    std::string name;
    int number_variables; /*!< Number of variables */
    int number_constraints; /*!< Number of constraints */
    FunctionType type;

    /* objective */
    double objective_sign; /*!< Sign of the objective function (1: minimization, -1: maximization) */
    std::string objective_name;
    FunctionType objective_type; /*!< Type of the objective (LINEAR, QUADRATIC, NONLINEAR) */
    [[nodiscard]] virtual double objective(const std::vector<double>& x) const = 0;
    virtual SparseGradient objective_gradient(std::vector<double>& x) const = 0;

    /* variables */
    std::vector<std::string> variable_name;
    //std::vector<bool> variable_discrete;
    std::vector<Range> variables_bounds;
    std::vector<ConstraintType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */

    /* constraints */
    std::vector<std::string> constraint_name;
    std::vector<Range> constraint_bounds;
    std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
    std::vector<ConstraintType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED) */
    virtual double evaluate_constraint(int j, std::vector<double>& x) const = 0;
    virtual std::vector<double> evaluate_constraints(std::vector<double>& x) const = 0;
    virtual void constraint_gradient(std::vector<double>& x, int j, SparseGradient& gradient) const = 0;
    virtual std::vector<SparseGradient> constraints_jacobian(std::vector<double>& x) const = 0;
    void determine_bounds_types(std::vector<Range>& variables_bounds, std::vector<ConstraintType>& status);
    std::map<int, int> equality_constraints; /*!< inequality constraints */
    std::map<int, int> inequality_constraints; /*!< inequality constraints */
    std::map<int, int> linear_constraints;

    /* Hessian */
    int hessian_maximum_number_nonzeros; /*!< Number of nonzero elements in the Hessian */
    [[nodiscard]] virtual CSCMatrix lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers) const = 0;

    virtual std::vector<double> primal_initial_solution() = 0;
    virtual std::vector<double> dual_initial_solution() = 0;

    [[nodiscard]] double compute_constraint_residual(const std::vector<double>& constraints, const std::string& norm_value) const;
    [[nodiscard]] double compute_constraint_residual(const std::vector<double>& constraints, const std::set<int>& constraint_set, const std::string& norm_value) const;

protected:
    void determine_constraints_();
};

class CppProblem : public Problem {
public:
    CppProblem(std::string name, int number_variables, int number_constraints, double (*objective)(std::vector<double> x), std::vector<double> (*objective_gradient)(std::vector<double> x));

    double objective(const std::vector<double>& x) const override;
    SparseGradient objective_gradient(std::vector<double>& x) const override;

    double evaluate_constraint(int j, std::vector<double>& x) const override;
    std::vector<double> evaluate_constraints(std::vector<double>& x) const override;
    void constraint_gradient(std::vector<double>& x, int j, SparseGradient& gradient) const override;
    std::vector<SparseGradient> constraints_jacobian(std::vector<double>& x) const override;

    CSCMatrix lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers) const override;

    /* variables */
    void add_variable(std::string name, Range& bounds);
    void add_variables(std::vector<std::string> names, std::vector<Range>& bounds);

    /* constraints */
    void add_constraint(std::string name, Range& bounds, FunctionType type);

    /* initial point */
    std::vector<double> primal_initial_solution() override;
    std::vector<double> dual_initial_solution() override;

private:
    double (*objective_)(std::vector<double> x);
    std::vector<double> (*objective_gradient_)(std::vector<double> x);
    std::vector<double (*)(std::vector<double> x) > constraints_;
    int current_variable_;
    int current_constraint_;
};

#endif // PROBLEM_H
