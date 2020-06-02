#ifndef PROBLEM_H
#define PROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "Constraint.hpp"
#include "Matrix.hpp"

enum FunctionType {
    LINEAR = 0, /*!< Linear function */
    QUADRATIC, /*!< Quadratic function */
    NONLINEAR /*!< Nonlinear function */
};

struct IEEE_Error : public std::exception {
    virtual const char* what() const throw () = 0;
};

struct IEEE_HEssianError : IEEE_Error {

    const char* what() const throw () {
        return "An IEEE error was encountered while evaluating the Hessian";
    }
};

struct IEEE_GradientError : IEEE_Error {

    const char* what() const throw () {
        return "An IEEE error was encountered while evaluating a gradient";
    }
};

struct IEEE_FunctionError : IEEE_Error {

    const char* what() const throw () {
        return "An IEEE error was encountered while evaluating a function";
    }
};

/*! \class Problem
 * \brief Optimization problem
 *
 *  Description of an optimization problem
 */
class Problem {
public:
    Problem(std::string name, int number_variables, int number_constraints);
    virtual ~Problem();

    std::string name;
    int number_variables; /*!< Number of variables */
    int number_constraints; /*!< Number of constraints */

    /* objective */
    double objective_sign; /*!< Sign of the objective function (1: minimization, -1: maximization) */
    std::string objective_name;
    FunctionType objective_type; /*!< Type of the objective (LINEAR, QUADRATIC, NONLINEAR) */
    virtual double objective(std::vector<double>& x) = 0;
    virtual std::vector<double> objective_dense_gradient(std::vector<double>& x) = 0;
    virtual std::map<int, double> objective_sparse_gradient(std::vector<double>& x) = 0;

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
    virtual double evaluate_constraint(int j, std::vector<double>& x) = 0;
    virtual std::vector<double> evaluate_constraints(std::vector<double>& x) = 0;
    virtual std::vector<double> constraint_dense_gradient(int j, std::vector<double>& x) = 0;
    virtual void constraint_sparse_gradient(std::vector<double>& x, int j, std::map<int, double>& gradient) = 0;
    virtual std::vector<std::map<int, double> > constraints_sparse_jacobian(std::vector<double>& x) = 0;
    void determine_bounds_types(std::vector<Range>& variables_bounds, std::vector<ConstraintType>& status);
    std::map<int, int> equality_constraints; /*!< inequality constraints */
    std::map<int, int> inequality_constraints; /*!< inequality constraints */
    std::map<int, int> linear_constraints;
    
    /* Hessian */
    int hessian_maximum_number_nonzeros; /*!< Number of nonzero elements in the Hessian */
    virtual CSCMatrix lagrangian_hessian(std::vector<double>& x, double objective_multiplier, std::vector<double>& multipliers) = 0;

    virtual std::vector<double> primal_initial_solution() = 0;
    virtual std::vector<double> dual_initial_solution() = 0;
    
    double compute_constraint_residual(std::vector<double>& constraints, std::string norm_value);
    double compute_constraint_residual(std::vector<double>& constraints, std::set<int> constraint_set, std::string norm_value);

    int number_eval_objective;
    int number_eval_constraints;
    int number_eval_jacobian;
    int number_eval_hessian;

protected:
    void determine_constraints_();
};

class CppProblem: public Problem {
public:
    CppProblem(std::string name, int number_variables, int number_constraints, double (*objective)(std::vector<double> x), std::vector<double> (*objective_gradient)(std::vector<double> x));
    
    double objective(std::vector<double>& x) override;
    std::vector<double> objective_dense_gradient(std::vector<double>& x) override;
    std::map<int, double> objective_sparse_gradient(std::vector<double>& x) override;
    
    double evaluate_constraint(int j, std::vector<double>& x) override;
    std::vector<double> evaluate_constraints(std::vector<double>& x) override;
    std::vector<double> constraint_dense_gradient(int j, std::vector<double>& x) override;
    void constraint_sparse_gradient(std::vector<double>& x, int j, std::map<int, double>& gradient) override;
    std::vector<std::map<int, double> > constraints_sparse_jacobian(std::vector<double>& x) override;
    
    CSCMatrix lagrangian_hessian(std::vector<double>& x, double objective_multiplier, std::vector<double>& multipliers) override;
    
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
    std::vector<double (*)(std::vector<double> x)> constraints_;
    int current_variable_;
    int current_constraint_;
};

#endif // PROBLEM_H
