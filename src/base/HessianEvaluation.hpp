#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <ostream>
#include <map>
#include <memory>
#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Matrix.hpp"
#include "LinearSolver.hpp"

// virtual (abstract) class
class HessianEvaluation {
public:
    HessianEvaluation(int dimension);
    virtual ~HessianEvaluation();

    int dimension;
    
    virtual void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) = 0;
    CSCMatrix modify_inertia(CSCMatrix& hessian, std::shared_ptr<LinearSolver> linear_solver);
};

class ExactHessianEvaluation : public HessianEvaluation {
public:
    ExactHessianEvaluation(int dimension);

    void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers);
};

class ExactHessianInertiaControlEvaluation : public HessianEvaluation {
public:
    ExactHessianInertiaControlEvaluation(int dimension, std::string linear_solve_name);

    void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers);
private:
    std::shared_ptr<LinearSolver> linear_solver_; /*!< Solver that solves the subproblem */
};

class BFGSHessianEvaluation : public HessianEvaluation {
public:
    BFGSHessianEvaluation(int dimension);

    void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers);
    
private:
    CSCMatrix previous_hessian_;
    std::vector<double> previous_x_;
};

class HessianEvaluationFactory {
    public:
		static std::shared_ptr<HessianEvaluation> create(std::string hessian_evaluation_method, int dimension, bool convexify);
};
#endif // HESSIANEVALUATION_H
