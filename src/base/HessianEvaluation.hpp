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
   explicit HessianEvaluation(int dimension, int hessian_maximum_number_nonzeros);
   virtual ~HessianEvaluation() = default;

   int dimension;
   CSCMatrix hessian;

   virtual void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) = 0;
   CSCMatrix modify_inertia(CSCMatrix& matrix, LinearSolver& linear_solver);
};

class ExactHessianEvaluation : public HessianEvaluation {
public:
   explicit ExactHessianEvaluation(int dimension, int hessian_maximum_number_nonzeros);
   ~ExactHessianEvaluation() override = default;

   void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) override;

protected:
   double objective_multiplier;
};

class ExactHessianInertiaControlEvaluation : public HessianEvaluation {
public:
   ExactHessianInertiaControlEvaluation(int dimension, int hessian_maximum_number_nonzeros, const std::string& linear_solve_name);
   ~ExactHessianInertiaControlEvaluation() override = default;

   void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) override;

protected:
   std::unique_ptr<LinearSolver> linear_solver_; /*!< Solver that solves the subproblem */
   double objective_multiplier;
};

class BFGSHessianEvaluation : public HessianEvaluation {
public:
   explicit BFGSHessianEvaluation(int dimension, int hessian_maximum_number_nonzeros);
   ~BFGSHessianEvaluation() override = default;

   void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) override;

private:
   CSCMatrix previous_hessian_;
   std::vector<double> previous_x_;
};

class HessianEvaluationFactory {
public:
   static std::unique_ptr<HessianEvaluation>
   create(const std::string& hessian_evaluation_method, int dimension, int hessian_maximum_number_nonzeros, bool convexify);
};

#endif // HESSIANEVALUATION_H
