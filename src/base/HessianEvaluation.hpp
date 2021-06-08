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
   explicit HessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros);
   virtual ~HessianEvaluation() = default;

   size_t dimension;
   CSCMatrix hessian;

   virtual void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) = 0;
   CSCMatrix modify_inertia(CSCMatrix& matrix, LinearSolver& linear_solver);
};

class ExactHessianEvaluation : public HessianEvaluation {
public:
   explicit ExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros);
   ~ExactHessianEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override;

protected:
   //double objective_multiplier_;
};

class ExactHessianInertiaControlEvaluation : public HessianEvaluation {
public:
   ExactHessianInertiaControlEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& linear_solve_name);
   ~ExactHessianInertiaControlEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override;

protected:
   std::unique_ptr<LinearSolver> linear_solver_; /*!< Solver that computes the inertia */
   //double objective_multiplier;
};

class BFGSHessianEvaluation : public HessianEvaluation {
public:
   explicit BFGSHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros);
   ~BFGSHessianEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override;

private:
   CSCMatrix previous_hessian_;
   std::vector<double> previous_x_;
};

class HessianEvaluationFactory {
public:
   static std::unique_ptr<HessianEvaluation>
   create(const std::string& hessian_evaluation_method, size_t dimension, size_t hessian_maximum_number_nonzeros, bool convexify);
};

#endif // HESSIANEVALUATION_H
