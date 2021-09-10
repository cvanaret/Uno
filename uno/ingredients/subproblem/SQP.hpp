#ifndef SQP_H
#define SQP_H

#include <cassert>
#include "Subproblem.hpp"
#include "HessianEvaluation.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

template<typename QPSolverType>
class SQP : public Subproblem {
public:
   SQP(const Problem& problem, size_t number_variables, size_t number_constraints, const std::string& hessian_evaluation_method,
         bool use_trust_region);

   void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   int get_hessian_evaluation_count() const override;

protected:
   /* use pointers to allow polymorphism */
   const std::unique_ptr <QPSolver<typename QPSolverType::matrix_type>> solver; /*!< Solver that solves the subproblem */
   const std::unique_ptr <HessianEvaluation<typename QPSolverType::matrix_type>> hessian_evaluation; /*!< Strategy to compute or approximate the
 * Hessian */
   std::vector<double> initial_point;
};

template<typename QPSolverType>
inline SQP<QPSolverType>::SQP(const Problem& problem, size_t number_variables, size_t number_constraints, const std::string&
hessian_evaluation_method, bool use_trust_region) :
      Subproblem(number_variables, number_constraints),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory<QPSolverType>::create(number_variables, number_constraints,
            problem.hessian_maximum_number_nonzeros + number_variables, true)),
      /* if no trust region is used, the problem should be convexified by controlling the inertia of the Hessian */
      hessian_evaluation(HessianEvaluationFactory<typename QPSolverType::matrix_type>::create(hessian_evaluation_method, number_variables,
            problem.hessian_maximum_number_nonzeros + problem.number_variables, !use_trust_region)),
      initial_point(number_variables) {
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   for (auto& row: this->constraints_jacobian) {
      row.clear();
   }
   problem.constraints_jacobian(current_iterate.x, this->constraints_jacobian);
   this->objective_gradient.clear();
   problem.evaluate_objective_gradient(current_iterate.x, this->objective_gradient);
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

   /* set the initial point */
   clear(this->initial_point);
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_evaluation->compute(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // scale objective gradient
   if (objective_multiplier == 0.) {
      this->objective_gradient.clear();
   }
   else if (objective_multiplier < 1.) {
      this->objective_gradient = current_iterate.objective_gradient;
      scale(this->objective_gradient, objective_multiplier);
   }
   clear(this->initial_point);
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

template<typename QPSolverType>
inline Direction SQP<QPSolverType>::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* compute QP direction */
   Direction direction = this->solver->solve_QP(this->variables_bounds, this->constraints_bounds, this->objective_gradient,
         this->constraints_jacobian, this->hessian_evaluation->hessian, this->initial_point);
   this->number_subproblems_solved++;

   // compute dual displacements (SQP methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   return direction;
}

template<typename QPSolverType>
inline PredictedReductionModel SQP<QPSolverType>::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture this and direction by reference
      // precompute expensive quantities
      double linear_term = dot(direction.x, this->objective_gradient);
      double quadratic_term = this->hessian_evaluation->hessian.quadratic_product(direction.x, direction.x) / 2.;
      std::cout << "EXPENSIVE STUFF\n";
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * (linear_term + step_length * quadratic_term);
      };
   });
}

template<typename QPSolverType>
inline int SQP<QPSolverType>::get_hessian_evaluation_count() const {
   return this->hessian_evaluation->evaluation_count;
}

#endif // SQP_H
