#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

template<typename LPSolverType>
class SLP : public Subproblem {
public:
   SLP(const Problem& problem, size_t number_variables, size_t number_constraints);

   void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   double compute_predicted_reduction(const Direction& direction, double step_length) const override;
   int get_hessian_evaluation_count() const override;

private:
   /* use pointers to allow polymorphism */
   const std::unique_ptr <LPSolverType> solver; /*!< Solver that solves the subproblem */
   std::vector<double> initial_point;
};

template<typename LPSolverType>
inline SLP<LPSolverType>::SLP(const Problem& problem, size_t number_variables, size_t number_constraints) : Subproblem(number_variables,
      number_constraints),
      solver(LPSolverFactory<LPSolverType>::create(number_variables, number_constraints)),
      initial_point(number_variables) {
   // register the original constraints bounds
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->constraints_bounds[j] = problem.constraint_bounds[j];
   }
}

template<typename LPSolverType>
inline void SLP<LPSolverType>::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
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

template<typename LPSolverType>
inline void SLP<LPSolverType>::update_objective_multiplier(const Problem& /*problem*/, const Iterate& current_iterate, double objective_multiplier) {
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

template<typename LPSolverType>
inline void SLP<LPSolverType>::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

template<typename LPSolverType>
inline Direction SLP<LPSolverType>::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* solve the LP */
   Direction direction = this->solver->solve_LP(variables_bounds, constraints_bounds, this->objective_gradient, this->constraints_jacobian,
         this->initial_point);
   // compute dual displacements (SQP methods compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   this->number_subproblems_solved++;
   return direction;
}

template<typename LPSolverType>
inline double SLP<LPSolverType>::compute_predicted_reduction(const Direction& direction, double step_length) const {
   // the predicted reduction is linear in the step length
   return -step_length * direction.objective;
}

template<typename LPSolverType>
inline int SLP<LPSolverType>::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}

#endif // SLP_H
