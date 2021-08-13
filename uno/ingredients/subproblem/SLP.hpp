#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

template<typename LPSolverType>
class SLP : public Subproblem {
public:
   SLP(size_t number_variables, size_t number_constraints) : Subproblem(number_variables, number_constraints),
         solver(LPSolverFactory<LPSolverType>::create(number_variables, number_constraints)),
         initial_point(number_variables) {
   }

   void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override {
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

   void update_objective_multiplier(const Problem& /*problem*/, const Iterate& current_iterate, double objective_multiplier) override {
      // scale objective gradient
      if (objective_multiplier == 0.) {
         clear(this->objective_gradient);
      }
      else if (objective_multiplier < 1.) {
         this->objective_gradient = current_iterate.objective_gradient;
         scale(this->objective_gradient, objective_multiplier);
      }
      clear(this->initial_point);
   }

   void set_initial_point(const std::vector<double>& point) override {
      copy_from(this->initial_point, point);
   }

   Direction solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) override {
      /* solve the LP */
      Direction direction = this->solver->solve_LP(variables_bounds, constraints_bounds, this->objective_gradient,
            this->constraints_jacobian,this->initial_point);
      // compute dual displacements (SQP methods compute the new duals, not the displacements)
      for (size_t j = 0; j < problem.number_constraints; j++) {
         direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
      }
      this->number_subproblems_solved++;
      return direction;
   }

   double compute_predicted_reduction(const Direction& direction, double step_length) const override {
      // the predicted reduction is linear in the step length
      return -step_length * direction.objective;
   }

   int get_hessian_evaluation_count() const override {
      // no second order evaluation is used
      return 0;
   }

private:
   /* use references to allow polymorphism */
   const std::unique_ptr<LPSolverType> solver; /*!< Solver that solves the subproblem */
   std::vector<double> initial_point;
};

#endif // SLP_H
