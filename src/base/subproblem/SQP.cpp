#include <cmath>
#include <map>
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Vector.hpp"
#include "Logger.hpp"
#include "QPSolverFactory.hpp"

SQP::SQP(Problem& problem, std::string QP_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals) :
// maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      ActiveSetMethod(problem, scale_residuals), solver(
      QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints,
            problem.hessian_maximum_number_nonzeros + problem.number_variables, true)),
/* if no trust region is used, the problem should be convexified by controlling the inertia of the Hessian */
      hessian_evaluation(
            HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, problem.hessian_maximum_number_nonzeros,
                  !use_trust_region)) {
}

Direction SQP::compute_qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius) {
   DEBUG << "Current point: ";
   print_vector(DEBUG, current_iterate.x);
   DEBUG << "Current constraint multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.constraints);
   DEBUG << "Current lb multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
   DEBUG << "Current ub multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.upper_bounds);

   /* bounds of the variables */
   std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

   /* generate the initial point */
   std::vector<double> d0(variables_bounds.size()); // = {0.}

   /* solve the QP */
   Direction direction =
         solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian,
               this->hessian_evaluation->hessian, d0);
   this->number_subproblems_solved++;
   DEBUG << direction;
   return direction;
}

std::vector<Direction>
SQP::compute_directions(Problem& problem, Iterate& current_iterate, double /*objective_multiplier*/, double trust_region_radius) {
   /* compute optimality step */
   this->evaluate_optimality_iterate_(problem, current_iterate);
   Direction direction = this->compute_qp_step_(problem, *this->solver, current_iterate, trust_region_radius);

   if (direction.status != INFEASIBLE) {
      direction.phase = OPTIMALITY;
      direction.objective_multiplier = problem.objective_sign;
      direction.predicted_reduction = [&](Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) {
         return this->compute_qp_predicted_reduction_(problem, current_iterate, direction, step_length);
      };
      return std::vector<Direction>{direction};
   }
   else {
      /* infeasible subproblem during optimality phase */
      return this->restore_feasibility(problem, current_iterate, direction, trust_region_radius);
   }
}

Direction SQP::compute_l1qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, ConstraintPartition& constraint_partition,
      std::vector<double>& initial_solution, double trust_region_radius) {
   /* compute the objective */
   this->compute_l1_linear_objective_(current_iterate, constraint_partition);

   /* bounds of the variables */
   std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   std::vector<Range> constraints_bounds = this->generate_feasibility_bounds_(problem, current_iterate.constraints, constraint_partition);

   /* generate the initial point */
   std::vector<double>& d0 = initial_solution;

   /* solve the QP */
   Direction direction =
         solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian,
               this->hessian_evaluation->hessian, d0);
   direction.objective_multiplier = 0.;
   direction.constraint_partition = constraint_partition;
   this->number_subproblems_solved++;
   DEBUG << direction;
   return direction;
}

std::vector<Direction>
SQP::restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
   DEBUG << "\nCreating the restoration problem with " << phase_2_direction.constraint_partition.infeasible.size()
         << " infeasible constraints\n";
   this->evaluate_feasibility_iterate_(problem, current_iterate, phase_2_direction.constraint_partition);
   Direction direction =
         this->compute_l1qp_step_(problem, *this->solver, current_iterate, phase_2_direction.constraint_partition, phase_2_direction.x,
               trust_region_radius);
   direction.phase = RESTORATION;
   direction.objective_multiplier = problem.objective_sign;
   direction.predicted_reduction = [&](Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) {
      return this->compute_qp_predicted_reduction_(problem, current_iterate, direction, step_length);
   };
   return std::vector<Direction>{direction};
}

/* private methods */

void SQP::evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate) const {
   /* compute first- and second-order information */
   current_iterate.compute_objective_gradient(problem);
   current_iterate.compute_constraints_jacobian(problem);
   this->hessian_evaluation->compute(problem, current_iterate, problem.objective_sign, current_iterate.multipliers.constraints);
}

void SQP::evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition) {
   /* update the multipliers of the general constraints */
   std::vector<double>
         constraint_multipliers = this->generate_l1_multipliers_(problem, current_iterate.multipliers.constraints, constraint_partition);
   /* compute first- and second-order information */
   current_iterate.compute_constraints_jacobian(problem);
   //current_iterate.is_hessian_computed = false;
   double objective_multiplier = 0.;
   this->hessian_evaluation->compute(problem, current_iterate, objective_multiplier, constraint_multipliers);
}

double SQP::compute_qp_predicted_reduction_(Problem& /*problem*/, Iterate& current_iterate, Direction& direction, double step_length) {
   // the predicted reduction is quadratic in the step length
   if (step_length == 1.) {
      return -direction.objective;
   }
   else {
      double linear_term = dot(direction.x, current_iterate.objective_gradient);
      double quadratic_term = this->hessian_evaluation->hessian.quadratic_product(direction.x, direction.x) / 2.;
      return -step_length * (linear_term + step_length * quadratic_term);
   }
}
