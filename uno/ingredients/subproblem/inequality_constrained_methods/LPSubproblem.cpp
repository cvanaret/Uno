// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <linear_algebra/COOSymmetricMatrix.hpp>
#include "LPSubproblem.hpp"
#include "solvers/LP/LPSolverFactory.hpp"

LPSubproblem::LPSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_objective_gradient_nonzeros,
      size_t max_number_jacobian_nonzeros, const Options& options) :
      InequalityConstrainedMethod(max_number_variables, max_number_constraints),
      solver(LPSolverFactory::create(options.get_string("LP_solver"), max_number_variables, max_number_constraints,
            max_number_objective_gradient_nonzeros, max_number_jacobian_nonzeros, options)),
      zero_hessian(COOSymmetricMatrix<double>::zero(max_number_variables)) {
}

bool LPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   return true;
}

void LPSubproblem::evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate, const WarmstartInformation& warmstart_information) {
   // objective gradient
   if (warmstart_information.objective_changed) {
      problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   }
   // constraints and constraint Jacobian
   if (warmstart_information.constraints_changed) {
      problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
      problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
   }
}

void LPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate, Direction& direction,
      const WarmstartInformation& warmstart_information) {
   // evaluate the functions at the current iterate
   this->evaluate_functions(problem, current_iterate, warmstart_information);

   // set bounds of the variable displacements
   if (warmstart_information.variable_bounds_changed) {
      this->set_direction_bounds(problem, current_iterate);
   }

   // set bounds of the linearized constraints
   if (warmstart_information.constraint_bounds_changed) {
      this->set_linearized_constraint_bounds(problem, this->evaluations.constraints);
   }

   // solve the LP
   this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->direction_lower_bounds, this->direction_upper_bounds,
         this->linearized_constraints_lower_bounds, this->linearized_constraints_upper_bounds, this->evaluations.objective_gradient,
         this->evaluations.constraint_jacobian, this->initial_point, direction, warmstart_information);
   InequalityConstrainedMethod::compute_dual_displacements(problem, current_iterate, direction);
   this->number_subproblems_solved++;
   // reset the initial point
   initialize_vector(this->initial_point, 0.);
}

const SymmetricMatrix<double>& LPSubproblem::get_lagrangian_hessian() const {
   return this->zero_hessian;
}

size_t LPSubproblem::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}
