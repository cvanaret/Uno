// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSubproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/LP/LPSolverFactory.hpp"
#include "tools/Options.hpp"

LPSubproblem::LPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
      size_t number_jacobian_nonzeros, const Options& options) :
      InequalityConstrainedMethod(number_variables, number_constraints),
      solver(LPSolverFactory::create(options.get_string("LP_solver"), number_variables, number_constraints,
            number_objective_gradient_nonzeros, number_jacobian_nonzeros, options)),
      zero_hessian(COOSymmetricMatrix<double>::zero(number_variables)) {
}

bool LPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   return true;
}

void LPSubproblem::evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate, const WarmstartInformation& warmstart_information) {
   // objective gradient
   if (warmstart_information.objective_changed) {
      problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);
   }
   // constraints and constraint Jacobian
   if (warmstart_information.constraints_changed) {
      problem.evaluate_constraints(current_iterate, this->constraints);
      problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
   }
}

void LPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
      Direction& direction, const WarmstartInformation& warmstart_information) {
   // evaluate the functions at the current iterate
   this->evaluate_functions(problem, current_iterate, warmstart_information);

   // set bounds of the variable displacements
   if (warmstart_information.variable_bounds_changed) {
      this->set_direction_bounds(problem, current_iterate);
   }

   // set bounds of the linearized constraints
   if (warmstart_information.constraint_bounds_changed) {
      this->set_linearized_constraint_bounds(problem, this->constraints);
   }

   // solve the LP
   this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->direction_lower_bounds, this->direction_upper_bounds,
         this->linearized_constraints_lower_bounds, this->linearized_constraints_upper_bounds, this->objective_gradient,
         this->constraint_jacobian, this->initial_point, direction, warmstart_information);
   InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
   this->number_subproblems_solved++;
   // reset the initial point
   this->initial_point.fill(0.);
}

const SymmetricMatrix<double>& LPSubproblem::get_lagrangian_hessian() const {
   return this->zero_hessian;
}

size_t LPSubproblem::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}
