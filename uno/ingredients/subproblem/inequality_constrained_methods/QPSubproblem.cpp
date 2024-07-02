// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPSubproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/HessianModelFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

QPSubproblem::QPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
      size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) :
      InequalityConstrainedMethod(number_variables, number_constraints),
      use_regularization(options.get_string("globalization_mechanism") != "TR" || options.get_bool("convexify_QP")),
      enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
      // if no trust region is used, the problem should be convexified to guarantee boundedness
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), number_variables,
            number_hessian_nonzeros + number_variables, this->use_regularization, options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), number_variables, number_constraints,
            number_objective_gradient_nonzeros, number_jacobian_nonzeros,
            // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
            std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, hessian_model->hessian->capacity),
            options)) {
}

void QPSubproblem::initialize_statistics(Statistics& statistics, const Options& options) {
   if (this->use_regularization) {
      statistics.add_column("regularization", Statistics::double_width, options.get_int("statistics_regularization_column_order"));
   }
}

bool QPSubproblem::generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) {
   if (this->enforce_linear_constraints_at_initial_iterate) {
      const bool is_feasible = Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers,
            *this->solver);
      if (not is_feasible) {
         return false;
      }
   }
   return true;
}

void QPSubproblem::evaluate_functions(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
      const Multipliers& current_multipliers, const WarmstartInformation& warmstart_information) {
   // Lagrangian Hessian
   if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
      this->hessian_model->evaluate(statistics, problem, current_iterate.primals, current_multipliers.constraints);
   }
   // objective gradient, constraints and constraint Jacobian
   if (warmstart_information.objective_changed) {
      problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);
   }
   if (warmstart_information.constraints_changed) {
      problem.evaluate_constraints(current_iterate, this->constraints);
      problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
   }
}

void QPSubproblem::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
      Direction& direction, const WarmstartInformation& warmstart_information) {
   // evaluate the functions at the current iterate
   this->evaluate_functions(statistics, problem, current_iterate, current_multipliers, warmstart_information);

   // set bounds of the variable displacements
   if (warmstart_information.variable_bounds_changed) {
      this->set_direction_bounds(problem, current_iterate);
   }

   // set bounds of the linearized constraints
   if (warmstart_information.constraint_bounds_changed) {
      this->set_linearized_constraint_bounds(problem, this->constraints);
   }

   // solve the QP
   this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->direction_lower_bounds, this->direction_upper_bounds,
         this->linearized_constraints_lower_bounds, this->linearized_constraints_upper_bounds, this->objective_gradient,
         this->constraint_jacobian, *this->hessian_model->hessian, this->initial_point, direction, warmstart_information);
   InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
   this->number_subproblems_solved++;
   // reset the initial point
   this->initial_point.fill(0.);
}

const SymmetricMatrix<double>& QPSubproblem::get_lagrangian_hessian() const {
   return *this->hessian_model->hessian;
}

size_t QPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}
