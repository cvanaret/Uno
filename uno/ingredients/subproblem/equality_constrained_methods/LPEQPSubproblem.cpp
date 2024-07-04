// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPEQPSubproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/HessianModelFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"
#include "tools/Statistics.hpp"

LPEQPSubproblem::LPEQPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
      size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) :
      InequalityConstrainedMethod(number_variables, number_constraints),
      use_regularization(true),
      enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
      // if no trust region is used, the problem should be convexified to guarantee boundedness
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), number_variables,
            number_hessian_nonzeros + number_variables, this->use_regularization, options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), number_variables, number_constraints,
            number_objective_gradient_nonzeros, number_jacobian_nonzeros,
            // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
            std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, this->hessian_model->hessian->capacity),
            options)) {
}

void LPEQPSubproblem::initialize_statistics(Statistics& statistics, const Options& options) {
   if (this->use_regularization) {
      statistics.add_column("regularization", Statistics::double_width, options.get_int("statistics_regularization_column_order"));
   }
}

bool LPEQPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*first_iterate*/) {
   return true;
}

void LPEQPSubproblem::evaluate_functions(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
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

void LPEQPSubproblem::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
      const Multipliers& current_multipliers, Direction& direction, const WarmstartInformation& warmstart_information) {
    
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

   // solve LP subproblem
   Direction LP_direction(problem.number_variables, problem.number_constraints);
   this->solve_LP(problem, current_multipliers, LP_direction, warmstart_information);

   DEBUG << "d^*(LP) = "; print_vector(DEBUG, view(LP_direction.primals, 0, LP_direction.number_variables));
   DEBUG << "bound constraints active at lower bound =";
   for (size_t i: LP_direction.active_set.bounds.at_lower_bound) {
      DEBUG << " x" << i;
   }
   DEBUG << '\n';
   DEBUG << "bound constraints active at upper bound =";
   for (size_t i: LP_direction.active_set.bounds.at_upper_bound) {
      DEBUG << " x" << i;
   }
   DEBUG << '\n';
   DEBUG << "constraints at lower bound =";
   for (size_t j: LP_direction.active_set.constraints.at_lower_bound) {
      DEBUG << " c" << j;
   }
   DEBUG << '\n';
   DEBUG << "constraints at upper bound =";
   for (size_t j: LP_direction.active_set.constraints.at_upper_bound) {
      DEBUG << " c" << j;
   }
   DEBUG << '\n';

   // set up EQP subproblem: set inactive bounds +/- INF and others as equations
   this->set_variable_EQP_bounds(problem, current_iterate, LP_direction);
   this->set_linearized_EQP_bounds(problem, this->constraints, LP_direction);

   // compute EQP direction
   this->solve_QP(problem, current_multipliers, direction, warmstart_information);
}

void LPEQPSubproblem::solve_LP(const OptimizationProblem& problem, const Multipliers& current_multipliers, Direction& direction_LP,
      const WarmstartInformation& warmstart_information) {
   this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->direction_lower_bounds, this->direction_upper_bounds,
         this->linearized_constraints_lower_bounds, this->linearized_constraints_upper_bounds, this->objective_gradient,
         this->constraint_jacobian, this->initial_point, direction_LP, warmstart_information);
   InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction_LP.multipliers);
   this->number_subproblems_solved++;
}

void LPEQPSubproblem::solve_QP(const OptimizationProblem& problem, const Multipliers& current_multipliers, Direction& direction,
      const WarmstartInformation& warmstart_information) {
   this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->direction_lower_bounds,
         this->direction_upper_bounds, this->linearized_constraints_lower_bounds, this->linearized_constraints_upper_bounds, this->objective_gradient,
         this->constraint_jacobian, *this->hessian_model->hessian, this->initial_point, direction, warmstart_information);
   InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
   this->number_subproblems_solved++;
}

void LPEQPSubproblem::set_variable_EQP_bounds(const OptimizationProblem& problem, const Iterate& current_iterate, Direction& LP_direction) {
   // initialize all bounds to {-INF<double>,+INF<double>}
   for (size_t i: Range(problem.number_variables)) {
      this->direction_lower_bounds[i] = -INF<double>;
      this->direction_upper_bounds[i] = INF<double>;
   }
   // set active lower bounds as equality constraints
   for (size_t i: LP_direction.active_set.bounds.at_lower_bound) {
      const double lb = problem.variable_lower_bound(i) - current_iterate.primals[i];
      if (lb >= -(this->trust_region_radius)) {
         this->direction_lower_bounds[i] = this->direction_upper_bounds[i] = lb;
      }
   }
   // set active lower bounds as equality constraints
   for (size_t i: LP_direction.active_set.bounds.at_upper_bound) {
      const double ub = problem.variable_upper_bound(i) - current_iterate.primals[i];
      if (ub <= this->trust_region_radius) {
         this->direction_lower_bounds[i] = this->direction_upper_bounds[i] = ub;
      }
   }
}

void LPEQPSubproblem::set_linearized_EQP_bounds(const OptimizationProblem& problem, const std::vector<double>& current_constraints, Direction& LP_direction) {
   // initialize all constraint bounds to {-INF<double>,+INF<double>}
   for (size_t j: Range(problem.number_constraints)) {
      this->linearized_constraints_lower_bounds[j] = -INF<double>;
      this->linearized_constraints_upper_bounds[j] = INF<double>;
   }
   // set active lower constraint bounds as equality constraints
   for (size_t j: LP_direction.active_set.constraints.at_lower_bound) {
      const double lb = problem.constraint_lower_bound(j) - current_constraints[j];
      this->linearized_constraints_lower_bounds[j] = this->linearized_constraints_upper_bounds[j] = lb;
   }
   // set active upper constraint bounds as equality constraints
   for (size_t j: LP_direction.active_set.constraints.at_upper_bound) {
      const double ub = problem.constraint_upper_bound(j) - current_constraints[j];
      this->linearized_constraints_lower_bounds[j] = this->linearized_constraints_upper_bounds[j] = ub;
   }
}

const SymmetricMatrix<double>& LPEQPSubproblem::get_lagrangian_hessian() const {
   return *this->hessian_model->hessian;
}

size_t LPEQPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}