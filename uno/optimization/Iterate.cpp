// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Model.hpp"
#include "tools/Logger.hpp"

size_t Iterate::number_eval_objective = 0;
size_t Iterate::number_eval_constraints = 0;
size_t Iterate::number_eval_objective_gradient = 0;
size_t Iterate::number_eval_jacobian = 0;

Iterate::Iterate(size_t max_number_variables, size_t max_number_constraints) :
      number_variables(max_number_variables), number_constraints(max_number_constraints),
      primals(max_number_variables), multipliers(max_number_variables, max_number_constraints),
      evaluations(max_number_variables, max_number_constraints),
      lagrangian_gradient(max_number_variables) {
}

void Iterate::evaluate_objective(const Model& model) {
   if (not this->is_objective_computed) {
      // evaluate the objective
      this->evaluations.objective = model.evaluate_objective(this->primals);
      // check finiteness
      if (this->evaluations.objective == INF<double>) {
         throw FunctionEvaluationError();
      }
      this->is_objective_computed = true;
      Iterate::number_eval_objective++;
   }
}

void Iterate::evaluate_constraints(const Model& model) {
   if (not this->are_constraints_computed) {
      // evaluate the constraints
      model.evaluate_constraints(this->primals, this->evaluations.constraints);
      // check finiteness
      if (std::any_of(this->evaluations.constraints.cbegin(), this->evaluations.constraints.cend(), [](double constraint_j) {
         return constraint_j == INF<double>;
      })) {
         throw FunctionEvaluationError();
      }
      this->are_constraints_computed = true;
      Iterate::number_eval_constraints++;
   }
}

void Iterate::evaluate_objective_gradient(const Model& model) {
   if (not this->is_objective_gradient_computed) {
      this->evaluations.objective_gradient.clear();
      // evaluate the objective gradient
      model.evaluate_objective_gradient(this->primals, this->evaluations.objective_gradient);
      this->is_objective_gradient_computed = true;
      Iterate::number_eval_objective_gradient++;
   }
}

void Iterate::evaluate_constraint_jacobian(const Model& model) {
   if (not this->is_constraint_jacobian_computed) {
      for (auto& row: this->evaluations.constraint_jacobian) {
         row.clear();
      }
      // evaluate the constraint Jacobian
      model.evaluate_constraint_jacobian(this->primals, this->evaluations.constraint_jacobian);
      this->is_constraint_jacobian_computed = true;
      Iterate::number_eval_jacobian++;
   }
}

void Iterate::set_number_variables(size_t new_number_variables) {
   this->number_variables = new_number_variables;
   this->primals.resize(new_number_variables);
   this->multipliers.lower_bounds.resize(new_number_variables);
   this->multipliers.upper_bounds.resize(new_number_variables);
   this->evaluations.objective_gradient.reserve(new_number_variables);
   this->lagrangian_gradient.resize(new_number_variables);
}

std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
   stream << "Primal variables: "; print_vector(stream, iterate.primals);
   stream << "            ┌ Constraint: "; print_vector(stream, iterate.multipliers.constraints);
   stream << "Multipliers │ Lower bound: "; print_vector(stream, iterate.multipliers.lower_bounds);
   stream << "            └ Upper bound: "; print_vector(stream, iterate.multipliers.upper_bounds);

   stream << "Objective value: " << iterate.evaluations.objective << '\n';

   stream << "          ┌ Optimality stationarity: " << iterate.residuals.optimality_stationarity << '\n';
   stream << "          │ Feasibility stationarity: " << iterate.residuals.feasibility_stationarity << '\n';
   stream << "Residuals │ Constraint violation: " << iterate.residuals.infeasibility << '\n';
   stream << "          │ Optimality complementarity: " << iterate.residuals.optimality_complementarity << '\n';
   stream << "          └ Feasibility complementarity: " << iterate.residuals.feasibility_complementarity << '\n';

   stream << "                  ┌ Infeasibility: " << iterate.progress.infeasibility << '\n';
   stream << "Progress measures │ Optimality: " << iterate.progress.optimality(1.) << '\n';
   stream << "                  └ Auxiliary terms: " << iterate.progress.auxiliary_terms << '\n';

   stream << "Lagrangian gradient: " << iterate.lagrangian_gradient;
   return stream;
}