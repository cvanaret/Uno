// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Model.hpp"
#include "tools/Logger.hpp"

size_t Iterate::number_eval_objective = 0;
size_t Iterate::number_eval_constraints = 0;
size_t Iterate::number_eval_jacobian = 0;

Iterate::Iterate(size_t max_number_variables, size_t max_number_constraints) :
      number_variables(max_number_variables), number_constraints(max_number_constraints),
      primals(max_number_variables), multipliers(max_number_variables, max_number_constraints),
      model_evaluations(max_number_variables, max_number_constraints),
      reformulation_evaluations(max_number_variables, max_number_constraints),
      lagrangian_gradient(max_number_variables) {
}

void Iterate::evaluate_objective(const Model& model) {
   if (!this->is_objective_computed) {
      // evaluate the objective
      this->model_evaluations.objective = model.evaluate_objective(this->primals);
      // check finiteness
      if (this->model_evaluations.objective == INF<double>) {
         throw FunctionEvaluationError();
      }
      this->is_objective_computed = true;
      Iterate::number_eval_objective++;
   }
}

void Iterate::evaluate_constraints(const Model& model) {
   if (!this->are_constraints_computed) {
      // evaluate the constraints
      model.evaluate_constraints(this->primals, this->model_evaluations.constraints);
      this->are_constraints_computed = true;
      Iterate::number_eval_constraints++;
   }
}

void Iterate::evaluate_objective_gradient(const Model& model) {
   if (!this->is_objective_gradient_computed) {
      this->model_evaluations.objective_gradient.clear();
      // evaluate the objective gradient
      model.evaluate_objective_gradient(this->primals, this->model_evaluations.objective_gradient);
      this->is_objective_gradient_computed = true;
   }
}

void Iterate::evaluate_constraint_jacobian(const Model& model) {
   if (!this->is_constraint_jacobian_computed) {
      for (auto& row: this->model_evaluations.constraint_jacobian) {
         row.clear();
      }
      // evaluate the constraint Jacobian
      model.evaluate_constraint_jacobian(this->primals, this->model_evaluations.constraint_jacobian);
      this->is_constraint_jacobian_computed = true;
      Iterate::number_eval_jacobian++;
   }
}

void Iterate::evaluate_lagrangian_gradient(const Model& model, double objective_multiplier, const std::vector<double>& constraint_multipliers,
      const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers) {
   /*
   std::cout << "Primals: "; print_vector(std::cout, this->primals);
   std::cout << "Duals: "; print_vector(std::cout, this->multipliers.constraints);
   std::cout << "Duals lb: "; print_vector(std::cout, this->multipliers.lower_bounds);
   std::cout << "Duals ub: "; print_vector(std::cout, this->multipliers.upper_bounds);
    */
   initialize_vector(this->lagrangian_gradient, 0.);

   // objective gradient: scale it with the objective multiplier
   if (objective_multiplier != 0.) {
      this->evaluate_objective_gradient(model);
      this->model_evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
         this->lagrangian_gradient[i] += objective_multiplier * derivative;
      });
   }

   // constraints
   this->evaluate_constraint_jacobian(model);
   for (size_t j = 0; j < model.number_constraints; j++) {
      // std::cout << "Grad c" << j << ": " << this->original_evaluations.constraint_jacobian[j] << "\n";
      const double multiplier_j = constraint_multipliers[j];
      if (multiplier_j != 0.) {
         this->model_evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            this->lagrangian_gradient[i] -= multiplier_j * derivative;
         });
      }
   }

   // std::cout << "Iterate has " << this->number_variables << " variables\n";
   // bound constraints
   for (size_t i = 0; i < this->number_variables; i++) {
      this->lagrangian_gradient[i] -= lower_bounds_multipliers[i] + upper_bounds_multipliers[i];
   }
}

void Iterate::set_number_variables(size_t new_number_variables) {
   this->number_variables = new_number_variables;
   this->primals.resize(new_number_variables);
   this->multipliers.lower_bounds.resize(new_number_variables);
   this->multipliers.upper_bounds.resize(new_number_variables);
   this->model_evaluations.objective_gradient.reserve(new_number_variables);
   this->lagrangian_gradient.resize(new_number_variables);
}

void Iterate::reset_evaluations() {
   this->is_objective_computed = false;
   this->is_objective_gradient_computed = false;
   this->are_constraints_computed = false;
   this->is_constraint_jacobian_computed = false;
}

std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
   stream << "Primal variables: ";
   print_vector(stream, iterate.primals);
   stream << "Lower bound multipliers: ";
   print_vector(stream, iterate.multipliers.lower_bounds);
   stream << "Upper bound multipliers: ";
   print_vector(stream, iterate.multipliers.upper_bounds);
   stream << "Constraint multipliers: ";
   print_vector(stream, iterate.multipliers.constraints);

   stream << "Objective value: " << iterate.model_evaluations.objective << '\n';
   stream << "Constraint violation: " << iterate.primal_constraint_violation << '\n';
   stream << "Stationarity (KKT/FJ) error: " << iterate.stationarity_error << '\n';
   stream << "Complementarity error: " << iterate.complementarity_error << '\n';
   stream << "Infeasibility measure: " << iterate.nonlinear_progress.infeasibility << '\n';
   stream << "Optimality measure: " << iterate.nonlinear_progress.scaled_optimality << '\n';
   return stream;
}