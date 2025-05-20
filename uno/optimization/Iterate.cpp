// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Iterate.hpp"
#include "Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   size_t Iterate::number_eval_objective = 0;
   size_t Iterate::number_eval_constraints = 0;
   size_t Iterate::number_eval_objective_gradient = 0;
   size_t Iterate::number_eval_jacobian = 0;

   Iterate::Iterate(size_t number_variables, size_t number_constraints) :
         number_variables(number_variables), number_constraints(number_constraints),
         primals(number_variables), multipliers(number_variables, number_constraints),
         model_evaluations(number_variables, number_constraints), residuals(number_variables) {
   }

   void Iterate::evaluate_objective(const Model& model) {
      if (!this->is_objective_computed) {
         // evaluate the objective
         this->model_evaluations.objective = model.evaluate_objective(this->primals);
         Iterate::number_eval_objective++;
         if (!is_finite(this->model_evaluations.objective)) {
            throw FunctionEvaluationError();
         }
         this->is_objective_computed = true;
      }
   }

   void Iterate::evaluate_constraints(const Model& model) {
      if (!this->are_constraints_computed) {
         if (model.is_constrained()) {
            // evaluate the constraints
            model.evaluate_constraints(this->primals, this->model_evaluations.constraints);
            Iterate::number_eval_constraints++;
            // check finiteness
            if (std::any_of(this->model_evaluations.constraints.cbegin(), this->model_evaluations.constraints.cend(), [](double constraint_j) {
               return !is_finite(constraint_j);
            })) {
               throw FunctionEvaluationError();
            }
         }
         this->are_constraints_computed = true;
      }
   }

   void Iterate::evaluate_objective_gradient(const Model& model) {
      if (!this->is_objective_gradient_computed) {
         this->model_evaluations.objective_gradient.clear();
         // evaluate the objective gradient
         model.evaluate_objective_gradient(this->primals, this->model_evaluations.objective_gradient);
         this->is_objective_gradient_computed = true;
         Iterate::number_eval_objective_gradient++;
      }
   }

   void Iterate::evaluate_constraint_jacobian(const Model& model) {
      if (!this->is_constraint_jacobian_computed) {
         this->model_evaluations.constraint_jacobian.clear();
         if (model.is_constrained()) {
            model.evaluate_constraint_jacobian(this->primals, this->model_evaluations.constraint_jacobian);
            Iterate::number_eval_jacobian++;
         }
         this->is_constraint_jacobian_computed = true;
      }
   }

   void Iterate::set_number_variables(size_t new_number_variables) {
      this->number_variables = new_number_variables;
   }

   std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
      stream << "Primal variables: " << view(iterate.primals, 0, iterate.number_variables) << '\n';
      stream << "            ┌ Constraint: " << view(iterate.multipliers.constraints, 0, iterate.number_constraints) << '\n';
      stream << "Multipliers │ Lower bound: " << view(iterate.multipliers.lower_bounds, 0, iterate.number_variables) << '\n';
      stream << "            └ Upper bound: " << view(iterate.multipliers.upper_bounds, 0, iterate.number_variables) << '\n';
      stream << "Objective value: " << iterate.model_evaluations.objective << '\n';
      stream << "Primal feasibility: " << iterate.primal_feasibility << '\n';
      stream << "          ┌ Stationarity: " << iterate.residuals.stationarity << '\n';
      stream << "Residuals │ Complementarity: " << iterate.residuals.complementarity << '\n';
      stream << "          └ Lagrangian gradient: " << iterate.residuals.lagrangian_gradient;
      stream << "                  ┌ Infeasibility: " << iterate.progress.infeasibility << '\n';
      stream << "Progress measures │ Optimality: " << iterate.progress.objective(1.) << '\n';
      stream << "                  └ Auxiliary terms: " << iterate.progress.auxiliary << '\n';
      return stream;
   }
} // namespace