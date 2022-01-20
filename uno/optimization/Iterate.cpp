#include "Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

size_t Iterate::number_eval_objective = 0;
size_t Iterate::number_eval_constraints = 0;
size_t Iterate::number_eval_jacobian = 0;

Iterate::Iterate(size_t max_number_variables, size_t max_number_constraints) :
   number_variables(max_number_variables), number_constraints(max_number_constraints),
   x(max_number_variables), multipliers(max_number_variables, max_number_constraints),
   constraints(multipliers.constraints.size()),
   subproblem_constraints(multipliers.constraints.size()),
   objective_gradient(max_number_variables),
   constraint_jacobian(max_number_constraints),
   lagrangian_gradient(max_number_variables) {
}

void Iterate::evaluate_objective(const Problem& problem) {
   if (!this->is_objective_computed) {
      // evaluate the objective
      this->objective = problem.evaluate_objective(this->x);
      this->is_objective_computed = true;
      Iterate::number_eval_objective++;
   }
}

void Iterate::evaluate_constraints(const Problem& problem) {
   if (!this->are_constraints_computed) {
      // evaluate the constraints
      problem.evaluate_constraints(this->x, this->constraints);
      this->are_constraints_computed = true;
      Iterate::number_eval_constraints++;
   }
}

void Iterate::evaluate_objective_gradient(const Problem& problem) {
   if (!this->is_objective_gradient_computed) {
      this->objective_gradient.clear();
      // evaluate the objective gradient
      problem.evaluate_objective_gradient(this->x, this->objective_gradient);
      this->is_objective_gradient_computed = true;
   }
}

void Iterate::evaluate_constraint_jacobian(const Problem& problem) {
   if (!this->is_constraint_jacobian_computed) {
      for (auto& row: this->constraint_jacobian) {
         row.clear();
      }
      // evaluate the constraint Jacobian
      problem.evaluate_constraint_jacobian(this->x, this->constraint_jacobian);
      this->is_constraint_jacobian_computed = true;
      Iterate::number_eval_jacobian++;
   }
}

void Iterate::evaluate_lagrangian_gradient(const Problem& problem, double objective_multiplier, const std::vector<double>& constraint_multipliers,
      const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers) {
   initialize_vector(this->lagrangian_gradient, 0.);

   // objective gradient
   if (objective_multiplier != 0.) {
      this->evaluate_objective_gradient(problem);

      // scale the objective gradient
      this->objective_gradient.for_each([&](size_t i, double derivative) {
         // in case there are additional variables, ignore them
         if (i < problem.number_variables) {
            this->lagrangian_gradient[i] += objective_multiplier * derivative;
         }
      });
   }
   // bound constraints
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->lagrangian_gradient[i] -= lower_bounds_multipliers[i] + upper_bounds_multipliers[i];
   }

   // constraints
   this->evaluate_constraint_jacobian(problem);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      if (multiplier_j != 0.) {
         this->constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            // in case there are additional variables, ignore them
            if (i < problem.number_variables) {
               this->lagrangian_gradient[i] -= multiplier_j * derivative;
            }
         });
      }
   }
}

void Iterate::set_number_variables(size_t new_number_variables) {
   this->x.resize(new_number_variables);
   this->multipliers.lower_bounds.resize(new_number_variables);
   this->multipliers.upper_bounds.resize(new_number_variables);
   this->objective_gradient.reserve(new_number_variables);
   this->lagrangian_gradient.resize(new_number_variables);
}

void Iterate::reset_evaluations() {
   this->is_objective_computed = false;
   this->is_objective_gradient_computed = false;
   this->are_constraints_computed = false;
   this->is_constraint_jacobian_computed = false;
}

std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
   stream << "x: ";
   print_vector(stream, iterate.x);
   stream << "Lower bound multipliers: ";
   print_vector(stream, iterate.multipliers.lower_bounds);
   stream << "Upper bound multipliers: ";
   print_vector(stream, iterate.multipliers.upper_bounds);
   stream << "Constraint multipliers: ";
   print_vector(stream, iterate.multipliers.constraints);
   stream << "Objective value: " << iterate.objective << "\n";

   stream << "Constraint residual: " << iterate.nonlinear_errors.constraints << "\n";
   stream << "KKT residual: " << iterate.nonlinear_errors.KKT << "\n";
   stream << "FJ residual: " << iterate.nonlinear_errors.KKT << "\n";
   stream << "Complementarity residual: " << iterate.nonlinear_errors.complementarity << "\n";

   stream << "Optimality measure: " << iterate.progress.objective << "\n";
   stream << "Feasibility measure: " << iterate.progress.infeasibility << "\n";
   return stream;
}