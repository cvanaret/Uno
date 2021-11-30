#include "Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

size_t Iterate::number_eval_objective = 0;
size_t Iterate::number_eval_constraints = 0;
size_t Iterate::number_eval_jacobian = 0;

Iterate::Iterate(size_t number_variables, size_t number_constraints) :
   x(number_variables),
   multipliers(number_variables, number_constraints),
   constraints(multipliers.constraints.size()),
   objective_gradient(number_variables),
   constraints_jacobian(number_constraints),
   lagrangian_gradient(number_variables) {
}

void Iterate::evaluate_objective(const Problem& problem, const Scaling& scaling) {
   if (!this->is_objective_computed) {
      this->objective = problem.evaluate_objective(this->x);
      this->objective *= scaling.get_objective_scaling();
      this->is_objective_computed = true;
      Iterate::number_eval_objective++;
   }
}

void Iterate::evaluate_constraints(const Problem& problem, const Scaling& scaling) {
   if (!this->are_constraints_computed) {
      problem.evaluate_constraints(this->x, this->constraints);
      for (size_t j = 0; j < problem.number_constraints; j++) {
         this->constraints[j] *= scaling.get_constraint_scaling(j);
      }
      this->are_constraints_computed = true;
      Iterate::number_eval_constraints++;
   }
}

void Iterate::evaluate_objective_gradient(const Problem& problem, const Scaling& scaling) {
   if (!this->is_objective_gradient_computed) {
      this->objective_gradient.clear();
      problem.evaluate_objective_gradient(this->x, this->objective_gradient);
      scale(this->objective_gradient, scaling.get_objective_scaling());
      this->is_objective_gradient_computed = true;
   }
}

void Iterate::evaluate_constraints_jacobian(const Problem& problem, const Scaling& scaling) {
   if (!this->is_constraints_jacobian_computed) {
      for (auto& row: this->constraints_jacobian) {
         row.clear();
      }
      problem.evaluate_constraint_jacobian(this->x, this->constraints_jacobian);
      for (size_t j = 0; j < problem.number_constraints; j++) {
         scale(this->constraints_jacobian[j], scaling.get_constraint_scaling(j));
      }
      this->is_constraints_jacobian_computed = true;
      Iterate::number_eval_jacobian++;
   }
}

void Iterate::evaluate_lagrangian_gradient(const Problem& problem, const Scaling& scaling, double objective_multiplier, const Multipliers& multipliers) {
   initialize_vector(this->lagrangian_gradient, 0.);

   // objective gradient
   if (objective_multiplier != 0.) {
      this->evaluate_objective_gradient(problem, scaling);

      // scale the objective gradient
      this->objective_gradient.for_each([&](size_t i, double derivative) {
         if (i < problem.number_variables) {
            this->lagrangian_gradient[i] += objective_multiplier * derivative;
         }
      });
   }
   // bound constraints
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->lagrangian_gradient[i] -= multipliers.lower_bounds[i] + multipliers.upper_bounds[i];
   }

   // constraints
   this->evaluate_constraints_jacobian(problem, scaling);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double multiplier_j = multipliers.constraints[j];
      if (multiplier_j != 0.) {
         this->constraints_jacobian[j].for_each([&](size_t i, double derivative) {
            if (i < problem.number_variables) {
               this->lagrangian_gradient[i] -= multiplier_j * derivative;
            }
         });
      }
   }
}

void Iterate::change_number_variables(size_t number_variables) {
   this->x.resize(number_variables);
   this->multipliers.lower_bounds.resize(number_variables);
   this->multipliers.upper_bounds.resize(number_variables);
   this->objective_gradient.reserve(number_variables);
   this->lagrangian_gradient.resize(number_variables);
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

   stream << "Constraint residual: " << iterate.errors.constraints << "\n";
   stream << "KKT residual: " << iterate.errors.KKT << "\n";
   stream << "FJ residual: " << iterate.errors.KKT << "\n";
   stream << "Complementarity residual: " << iterate.errors.complementarity << "\n";

   stream << "Optimality measure: " << iterate.progress.objective << "\n";
   stream << "Feasibility measure: " << iterate.progress.infeasibility << "\n";
   return stream;
}