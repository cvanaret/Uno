#include <limits>
#include "Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

int Iterate::number_eval_objective = 0;
int Iterate::number_eval_constraints = 0;
int Iterate::number_eval_jacobian = 0;

Iterate::Iterate(const std::vector<double>& x, const Multipliers& multipliers) :
   x(x),
   multipliers(multipliers),
   constraints(multipliers.constraints.size()),
   objective_gradient(x.size()),
   constraints_jacobian(multipliers.constraints.size(), x.size()) {
}

Iterate::Iterate(size_t number_variables, size_t number_constraints) :
   x(number_variables),
   multipliers(number_variables, number_constraints),
   constraints(multipliers.constraints.size()),
   objective_gradient(number_variables),
   constraints_jacobian(number_constraints, number_variables) {
}

void Iterate::compute_objective(const Problem& problem) {
   if (!this->is_objective_computed) {
      this->objective = problem.evaluate_objective(this->x);
      this->is_objective_computed = true;
      Iterate::number_eval_objective++;
   }
}

void Iterate::compute_constraints(const Problem& problem) {
   if (!this->are_constraints_computed) {
      problem.evaluate_constraints(this->x, this->constraints);
      this->are_constraints_computed = true;
      Iterate::number_eval_constraints++;
   }
}

void Iterate::compute_objective_gradient(const Problem& problem) {
   if (!this->is_objective_gradient_computed) {
      this->objective_gradient.clear();
      problem.evaluate_objective_gradient(this->x, this->objective_gradient);
      this->is_objective_gradient_computed = true;
   }
}

void Iterate::compute_constraints_jacobian(const Problem& problem) {
   if (!this->is_constraints_jacobian_computed) {
      for (auto& row: this->constraints_jacobian) {
         row.clear();
      }
      problem.constraints_jacobian(this->x, this->constraints_jacobian);
      this->is_constraints_jacobian_computed = true;
      Iterate::number_eval_jacobian++;
   }
}

std::vector<double> Iterate::lagrangian_gradient(const Problem& problem, double objective_multiplier, const Multipliers& multipliers) {
   std::vector<double> lagrangian_gradient(problem.number_variables);

   /* objective gradient */
   if (objective_multiplier != 0.) {
      this->compute_objective_gradient(problem);

      /* scale the objective gradient */
      this->objective_gradient.for_each([&](size_t i, double derivative) {
         if (i < problem.number_variables) {
            lagrangian_gradient[i] += objective_multiplier * derivative;
         }
      });
   }
   /* bound constraints */
   for (size_t i = 0; i < problem.number_variables; i++) {
      lagrangian_gradient[i] -= multipliers.lower_bounds[i] + multipliers.upper_bounds[i];
   }

   /* constraints */
   this->compute_constraints_jacobian(problem);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double multiplier_j = multipliers.constraints[j];
      if (multiplier_j != 0.) {
         this->constraints_jacobian[j].for_each([&](size_t i, double derivative) {
            if (i < problem.number_variables) {
               lagrangian_gradient[i] -= multiplier_j * derivative;
            }
         });
      }
   }
   return lagrangian_gradient;
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

void reset(Iterate& iterate) {
   iterate.is_objective_computed = false;
   iterate.are_constraints_computed = false;
   iterate.is_objective_gradient_computed = false;
   iterate.is_constraints_jacobian_computed = false;
   iterate.errors = {0., 0., 0., 0.};
   iterate.progress = {0., 0.};
}