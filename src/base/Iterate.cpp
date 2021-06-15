#include <limits>
#include <utility>
#include "Iterate.hpp"
#include "Vector.hpp"
#include "Logger.hpp"

int Iterate::number_eval_objective = 0;
int Iterate::number_eval_constraints = 0;
int Iterate::number_eval_jacobian = 0;

Iterate::Iterate(size_t number_variables, size_t number_constraints):
   x(number_variables), multipliers(number_variables, number_constraints), constraints(number_constraints) {
}

Iterate::Iterate(const std::vector<double>& x, const Multipliers& multipliers) : x(x), multipliers(multipliers),
      objective(std::numeric_limits<double>::infinity()), is_objective_computed(false),
      constraints(multipliers.constraints.size()), are_constraints_computed(false),
      is_objective_gradient_computed(false), is_constraints_jacobian_computed(false),
      //hessian(x.size(), 1), is_hessian_computed(false),
      residuals({0., 0., 0., 0.}),
      progress({0., 0.}) {
}

void Iterate::compute_objective(const Problem& problem) {
   if (!this->is_objective_computed) {
      this->objective = problem.objective(this->x);
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
      this->objective_gradient = problem.objective_gradient(this->x);
      this->is_objective_gradient_computed = true;
   }
}

void Iterate::set_objective_gradient(const SparseVector& objective_gradient) {
   this->objective_gradient = objective_gradient;
   this->is_objective_gradient_computed = true;
}

void Iterate::compute_constraints_jacobian(const Problem& problem) {
   if (!this->is_constraints_jacobian_computed) {
      this->constraints_jacobian = problem.constraints_jacobian(this->x);
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
      for (const auto[i, derivative]: this->objective_gradient) {
         if (i < problem.number_variables) {
            lagrangian_gradient[i] += objective_multiplier * derivative;
         }
      }
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
         for (const auto[i, derivative]: this->constraints_jacobian[j]) {
            if (i < problem.number_variables) {
               lagrangian_gradient[i] -= multiplier_j * derivative;
            }
         }
      }
   }
   return lagrangian_gradient;
}

void Iterate::clear() {
   this->is_objective_computed = false;
   this->is_objective_gradient_computed = false;
   this->is_constraints_jacobian_computed = false;
   this->residuals = {0., 0., 0., 0.};
   this->progress = {0., 0.};
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

   stream << "Constraint residual: " << iterate.residuals.constraints << "\n";
   stream << "KKT residual: " << iterate.residuals.KKT << "\n";
   stream << "FJ residual: " << iterate.residuals.KKT << "\n";
   stream << "Complementarity residual: " << iterate.residuals.complementarity << "\n";

   stream << "Optimality measure: " << iterate.progress.objective << "\n";
   stream << "Feasibility measure: " << iterate.progress.feasibility << "\n";
   return stream;
}