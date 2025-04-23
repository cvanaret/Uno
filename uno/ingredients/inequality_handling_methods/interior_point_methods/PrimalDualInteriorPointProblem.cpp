// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem, double barrier_parameter):
         OptimizationProblem(problem.model, problem.number_variables, problem.number_constraints),
         first_reformulation(problem), barrier_parameter(barrier_parameter) { }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->first_reformulation.get_objective_multiplier();
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
      this->first_reformulation.evaluate_objective_gradient(iterate, objective_gradient);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(first_reformulation.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - first_reformulation.variable_lower_bound(variable_index));
            // damping
            if (!is_finite(first_reformulation.variable_upper_bound(variable_index))) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(first_reformulation.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - first_reformulation.variable_upper_bound(variable_index));
            // damping
            if (!is_finite(first_reformulation.variable_lower_bound(variable_index))) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         objective_gradient.insert(variable_index, barrier_term);
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      this->first_reformulation.evaluate_constraints(iterate, constraints);
   }

   void PrimalDualInteriorPointProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      this->first_reformulation.evaluate_constraint_jacobian(iterate, constraint_jacobian);
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const {
      // original Lagrangian Hessian
      this->first_reformulation.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian);
      hessian.set_dimension(this->number_variables);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool finite_lower_bound = is_finite(this->first_reformulation.variable_lower_bound(variable_index));
         const bool finite_upper_bound = is_finite(this->first_reformulation.variable_upper_bound(variable_index));
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->first_reformulation.variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->first_reformulation.variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            hessian.insert(diagonal_barrier_term, variable_index, variable_index);
         }
      }
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, const Vector<double>& vector,
         const Multipliers& multipliers, Vector<double>& result) const {
      // original Lagrangian Hessian
      this->first_reformulation.compute_hessian_vector_product(hessian_model, vector, multipliers, result);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool finite_lower_bound = is_finite(this->first_reformulation.variable_lower_bound(variable_index));
         const bool finite_upper_bound = is_finite(this->first_reformulation.variable_upper_bound(variable_index));
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) { // lower bounded
               const double distance_to_bound = vector[variable_index] - this->first_reformulation.variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) { // upper bounded
               const double distance_to_bound = vector[variable_index] - this->first_reformulation.variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            result[variable_index] += diagonal_barrier_term * vector[variable_index];
         }
      }
   }

   double PrimalDualInteriorPointProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalDualInteriorPointProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   double PrimalDualInteriorPointProblem::constraint_lower_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   double PrimalDualInteriorPointProblem::constraint_upper_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_lower_bounded_variables() const {
      return this->first_reformulation.get_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_upper_bounded_variables() const {
      return this->first_reformulation.get_upper_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_lower_bounded_variables() const {
      return this->first_reformulation.get_single_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_upper_bounded_variables() const {
      return this->first_reformulation.get_single_upper_bounded_variables();
   }

   size_t PrimalDualInteriorPointProblem::number_objective_gradient_nonzeros() const {
      size_t number_nonzeros = this->first_reformulation.number_objective_gradient_nonzeros();
      // barrier contribution
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->first_reformulation.variable_lower_bound(variable_index)) || is_finite(this->first_reformulation.variable_upper_bound(variable_index))) {
            number_nonzeros++;
         }
      }
      return number_nonzeros;
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->first_reformulation.number_jacobian_nonzeros();
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros() const {
      size_t number_zeros = this->first_reformulation.number_hessian_nonzeros();
      // barrier contribution
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->first_reformulation.variable_lower_bound(variable_index)) || is_finite(this->first_reformulation.variable_upper_bound(variable_index))) {
            number_zeros++;
         }
      }
      return number_zeros;
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const {
      this->first_reformulation.evaluate_lagrangian_gradient(lagrangian_gradient, iterate, multipliers);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(first_reformulation.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - first_reformulation.variable_lower_bound(variable_index));
            // damping
            if (!is_finite(first_reformulation.variable_upper_bound(variable_index))) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(first_reformulation.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - first_reformulation.variable_upper_bound(variable_index));
            // damping
            if (!is_finite(first_reformulation.variable_lower_bound(variable_index))) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient.constraints_contribution[variable_index] += barrier_term;
      }
   }

   double PrimalDualInteriorPointProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->first_reformulation.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }
} // namespace