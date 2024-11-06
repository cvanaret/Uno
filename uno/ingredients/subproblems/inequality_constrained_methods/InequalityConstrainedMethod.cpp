// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedMethod.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   InequalityConstrainedMethod::InequalityConstrainedMethod(const std::string& hessian_model, size_t number_variables, size_t number_constraints,
         size_t number_hessian_nonzeros, bool convexify, const Options& options):
         Subproblem(hessian_model, number_variables, number_hessian_nonzeros, convexify, options),
         initial_point(number_variables),
         direction_lower_bounds(number_variables),
         direction_upper_bounds(number_variables),
         linearized_constraints_lower_bounds(number_constraints),
         linearized_constraints_upper_bounds(number_constraints),
         objective_gradient(number_variables),
         constraints(number_constraints),
         constraint_jacobian(number_constraints, number_variables) {
   }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) {
   }

   void InequalityConstrainedMethod::set_initial_point(const Vector<double>& point) {
      // copy the point into the member
      this->initial_point = point;
   }

   void InequalityConstrainedMethod::initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& /*current_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
      problem.set_elastic_variable_values(current_iterate, [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         iterate.primals[elastic_index] = 0.;
         iterate.feasibility_multipliers.lower_bounds[elastic_index] = 1.;
         iterate.feasibility_multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double InequalityConstrainedMethod::proximal_coefficient(const Iterate& /*current_iterate*/) const {
      return 0.;
   }

   void InequalityConstrainedMethod::exit_feasibility_problem(const OptimizationProblem& /*problem*/, Iterate& /*trial_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_direction_bounds(const OptimizationProblem& problem, const Iterate& current_iterate) {
      // bounds of original variables intersected with trust region
      for (size_t variable_index: Range(problem.get_number_original_variables())) {
         this->direction_lower_bounds[variable_index] = std::max(-this->trust_region_radius,
               problem.variable_lower_bound(variable_index) - current_iterate.primals[variable_index]);
         this->direction_upper_bounds[variable_index] = std::min(this->trust_region_radius,
               problem.variable_upper_bound(variable_index) - current_iterate.primals[variable_index]);
      }
      // bounds of additional variables (no trust region!)
      for (size_t variable_index: Range(problem.get_number_original_variables(), problem.number_variables)) {
         this->direction_lower_bounds[variable_index] = problem.variable_lower_bound(variable_index) - current_iterate.primals[variable_index];
         this->direction_upper_bounds[variable_index] = problem.variable_upper_bound(variable_index) - current_iterate.primals[variable_index];
      }
   }

   void InequalityConstrainedMethod::set_linearized_constraint_bounds(const OptimizationProblem& problem, const std::vector<double>& current_constraints) {
      for (size_t constraint_index: Range(problem.number_constraints)) {
         this->linearized_constraints_lower_bounds[constraint_index] = problem.constraint_lower_bound(constraint_index) -
               current_constraints[constraint_index];
         this->linearized_constraints_upper_bounds[constraint_index] = problem.constraint_upper_bound(constraint_index) -
               current_constraints[constraint_index];
      }
   }

   void InequalityConstrainedMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      // compute dual *displacements* (active-set methods usually compute the new duals, not the displacements)
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }

   // auxiliary measure is 0 in inequality-constrained methods
   void InequalityConstrainedMethod::set_auxiliary_measure(const Model& /*model*/, Iterate& iterate) {
      iterate.progress.auxiliary = 0.;
   }

   double InequalityConstrainedMethod::compute_predicted_auxiliary_reduction_model(const Model& /*model*/, const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/, double /*step_length*/) const {
      return 0.;
   }

   void InequalityConstrainedMethod::postprocess_iterate(const OptimizationProblem& /*problem*/, Iterate& /*iterate*/) {
   }
} // namespace
