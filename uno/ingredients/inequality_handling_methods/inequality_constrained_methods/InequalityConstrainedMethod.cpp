// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedMethod.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "ingredients/subproblem_solvers/InequalityQPSolverFactory.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   InequalityConstrainedMethod::InequalityConstrainedMethod(size_t number_variables, size_t number_constraints,
      size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options):
         InequalityHandlingMethod(options.get_string("hessian_model"), options.get_string("regularization_strategy"), options),
         initial_point(number_variables),
         enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
         // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
         solver(InequalityQPSolverFactory::create(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
            // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
            std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, number_hessian_nonzeros),
            options)) {
      // TODO fix number of Hessian nnz
   }

    InequalityConstrainedMethod::~InequalityConstrainedMethod() { }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& statistics, const Options& options) {
      this->regularization_strategy->initialize_statistics(statistics, options);
   }

   void InequalityConstrainedMethod::generate_initial_iterate(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         Iterate& /*initial_iterate*/) {
      if (this->enforce_linear_constraints_at_initial_iterate) {
         // Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, *this->solver);
      }
   }

   void InequalityConstrainedMethod::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      LagrangeNewtonSubproblem subproblem(problem, current_iterate, current_multipliers, *this->hessian_model, *this->regularization_strategy,
         trust_region_radius);
      this->solver->solve_inequality_constrained_QP(statistics, subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
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

   void InequalityConstrainedMethod::exit_feasibility_problem(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         Iterate& /*trial_iterate*/) {
      // do nothing
   }

   double InequalityConstrainedMethod::hessian_quadratic_product(const Vector<double>& primal_direction) const {
      return this->solver->hessian_quadratic_product(primal_direction);
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

   void InequalityConstrainedMethod::set_initial_point(const Vector<double>& point) {
      // copy the point into the member
      this->initial_point = point;
   }

   void InequalityConstrainedMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      // compute dual *displacements* (active-set methods usually compute the new duals, not the displacements)
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }
} // namespace