// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InequalityConstrainedMethod.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/BoxLPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   InequalityConstrainedMethod::InequalityConstrainedMethod(const Options& options):
         InequalityHandlingMethod(options), options(options) {
   }

   void InequalityConstrainedMethod::initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) {
      this->problem = &problem; // store the problem as is (no reformulation)
      assert(this->problem != nullptr);
      this->initial_point.resize(this->problem->number_variables);

      // allocate the LP/QP solver, depending on the presence of curvature in the subproblem
      const Subproblem subproblem{*this->problem, current_iterate, hessian_model, inertia_correction_strategy, trust_region_radius};
      if (!subproblem.has_curvature()) {
         if (subproblem.number_constraints == 0) {
            DEBUG << "No curvature and only bound constraints in the subproblems, allocating a box LP solver\n";
            this->solver = BoxLPSolverFactory::create();
         }
         else {
            DEBUG << "No curvature in the subproblems, allocating an LP solver\n";
            this->solver = LPSolverFactory::create(this->options);
         }
      }
      else {
         DEBUG << "Curvature in the subproblems, allocating a QP solver\n";
         this->solver = QPSolverFactory::create(this->options);
      }
      this->solver->initialize_memory(subproblem);
   }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::generate_initial_iterate(Iterate& initial_iterate) {
      this->problem->set_progress_measures(initial_iterate, this->progress_norm);
   }

   void InequalityConstrainedMethod::solve(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) {
      // create the subproblem and solve it
      const Subproblem subproblem{*this->problem, current_iterate, hessian_model, inertia_correction_strategy, trust_region_radius};
      this->solver->solve(statistics, subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_iterate.multipliers, direction.multipliers);
      ++this->number_subproblems_solved;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   void InequalityConstrainedMethod::initialize_feasibility_problem(Iterate& /*current_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
      problem.set_elastic_variable_values(current_iterate, [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         iterate.primals[elastic_index] = 0.;
         iterate.multipliers.lower_bounds[elastic_index] = 1.;
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double InequalityConstrainedMethod::proximal_coefficient() const {
      return 0.;
   }

   ProgressMeasures InequalityConstrainedMethod::compute_predicted_reductions(HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius, Iterate& current_iterate,
         const Direction& direction, double step_length) const {
      const Subproblem subproblem{*this->problem, current_iterate, hessian_model, inertia_correction_strategy, trust_region_radius};
      return subproblem.compute_predicted_reductions(this->get_evaluation_space(), direction, step_length, this->progress_norm);
   }

   bool InequalityConstrainedMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy,
         double trust_region_radius, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, UserCallbacks& user_callbacks) {
      // call the protected member function with the problem
      return InequalityHandlingMethod::is_iterate_acceptable(statistics, *this->problem, globalization_strategy,
         hessian_model, inertia_correction_strategy, trust_region_radius, current_iterate, trial_iterate, direction,
         step_length, user_callbacks);
   }

   void InequalityConstrainedMethod::evaluate_constraint_jacobian(Iterate& iterate) {
      auto& evaluation_space = this->solver->get_evaluation_space();
      evaluation_space.evaluate_constraint_jacobian(*this->problem, iterate);
   }

   void InequalityConstrainedMethod::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const auto& evaluation_space = this->solver->get_evaluation_space();
      evaluation_space.compute_constraint_jacobian_vector_product(vector, result);
   }

   void InequalityConstrainedMethod::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const auto& evaluation_space = this->solver->get_evaluation_space();
      evaluation_space.compute_constraint_jacobian_transposed_vector_product(vector, result);
   }

   double InequalityConstrainedMethod::compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const {
      const auto& evaluation_space = this->solver->get_evaluation_space();
      return evaluation_space.compute_hessian_quadratic_product(subproblem, vector);
   }

   // compute dual *displacements*
   // because of the way we form LPs/QPs, we get the new *multipliers* back from the solver. To get the dual displacements/direction,
   // we need to subtract the current multipliers
   void InequalityConstrainedMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }

   void InequalityConstrainedMethod::postprocess_iterate(Iterate& /*iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_initial_point(const Vector<double>& point) {
      // copy the point into the member
      this->initial_point = point;
   }

   std::string InequalityConstrainedMethod::get_name() const {
      return "inequality-constrained method";
   }

   EvaluationSpace& InequalityConstrainedMethod::get_evaluation_space() const {
      return this->solver->get_evaluation_space();
   }
} // namespace