// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "NoInequalityReformulation.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategy.hpp"
#include "ingredients/subproblem_solvers/SubproblemSolver.hpp"
#include "optimization/Direction.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   NoInequalityReformulation::NoInequalityReformulation(std::string name, const OptimizationProblem& problem,
         bool uses_trust_region, double objective_multiplier, Options& options):
      InequalityHandlingMethod(problem, options), name(std::move(name)) {
      // create the ingredients
      std::tie(this->inertia_correction_strategy, this->hessian_model, this->subproblem_solver) =
         HessianSubproblemSolverJointFactory::create(this->problem, uses_trust_region, objective_multiplier, options);
      this->subproblem = std::make_unique<Subproblem>(this->problem, *this->hessian_model, *this->inertia_correction_strategy);
      this->subproblem_solver->initialize_memory(*this->subproblem);
   }

   void NoInequalityReformulation::generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const {
      this->problem.generate_initial_iterate(initial_iterate, evaluations);
      this->subproblem_solver->generate_initial_iterate(*this->subproblem, initial_iterate, evaluations);
   }

   void NoInequalityReformulation::initialize_statistics(Statistics& statistics) {
      this->hessian_model->initialize_statistics(statistics);
      this->inertia_correction_strategy->initialize_statistics(statistics);
   }

   bool NoInequalityReformulation::update_parameterization(Statistics& /*statistics*/, const Iterate& /*current_iterate*/) {
      // the parameterization is not updated
      return false;
   }

   const Direction& NoInequalityReformulation::solve(Statistics& statistics, const Iterate& current_iterate, double trust_region_radius,
         const Vector<double>& initial_point, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      return this->subproblem_solver->solve(statistics, *this->subproblem, current_iterate, trust_region_radius,
         initial_point, current_evaluations, warmstart_information);
   }

   void NoInequalityReformulation::initialize_feasibility_problem(Iterate& /*current_iterate*/) {
      // do nothing
   }

   void NoInequalityReformulation::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate,
         Evaluations& /*evaluations*/) {
      // c(x) - p + n = 0
      // TODO set (one of) the elastic variables to +/- the value of the constraint if violated
      problem.set_elastic_variable_values(current_iterate, [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         iterate.primals[elastic_index] = 0.;
         iterate.multipliers.lower_bounds[elastic_index] = 1.;
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double NoInequalityReformulation::proximal_coefficient() const {
      return 0.;
   }

   bool NoInequalityReformulation::has_second_order_corrections() const {
      return this->subproblem_solver->has_second_order_corrections();
   }

   const Direction& NoInequalityReformulation::compute_second_order_correction(const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) {
      return this->subproblem_solver->compute_second_order_correction(*this->subproblem, current_iterate, constraints_SOC);
   }

   void NoInequalityReformulation::compute_least_squares_multipliers(Iterate& iterate, Evaluations& evaluations) {
      // no threshold on the multipliers
      this->subproblem_solver->compute_least_squares_multipliers(*this->subproblem, iterate, evaluations, INF<double>);
   }

   void NoInequalityReformulation::evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const {
      InequalityHandlingMethod::evaluate_progress_measures(this->problem, iterate, evaluations);
   }

   bool NoInequalityReformulation::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) const {
      return InequalityHandlingMethod::is_iterate_acceptable(statistics, globalization_strategy, *this->subproblem,
         this->subproblem_solver->get_workspace(), current_iterate, trial_iterate, direction, step_length, current_evaluations,
         trial_evaluations);
   }

   void NoInequalityReformulation::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      this->hessian_model->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations, trial_evaluations);
   }

   std::string NoInequalityReformulation::get_name() const {
      return this->name + " with " + this->hessian_model->name + " Hessian and " + this->inertia_correction_strategy->get_name()
         + " inertia correction";
   }
} // namespace