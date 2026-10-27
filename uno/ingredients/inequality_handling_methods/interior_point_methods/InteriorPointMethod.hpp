// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTERIORPOINTMETHOD_H
#define UNO_INTERIORPOINTMETHOD_H

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include "../InequalityHandlingMethod.hpp"
#include "BarrierParameterUpdateStrategy.hpp"
#include "ElasticInitialization.hpp"
#include "InteriorPointParameters.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/EQPSolver.hpp"
#include "options/Options.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename BarrierProblem>
   class InteriorPointMethod : public InequalityHandlingMethod {
   public:
      InteriorPointMethod(const OptimizationProblem& problem, bool uses_trust_region, double objective_multiplier,
         const Options& options, std::vector<OptionOverride>& option_overrides);
      InteriorPointMethod(const InteriorPointMethod&) = delete;
      InteriorPointMethod& operator=(const InteriorPointMethod&) = delete;

      void initialize_memory() override;
      void initialize_statistics(Statistics& statistics) override;
      void create_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations, double multipliers_threshold) const override;
      void update_multipliers_of_square_problem(Iterate& iterate, Evaluations& evaluations, double multipliers_threshold) override;
      [[nodiscard]] std::pair<size_t, size_t> get_dimensions() const override;
      void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const override;

      [[nodiscard]] bool update_parameterization(Statistics& statistics, const Iterate& current_iterate,
         Evaluations& evaluations) override;
      [[nodiscard]] bool force_parameterization_update(Statistics& statistics) override;
      [[nodiscard]] const Direction& solve(Statistics& statistics, const Iterate& current_iterate, double trust_region_radius,
         const Vector<double>& initial_point, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      void import_handoff(const PhaseHandoff& handoff, Iterate& current_iterate) override;
      [[nodiscard]] PhaseHandoff export_handoff() const override;
      void relax_bounds_around(const Iterate& iterate) const override;
      void set_elastic_variable_values(const l1RelaxedProblem& feasibility_problem, Iterate& iterate, Evaluations& evaluations) override;
      [[nodiscard]] double proximal_coefficient() const override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      [[nodiscard]] const Direction& compute_second_order_correction(const Iterate& current_iterate) override;
      void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) override;

      void compute_bound_dual_direction(const Vector<double>& current_primals, const Multipliers& current_multipliers,
         const Vector<double>& direction_primals, Multipliers& direction_multipliers, double& bound_dual_step_length) const override;

      void evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const override;
      [[nodiscard]] PredictedReductionModels build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const override;
      [[nodiscard]] Decision is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, Evaluations& trial_evaluations,
         const ProgressMeasures& predicted_reductions) const override;
      void accept_iterate_unconditionally(Iterate& trial_iterate, Evaluations& trial_evaluations) const override;
      void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;

      void compute_residuals(Iterate& iterate, Evaluations& evaluations) const override;
      [[nodiscard]] double compute_primal_infeasibility(const Iterate& iterate, Evaluations& evaluations) const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const InteriorPointParameters parameters;
      BarrierParameterUpdateStrategy<BarrierProblem> barrier_parameter_update_strategy;
      BarrierProblem barrier_problem;

      std::unique_ptr<InertiaCorrectionStrategy> inertia_correction_strategy;
      std::unique_ptr<HessianModel> hessian_model;
      std::unique_ptr<SubproblemSolver> subproblem_solver;
      std::unique_ptr<Subproblem> subproblem;

      bool first_feasibility_iteration{false};
      const bool is_square_problem;

      [[nodiscard]] const double& barrier_parameter() const;
   };

   // class template implementation

   // IPOPT (IsSquareProblem): as many free variables as equality constraints
   inline bool is_square(const OptimizationProblem& problem) {
      const auto& lower_bounds = problem.get_variables_lower_bounds();
      const auto& upper_bounds = problem.get_variables_upper_bounds();
      size_t number_free_variables = 0;
      for (size_t variable_index: Range(problem.number_variables)) {
         if (lower_bounds[variable_index] != upper_bounds[variable_index]) {
            ++number_free_variables;
         }
      }
      return (0 < number_free_variables && number_free_variables == problem.get_equality_constraints().size());
   }


   template <typename BarrierProblem>
   InteriorPointMethod<BarrierProblem>::InteriorPointMethod(const OptimizationProblem& problem, bool uses_trust_region,
      double objective_multiplier, const Options& options, std::vector<OptionOverride>& option_overrides):
         InequalityHandlingMethod(problem, options),
         parameters({
            options.get_double("barrier_tau_min"),
            options.get_double("barrier_k_sigma"),
            options.get_double("barrier_regularization_exponent"),
            options.get_double("barrier_push_variable_to_interior_k1"),
            options.get_double("barrier_push_variable_to_interior_k2"),
            options.get_double("barrier_damping_factor"),
            options.get_double("barrier_default_multiplier"),
            std::pow(std::numeric_limits<double>::epsilon(), 0.75) // slack_move
         }),
         barrier_parameter_update_strategy(options),
         barrier_problem(problem, this->parameters, this->barrier_parameter(), options.get_double("residual_scaling_threshold")),
         is_square_problem(is_square(problem)) {
      // create the ingredients
      std::tie(this->inertia_correction_strategy, this->hessian_model, this->subproblem_solver) =
         HessianSubproblemSolverJointFactory::create(this->barrier_problem, uses_trust_region, objective_multiplier,
            options, option_overrides);
      this->subproblem = std::make_unique<Subproblem>(this->barrier_problem, *this->hessian_model, *this->inertia_correction_strategy);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::update_multipliers_of_square_problem(Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) {
      if (this->is_square_problem) {
         this->subproblem_solver->compute_least_squares_multipliers(*this->subproblem, iterate, evaluations,
            multipliers_threshold);
      }
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_memory() {
      this->subproblem_solver->initialize_memory(*this->subproblem);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_statistics(Statistics& statistics) {
      this->hessian_model->initialize_statistics(statistics);
      this->inertia_correction_strategy->initialize_statistics(statistics);
      this->subproblem_solver->initialize_statistics(statistics);
      statistics.add_column("Barrier", Statistics::double_width, 2, /* is_extended = */ true);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::create_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations,
         double multipliers_threshold) const {
      this->barrier_problem.create_initial_iterate(initial_iterate, evaluations);
      this->subproblem_solver->compute_least_squares_multipliers(*this->subproblem, initial_iterate, evaluations,
         multipliers_threshold);
   }

   template <typename BarrierProblem>
   std::pair<size_t, size_t> InteriorPointMethod<BarrierProblem>::get_dimensions() const {
      return {this->barrier_problem.number_variables, this->barrier_problem.number_constraints};
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      this->barrier_problem.postprocess_solution(iterate, evaluations);
   }

   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::update_parameterization(Statistics& statistics, const Iterate& current_iterate,
         Evaluations& evaluations) {
      bool update = false;
      // possibly update the barrier parameter
      if (!this->first_feasibility_iteration) {
         update = this->barrier_parameter_update_strategy.update_barrier_parameter(this->barrier_problem, current_iterate,
            evaluations);
      }
      else {
         this->first_feasibility_iteration = false;
      }
      statistics.set("Barrier", this->barrier_parameter());
      return update;
   }

   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::force_parameterization_update(Statistics& statistics) {
      const bool update = this->barrier_parameter_update_strategy.force_barrier_parameter_decrease();
      statistics.set("Barrier", this->barrier_parameter());
      return update;
   }

   template <typename BarrierProblem>
   const Direction& InteriorPointMethod<BarrierProblem>::solve(Statistics& statistics, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      DEBUG << "Barrier parameter = " << this->barrier_parameter() << '\n';
      return this->subproblem_solver->solve(statistics, *this->subproblem, current_iterate, trust_region_radius,
         initial_point, current_evaluations, warmstart_information);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::import_handoff(const PhaseHandoff& handoff, Iterate& current_iterate) {
      this->first_feasibility_iteration = true;

      // update the barrier parameter
      const double new_barrier_parameter = std::max(handoff.barrier_parameter.value_or(this->barrier_parameter()),
         current_iterate.scaled_residuals.primal_infeasibility);
      this->barrier_parameter_update_strategy.set_barrier_parameter(new_barrier_parameter);
      DEBUG << "Barrier parameter = " << this->barrier_parameter() << '\n';
   }

   template <typename BarrierProblem>
   PhaseHandoff InteriorPointMethod<BarrierProblem>::export_handoff() const {
      return {this->barrier_parameter()};
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::relax_bounds_around(const Iterate& iterate) const {
      this->barrier_problem.relax_bounds_around(iterate);
   }

   template <typename BarrierProblem>
   // set the elastic variables of the current iterate
   void InteriorPointMethod<BarrierProblem>::set_elastic_variable_values(const l1RelaxedProblem& feasibility_problem,
         Iterate& iterate, Evaluations& evaluations) {
      DEBUG << "IPM: setting the elastic variables and their duals\n";
      const double rho = feasibility_problem.constraint_violation_coefficient;
      const double mu = this->barrier_parameter();

      // cap the bound multipliers to [-rho, +rho] (IPOPT: z <- min(rho, z))
      for (size_t variable_index: Range(this->barrier_problem.number_variables)) {
         if (iterate.multipliers.lower_bounds[variable_index] > 0.) {
            iterate.multipliers.lower_bounds[variable_index] = std::min(iterate.multipliers.lower_bounds[variable_index], rho);
         }
         if (iterate.multipliers.upper_bounds[variable_index] < 0.) {
            iterate.multipliers.upper_bounds[variable_index] = std::max(iterate.multipliers.upper_bounds[variable_index], -rho);
         }
      }

      // reset the elastics and register the positive (p) and negative (n) elastics of each constraint
      static constexpr size_t no_variable = std::numeric_limits<size_t>::max();
      const size_t number_constraints = this->barrier_problem.number_constraints;
      std::vector<size_t> positive_elastic(number_constraints, no_variable);
      std::vector<size_t> negative_elastic(number_constraints, no_variable);
      feasibility_problem.set_elastic_variable_values([&](size_t constraint_index, size_t elastic_index, ElasticType elastic_type) {
         iterate.primals[elastic_index] = 0.;
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
         ((elastic_type == ElasticType::POSITIVE) ? positive_elastic : negative_elastic)[constraint_index] = elastic_index;
      });
      std::vector<size_t> slack_of_constraint(number_constraints, no_variable);
      for (const auto [constraint_index, slack_index]: this->barrier_problem.get_slacks()) {
         slack_of_constraint[constraint_index] = slack_index;
      }

      // residuals of the barrier constraints without elastics: r = c(x) - s (inequality) or r = c(x) - b (equality)
      Vector<double> residuals(number_constraints);
      this->barrier_problem.evaluate_constraints(iterate, residuals.view(), evaluations);

      // the elastics p (Jacobian coefficient -1) and n (+1) of constraint j satisfy the restoration stationarity
      // rho + y - z_p = 0 and rho - y - z_n = 0 with z_p p = z_n n = mu. The barrier constraint is
      // - equality (b = l = u), p and n:                    c - b - p + n = 0
      // - two-sided inequality (l < u finite), p and n:     c - s - p + n = 0,  l <= s <= u
      // - upper-bounded inequality (u finite), p only:      c - s - p = 0,      s <= u
      // - lower-bounded inequality (l finite), n only:      c - s + n = 0,      l <= s
      // Linear constraints that are not relaxed have no elastics and are left untouched
      const auto& slack_lower_bounds = this->barrier_problem.get_internal_variables_lower_bounds();
      const auto& slack_upper_bounds = this->barrier_problem.get_internal_variables_upper_bounds();
      for (size_t constraint_index: Range(number_constraints)) {
         const size_t p_index = positive_elastic[constraint_index];
         const size_t n_index = negative_elastic[constraint_index];
         if (p_index != no_variable && n_index != no_variable) {
            // the slack (if any) is fixed: p - n = r and z_p + z_n = 2 rho (eliminating y)
            const auto [p, n] = solve_elastic_pair(residuals[constraint_index], mu / rho);
            iterate.primals[p_index] = p;
            iterate.primals[n_index] = n;
            iterate.multipliers.lower_bounds[p_index] = mu / p;
            iterate.multipliers.lower_bounds[n_index] = mu / n;
         }
         else if (p_index != no_variable || n_index != no_variable) {
            // single elastic: the constraint is an inequality and the slack must move with the elastic, otherwise the
            // constraint cannot hold when the sign of r is wrong. The slack stationarity y + z_s = 0 (z_s: multiplier of
            // its only bound) and z_s (distance of s to its bound) = mu give, with e the elastic and d the distance:
            //    e - d = violation of the bound by c,   z_e + |z_s| = rho
            const size_t slack_index = slack_of_constraint[constraint_index];
            if (slack_index == no_variable) {
               throw std::logic_error("IPM: a constraint with a single elastic variable has no slack");
            }
            const double constraint_value = residuals[constraint_index] + iterate.primals[slack_index];
            if (p_index != no_variable) { // c - p = s <= u
               const double upper_bound = slack_upper_bounds[slack_index];
               const auto [p, distance] = solve_elastic_pair(constraint_value - upper_bound, 2. * mu / rho);
               iterate.primals[p_index] = p;
               iterate.multipliers.lower_bounds[p_index] = mu / p;
               iterate.primals[slack_index] = upper_bound - distance;
               iterate.multipliers.upper_bounds[slack_index] = -mu / distance;
            }
            else { // c + n = s >= l
               const double lower_bound = slack_lower_bounds[slack_index];
               const auto [n, distance] = solve_elastic_pair(lower_bound - constraint_value, 2. * mu / rho);
               iterate.primals[n_index] = n;
               iterate.multipliers.lower_bounds[n_index] = mu / n;
               iterate.primals[slack_index] = lower_bound + distance;
               iterate.multipliers.lower_bounds[slack_index] = mu / distance;
            }
         }
      }
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::proximal_coefficient() const {
      return std::sqrt(this->barrier_parameter());
   }

   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::has_second_order_corrections() const {
      return this->subproblem_solver->has_second_order_corrections();
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_second_order_corrections(const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      this->subproblem_solver->initialize_second_order_corrections(*this->subproblem, current_iterate, trial_iterate,
         current_evaluations, trial_evaluations);
   }

   template <typename BarrierProblem>
   const Direction& InteriorPointMethod<BarrierProblem>::compute_second_order_correction(const Iterate& current_iterate) {
      return this->subproblem_solver->compute_second_order_correction(*this->subproblem, current_iterate);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::update_second_order_corrections(const Iterate& trial_iterate,
         Evaluations& trial_evaluations) {
      this->subproblem_solver->update_second_order_corrections(*this->subproblem, trial_iterate, trial_evaluations);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::compute_bound_dual_direction(const Vector<double>& current_primals,
         const Multipliers& current_multipliers, const Vector<double>& direction_primals, Multipliers& direction_multipliers,
         double& bound_dual_step_length) const {
      this->barrier_problem.compute_bound_dual_direction(current_primals, current_multipliers, direction_primals,
         direction_multipliers, bound_dual_step_length);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const {
      InequalityHandlingMethod::evaluate_progress_measures(this->barrier_problem, iterate, evaluations);
   }

   template <typename BarrierProblem>
   PredictedReductionModels InteriorPointMethod<BarrierProblem>::build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const {
      return this->subproblem->build_predicted_reduction_model(current_iterate, direction, this->progress_norm,
         current_evaluations, this->subproblem_solver->get_workspace());
   }

   template <typename BarrierProblem>
   Decision InteriorPointMethod<BarrierProblem>::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, Evaluations& trial_evaluations,
         const ProgressMeasures& predicted_reductions) const {
      return InequalityHandlingMethod::is_iterate_acceptable(statistics, globalization_strategy, *this->subproblem,
         current_iterate, trial_iterate, direction, trial_evaluations, predicted_reductions);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::accept_iterate_unconditionally(Iterate& trial_iterate, Evaluations& trial_evaluations) const {
      InequalityHandlingMethod::accept_iterate_unconditionally(*this->subproblem, trial_iterate, trial_evaluations);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      this->hessian_model->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations, trial_evaluations);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::compute_residuals(Iterate& iterate, Evaluations& evaluations) const {
      InequalityHandlingMethod::compute_residuals(this->barrier_problem, iterate, evaluations);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::compute_primal_infeasibility(const Iterate& iterate, Evaluations& evaluations) const {
      return this->barrier_problem.primal_infeasibility(iterate, evaluations, this->residual_norm);
   }

   // protected member functions

   template <typename BarrierProblem>
   const double& InteriorPointMethod<BarrierProblem>::barrier_parameter() const {
      return this->barrier_parameter_update_strategy.get_barrier_parameter();
   }

   template <typename BarrierProblem>
   std::string InteriorPointMethod<BarrierProblem>::get_name() const {
      return "primal-dual interior-point method with " + this->hessian_model->name + " Hessian and "
         + this->inertia_correction_strategy->get_name() + " inertia correction";
   }
} // namespace

#endif // UNO_INTERIORPOINTMETHOD_H