#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Problem& problem, const Options& options) :
      // create the phase-1 feasibility problem (objective multiplier = 0) with elastic variables
      ConstraintRelaxationStrategy(problem, 0., options),
      // create the globalization strategies (one for each phase)
      phase_1_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, 4);

   // initialize the subproblem
   this->subproblem->initialize(statistics, this->original_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   this->compute_optimality_progress_measures(first_iterate);
   this->subproblem->compute_nonlinear_residuals(this->original_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, double trust_region_radius) {
   if (this->current_phase == OPTIMALITY) {
      return this->solve_optimality_problem(statistics, current_iterate, trust_region_radius);
   }
   else {
      return this->solve_feasibility_problem(statistics, current_iterate, trust_region_radius, std::nullopt);
   }
}

Direction FeasibilityRestoration::solve_optimality_problem(Statistics& statistics, Iterate& current_iterate, double trust_region_radius) {
   this->subproblem->build_current_subproblem(this->original_problem, current_iterate, this->original_problem.objective_sign, trust_region_radius);

   // solve the subproblem
   DEBUG << "Current iterate\n" << current_iterate << "\n";
   DEBUG << "Solving the optimality subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->original_problem, current_iterate);
   direction.objective_multiplier = this->original_problem.objective_sign;
   direction.norm = norm_inf(direction.x);
   DEBUG << direction << "\n";

   // infeasible subproblem: form the feasibility problem
   if (direction.status == INFEASIBLE) {
      // try to minimize the constraint violation by solving the feasibility subproblem
      //current_iterate.reset_evaluations();
      current_iterate.is_objective_computed = false;
      current_iterate.is_objective_gradient_computed = false;
      // correct the evaluations
      this->relaxed_problem.add_elastics_to_constraints(current_iterate.x, current_iterate.constraints);
      this->relaxed_problem.add_elastics_to_constraint_jacobian(current_iterate.constraint_jacobian);
      direction = this->solve_feasibility_problem(statistics, current_iterate, trust_region_radius, direction.x);
   }
   return direction;
}

Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
      const std::optional<std::vector<double>>& optional_phase_2_solution) {
   // form and solve the feasibility problem (with or without constraint partition)
   // register the proximal coefficient and reference point
   this->relaxed_problem.set_proximal_coefficient(this->subproblem->get_proximal_coefficient());
   this->relaxed_problem.set_proximal_reference_point(current_iterate.x);

   // build the objective model of the feasibility problem
   //initialize_vector(current_iterate.multipliers.constraints, 0.);
   this->relaxed_problem.reset_elastic_variables(current_iterate);
   this->subproblem->build_current_subproblem(this->relaxed_problem, current_iterate, 0., trust_region_radius);

   // start from the phase-2 solution
   if (optional_phase_2_solution.has_value()) {
      const std::vector<double>& phase_2_primal_direction = optional_phase_2_solution.value();
      this->subproblem->set_initial_point(phase_2_primal_direction);
   }

   DEBUG << "Current iterate\n" << current_iterate << "\n";
   DEBUG << "Solving the feasibility subproblem\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   feasibility_direction.norm = norm_inf(feasibility_direction.x, 0, this->original_problem.number_variables);
   // create constraint partition
   ConstraintPartition constraint_partition(this->original_problem.number_constraints);
   constraint_partition.infeasible = this->relaxed_problem.get_violated_linearized_constraints(feasibility_direction.x);
   feasibility_direction.constraint_partition = constraint_partition;
   DEBUG << feasibility_direction << "\n";
   assert(feasibility_direction.status == OPTIMAL && "The subproblem was not solved to optimality");
   return feasibility_direction;
}

void FeasibilityRestoration::compute_restoration_progress_measures(Iterate& iterate, const std::vector<size_t>& constraint_set) {
   iterate.evaluate_constraints(this->original_problem);
   const double feasibility_measure = this->original_problem.compute_constraint_violation(iterate.constraints, L1_NORM);
   const double optimality_measure = this->original_problem.compute_constraint_violation(iterate.constraints, constraint_set, L1_NORM);
   iterate.nonlinear_progress = {feasibility_measure, optimality_measure};
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the progress measures are recomputed\n";
      this->subproblem->subproblem_definition_changed = false;
      this->phase_2_strategy->reset();
      this->compute_optimality_progress_measures(current_iterate);
   }

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      this->compute_optimality_progress_measures(trial_iterate);
      trial_iterate.evaluate_objective(this->original_problem);
      accept = true;
   }
   else {
      // possibly switch between phase 1 (restoration) and phase 2 (optimality)
      GlobalizationStrategy& current_phase_strategy = this->switch_phase(current_iterate, trial_iterate, direction);

      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.check_acceptance(statistics, current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress, direction.objective_multiplier,
            predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
      this->subproblem->compute_nonlinear_residuals(current_problem, trial_iterate);
      trial_iterate.evaluate_objective(this->original_problem);
   }
   return accept;
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) {
   if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->current_phase = FEASIBILITY_RESTORATION;
      DEBUG << "Switching from optimality to restoration phase\n";
      this->phase_2_strategy->notify(current_iterate);
      this->phase_1_strategy->reset();
      current_iterate.reset_evaluations();
      this->compute_restoration_progress_measures(current_iterate, direction.constraint_partition.value().infeasible);
      this->phase_1_strategy->notify(current_iterate);
   }
   // possibly go from 1 (restoration) to phase 2 (optimality)
   else if (this->current_phase == FEASIBILITY_RESTORATION && direction.constraint_partition.value().infeasible.empty()) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      trial_iterate.set_number_variables(this->original_problem.number_variables);
      current_iterate.reset_evaluations();
      this->compute_optimality_progress_measures(current_iterate);
   }

   // evaluate the progress measures of the trial iterate
   if (this->current_phase == OPTIMALITY) {
      this->compute_optimality_progress_measures(trial_iterate);
   }
   else {
      this->compute_restoration_progress_measures(trial_iterate, direction.constraint_partition.value().infeasible);
   }


   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

Direction FeasibilityRestoration::compute_second_order_correction(Iterate& trial_iterate) {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   return this->subproblem->compute_second_order_correction(current_problem, trial_iterate);
}

PredictedReductionModel FeasibilityRestoration::generate_predicted_reduction_model(const Iterate& current_iterate, const Direction& direction) const {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   return this->subproblem->generate_predicted_reduction_model(current_problem, current_iterate, direction);
}

void FeasibilityRestoration::register_accepted_iterate(Iterate& iterate) {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   this->subproblem->register_accepted_iterate(current_problem, iterate);
}