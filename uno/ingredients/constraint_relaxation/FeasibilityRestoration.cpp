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
   this->subproblem->compute_progress_measures(this->original_problem, first_iterate);
   this->subproblem->compute_residuals(this->original_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

void FeasibilityRestoration::create_current_subproblem(Iterate& current_iterate, double trust_region_radius) {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   this->subproblem->build_current_subproblem(current_problem, current_iterate, current_problem.objective_sign, trust_region_radius);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   // solve the subproblem
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   Direction direction = this->subproblem->solve(statistics, current_problem, current_iterate);
   direction.objective_multiplier = (this->current_phase == OPTIMALITY) ? this->original_problem.objective_sign : 0.;
   DEBUG << direction << "\n";

   // infeasible subproblem: form the feasibility problem
   if (this->current_phase == OPTIMALITY && direction.status == INFEASIBLE) {
      // try to minimize the constraint violation by solving the feasibility subproblem
      current_iterate.reset_evaluations();
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.x);
   }
   else {
      assert(direction.status == OPTIMAL && "The subproblem was not solved to optimality");
   }
   return direction;
}

Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_solution) {
   // make the current iterate larger
   current_iterate.set_number_variables(this->relaxed_problem.number_variables);
   // form and solve the feasibility problem (with or without constraint partition)
   this->create_current_feasibility_problem(current_iterate, optional_phase_2_solution);

   DEBUG << "Solving the feasibility subproblem at the current iterate:\n" << current_iterate << "\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   DEBUG << feasibility_direction << "\n";
   assert(feasibility_direction.status == OPTIMAL && "The subproblem was not solved to optimality");
   return feasibility_direction;
}

void FeasibilityRestoration::create_current_feasibility_problem(Iterate& current_iterate, const std::optional<std::vector<double>>& optional_phase_2_solution) {
   // register the proximal coefficient and reference point
   this->relaxed_problem.set_proximal_coefficient(this->subproblem->get_proximal_coefficient());
   this->relaxed_problem.set_proximal_reference_point(current_iterate.x);

   // build the objective model of the feasibility problem
   initialize_vector(current_iterate.multipliers.constraints, 0.);
   this->subproblem->build_current_subproblem(this->relaxed_problem, current_iterate, 0., 10.);
   //this->subproblem->build_objective_model(this->relaxed_problem, current_iterate, 0.);

   // start from the phase-2 solution
   if (optional_phase_2_solution.has_value()) {
      const std::vector<double>& phase_2_primal_direction = optional_phase_2_solution.value();
      this->subproblem->set_initial_point(phase_2_primal_direction);
   }
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the progress measures are recomputed\n";
      this->subproblem->subproblem_definition_changed = false;
      this->phase_2_strategy->reset();
      this->subproblem->compute_progress_measures(this->original_problem, current_iterate);
   }

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      this->subproblem->compute_progress_measures(this->original_problem, trial_iterate);
      accept = true;
   }
   else {
      // possibly switch between phase 1 (restoration) and phase 2 (optimality)
      GlobalizationStrategy& current_phase_strategy = this->switch_phase(current_iterate, trial_iterate, direction);

      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction.objective_multiplier,
            predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
      this->subproblem->compute_residuals(current_problem, trial_iterate);
   }
   return accept;
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) {
   // possibly go from 1 (restoration) to phase 2 (optimality)
   if (this->current_phase == FEASIBILITY_RESTORATION && this->relaxed_problem.compute_linearized_constraint_violation(current_iterate.x) == 0.) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      current_iterate.reset_evaluations();
      this->subproblem->compute_progress_measures(this->original_problem, current_iterate);
   }
   // possibly go from phase 2 (optimality) to 1 (restoration)
   else if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->current_phase = FEASIBILITY_RESTORATION;
      DEBUG << "Switching from optimality to restoration phase\n";
      this->phase_2_strategy->notify(current_iterate);
      this->phase_1_strategy->reset();
      current_iterate.reset_evaluations();
      this->subproblem->compute_progress_measures(this->original_problem, current_iterate);
      this->phase_1_strategy->notify(current_iterate);
   }

   // evaluate the progress measures of the trial iterate
   //const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   this->subproblem->compute_progress_measures(this->original_problem, trial_iterate);

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