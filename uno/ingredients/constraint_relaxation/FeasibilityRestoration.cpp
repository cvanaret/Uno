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
   first_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(first_iterate);
   first_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->original_problem, first_iterate);
   this->subproblem->compute_nonlinear_residuals(this->original_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Current iterate\n" << current_iterate << "\n";
   if (this->current_phase == OPTIMALITY) {
      return this->solve_optimality_problem(statistics, current_iterate);
   }
   else {
      return this->solve_feasibility_problem(statistics, current_iterate, std::nullopt);
   }
}

Direction FeasibilityRestoration::solve_optimality_problem(Statistics& statistics, Iterate& current_iterate) {
   this->subproblem->build_objective_model(this->original_problem, current_iterate, this->original_problem.objective_sign);
   this->subproblem->build_constraint_model(this->original_problem, current_iterate);
   current_iterate.subproblem_evaluations.objective_gradient = current_iterate.problem_evaluations.objective_gradient;
   current_iterate.subproblem_evaluations.constraints = current_iterate.problem_evaluations.constraints;
   current_iterate.subproblem_evaluations.constraint_jacobian = current_iterate.problem_evaluations.constraint_jacobian;

   // solve the subproblem
   DEBUG << "Solving the optimality subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->original_problem, current_iterate);
   direction.objective_multiplier = this->original_problem.objective_sign;
   direction.norm = norm_inf(direction.x, 0, this->original_problem.number_variables);
   DEBUG << direction << "\n";

   // infeasible subproblem: try to minimize the constraint violation by solving the feasibility subproblem
   if (direction.status == INFEASIBLE) {
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.x);
   }
   return direction;
}

// form and solve the feasibility problem (with or without constraint partition)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_solution) {
   // register the proximal coefficient and reference point
   this->relaxed_problem.set_proximal_coefficient(this->subproblem->get_proximal_coefficient());
   this->relaxed_problem.set_proximal_reference_point(current_iterate.x);

   // build the objective model of the feasibility problem
   //if (this->current_phase == OPTIMALITY) {
      this->subproblem->set_elastic_variables(this->relaxed_problem, current_iterate);
   //}
   this->relaxed_problem.evaluate_constraints(current_iterate);
   this->relaxed_problem.evaluate_constraint_jacobian(current_iterate);
   this->subproblem->build_objective_model(this->relaxed_problem, current_iterate, 0.);
   this->subproblem->build_constraint_model(this->relaxed_problem, current_iterate);

   // start from the phase-2 solution
   this->subproblem->set_initial_point(optional_phase_2_solution);

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

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the optimality measure is recomputed\n";
      current_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->original_problem, current_iterate);
      this->phase_2_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }

   // possibly switch between phase 1 (restoration) and phase 2 (optimality)
   GlobalizationStrategy& current_phase_strategy = this->switch_phase(current_iterate, trial_iterate, direction);

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.check_acceptance(statistics, current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
      trial_iterate.evaluate_objective(this->original_problem);
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
      // update the measure of optimality
      const std::vector<size_t>& infeasible_constraints = direction.constraint_partition.value().infeasible;
      current_iterate.nonlinear_progress.objective = this->compute_optimality_measure(current_iterate, infeasible_constraints);
      this->phase_1_strategy->notify(current_iterate);
   }
   // possibly go from 1 (restoration) to phase 2 (optimality)
   else if (this->current_phase == FEASIBILITY_RESTORATION && direction.constraint_partition.value().infeasible.empty()) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      current_iterate.set_number_variables(this->original_problem.number_variables);
      current_iterate.reset_evaluations();
      current_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->original_problem, current_iterate);
      trial_iterate.set_number_variables(this->original_problem.number_variables);
   }

   // evaluate the progress measures of the trial iterate
   trial_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(trial_iterate);
   if (this->current_phase == OPTIMALITY) {
      trial_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->original_problem, trial_iterate);
   }
   else {
      trial_iterate.nonlinear_progress.objective = this->compute_optimality_measure(trial_iterate, direction.constraint_partition.value().infeasible);
   }

   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) {
   // set the bounds of all the variables (primal + elastics)
   this->subproblem->set_variable_bounds(this->relaxed_problem, current_iterate, trust_region_radius);
}

Direction FeasibilityRestoration::compute_second_order_correction(Iterate& trial_iterate) {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   return this->subproblem->compute_second_order_correction(current_problem, trial_iterate);
}

PredictedReductionModel FeasibilityRestoration::generate_predicted_reduction_model(const Iterate& current_iterate, const Direction& direction) const {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   return this->subproblem->generate_predicted_reduction_model(current_problem, current_iterate, direction);
}

double FeasibilityRestoration::compute_infeasibility_measure(Iterate& iterate) {
   iterate.evaluate_constraints(this->original_problem);
   return this->original_problem.compute_constraint_violation(iterate.problem_evaluations.constraints, L1_NORM);
}

double FeasibilityRestoration::compute_optimality_measure(Iterate& iterate, const std::vector<size_t>& infeasible_constraints) {
   iterate.evaluate_constraints(this->original_problem);
   return this->original_problem.compute_constraint_violation(iterate.problem_evaluations.constraints, infeasible_constraints, L1_NORM);
}

void FeasibilityRestoration::register_accepted_iterate(Iterate& iterate) {
   const Problem& current_problem = (this->current_phase == OPTIMALITY) ? this->original_problem : this->relaxed_problem;
   this->subproblem->postprocess_accepted_iterate(current_problem, iterate);
}