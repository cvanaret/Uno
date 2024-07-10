// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      max_iterations(options.get_unsigned_int("max_iterations")),
      time_limit(options.get_double("time_limit")) {
}

Result Uno::solve(const Model& model, Iterate& current_iterate, const Options& options) {
   std::cout << "Problem " << model.name << '\n' << model.number_variables << " variables, " << model.number_constraints << " constraints\n\n";
   
   Timer timer{};
   Statistics statistics = Uno::create_statistics(model, options);
   // use the initial primal-dual point to initialize the strategies and generate the initial iterate
   this->initialize(statistics, current_iterate, options);

   bool termination = false;
   current_iterate.status = TerminationStatus::NOT_OPTIMAL;
   // allocate the trial iterate once and for all here
   Iterate trial_iterate(current_iterate.number_variables, current_iterate.number_constraints);
   size_t major_iterations = 0;
   try {
      // check for termination
      while (not termination) {
         major_iterations++;
         statistics.start_new_line();
         statistics.set("iter", major_iterations);
         DEBUG << "### Outer iteration " << major_iterations << '\n';

         // compute an acceptable iterate by solving a subproblem at the current point
         this->globalization_mechanism.compute_next_iterate(statistics, model, current_iterate, trial_iterate);
         // determine if Uno can terminate
         termination = this->termination_criteria(trial_iterate.status, major_iterations, timer.get_duration());
         // the trial iterate becomes the current iterate for the next iteration
         std::swap(current_iterate, trial_iterate);
      }
   }
   catch (std::exception& exception) {
      statistics.start_new_line();
      statistics.set("status", exception.what());
      if (Logger::level == INFO) statistics.print_current_line();
      DEBUG << exception.what();
   }
   if (Logger::level == INFO) statistics.print_footer();

   Uno::postprocess_iterate(model, current_iterate, current_iterate.status);
   return this->create_result(model, current_iterate, major_iterations, timer);
}

void Uno::initialize(Statistics& statistics, Iterate& current_iterate, const Options& options) {
   try {
      statistics.start_new_line();
      statistics.set("iter", 0);
      statistics.set("status", "initial point");
      this->globalization_mechanism.initialize(statistics, current_iterate, options);
      options.print(true);
      if (Logger::level == INFO) statistics.print_current_line();
   }
   catch (const std::exception& e) {
      ERROR << RED << "An error occurred at the initial iterate: " << e.what() << RESET;
      throw;
   }
}

Statistics Uno::create_statistics(const Model& model, const Options& options) {
   Statistics statistics(options);
   statistics.add_column("iter", Statistics::int_width, options.get_int("statistics_major_column_order"));
   statistics.add_column("step norm", Statistics::double_width - 4, options.get_int("statistics_step_norm_column_order"));
   statistics.add_column("objective", Statistics::double_width - 4, options.get_int("statistics_objective_column_order"));
   if (model.is_constrained()) {
      statistics.add_column("primal feas.", Statistics::double_width - 3, options.get_int("statistics_primal_feasibility_column_order"));
   }
   statistics.add_column("stationarity", Statistics::double_width - 3, options.get_int("statistics_stationarity_column_order"));
   statistics.add_column("dual feas.", Statistics::double_width - 4, options.get_int("statistics_dual_feasibility_column_order"));
   statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
   statistics.add_column("status", Statistics::string_width, options.get_int("statistics_status_column_order"));
   return statistics;
}

bool Uno::termination_criteria(TerminationStatus current_status, size_t iteration, double current_time) const {
   return current_status != TerminationStatus::NOT_OPTIMAL || this->max_iterations <= iteration || this->time_limit <= current_time;
}

void Uno::postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status) {
   // in case the objective was not yet evaluated, evaluate it
   iterate.evaluate_objective(model);
   model.postprocess_solution(iterate, termination_status);
   DEBUG2 << "Final iterate:\n" << iterate;
}

Result Uno::create_result(const Model& model, Iterate& current_iterate, size_t major_iterations, const Timer& timer) {
   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t number_hessian_evaluations = this->globalization_mechanism.get_hessian_evaluation_count();
   return {std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations, timer.get_duration(),
         Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_objective_gradient,
         Iterate::number_eval_jacobian, number_hessian_evaluations, number_subproblems_solved};
}

void join(const std::vector<std::string>& vector, char separator) {
   if (not vector.empty()) {
      std::cout << vector[0];
      for (size_t variable_index: Range(1, vector.size())) {
         std::cout << separator << ' ' << vector[variable_index];
      }
   }
}

void Uno::print_uno_version() {
   std::cout << "Welcome in Uno 1.0\n";
   std::cout << "To solve an AMPL model, type ./uno_ampl path_to_file/file.nl\n";
   std::cout << "To choose a constraint relaxation strategy, use the argument -constraint_relaxation_strategy "
                "[feasibility_restoration|l1_relaxation]\n";
   std::cout << "To choose a subproblem method, use the argument -subproblem [QP|LP|primal_dual_interior_point]\n";
   std::cout << "To choose a globalization mechanism, use the argument -globalization_mechanism [LS|TR]\n";
   std::cout << "To choose a globalization strategy, use the argument -globalization_strategy "
                "[l1_merit|fletcher_filter_method|waechter_filter_method]\n";
   std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|byrd]\n";
   std::cout << "The options can be combined in the same command line. Autocompletion is possible (see README).\n";
}

void Uno::print_available_strategies() {
   std::cout << "Available strategies:\n";
   std::cout << "Constraint relaxation strategies: ";
   join(ConstraintRelaxationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization mechanisms: ";
   join(GlobalizationMechanismFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization strategies: ";
   join(GlobalizationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Subproblems: ";
   join(SubproblemFactory::available_strategies(), ',');
   std::cout << '\n';
}

void Uno::print_strategy_combination(const Options& options) {
   std::string combination = options.get_string("globalization_mechanism") + " " + options.get_string("constraint_relaxation_strategy") + " " +
                             options.get_string("globalization_strategy") + " " + options.get_string("subproblem");
   std::cout << "\nUno (" << combination << ")\n";
}

void Uno::print_optimization_summary(const Options& options, const Result& result) {
   Uno::print_strategy_combination(options);
   std::cout << Timer::get_current_date();
   std::cout << "────────────────────────────────────────\n";
   const bool print_solution = options.get_bool("print_solution");
   result.print(print_solution);
}