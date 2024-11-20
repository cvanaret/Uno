// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblems/SubproblemFactory.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "solvers/QPSolverFactory.hpp"
#include "solvers/LPSolverFactory.hpp"
#include "solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

namespace uno {
   Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
         globalization_mechanism(globalization_mechanism),
         max_iterations(options.get_unsigned_int("max_iterations")),
         time_limit(options.get_double("time_limit")),
         print_solution(options.get_bool("print_solution")),
         strategy_combination(Uno::get_strategy_combination(options)) { }
   
   Level Logger::level = INFO;

   std::optional<Result> Uno::solve(const Model& model, Iterate& current_iterate, const Options& options) {
      Timer timer{};
      Statistics statistics = Uno::create_statistics(model, options);
      WarmstartInformation warmstart_information{};
      warmstart_information.whole_problem_changed();

      try {
         // use the initial primal-dual point to initialize the strategies and generate the initial iterate
         this->initialize(statistics, current_iterate, options);
         // allocate the trial iterate once and for all here
         Iterate trial_iterate(current_iterate);

         size_t major_iterations = 0;
         try {
            bool termination = false;
            // check for termination
            while (not termination) {
               major_iterations++;
               statistics.start_new_line();
               statistics.set("iter", major_iterations);
               DEBUG << "### Outer iteration " << major_iterations << '\n';

               // compute an acceptable iterate by solving a subproblem at the current point
               warmstart_information.iterate_changed();
               this->globalization_mechanism.compute_next_iterate(statistics, model, current_iterate, trial_iterate, warmstart_information);
               termination = this->termination_criteria(trial_iterate.status, major_iterations, timer.get_duration());
               // the trial iterate becomes the current iterate for the next iteration
               std::swap(current_iterate, trial_iterate);
            }
         }
         catch (std::exception& exception) {
            statistics.start_new_line();
            statistics.set("status", exception.what());
            if (Logger::level == INFO) statistics.print_current_line();
            DEBUG << exception.what() << '\n';
         }
         if (Logger::level == INFO) statistics.print_footer();

         Uno::postprocess_iterate(model, current_iterate, current_iterate.status);
         Result result = this->create_result(model, current_iterate, major_iterations, timer);
         this->print_optimization_summary(result);
         return result;
      }
      catch (const std::exception& e) {
         DISCRETE  << "An error occurred at the initial iterate: " << e.what()  << '\n';
         return std::nullopt;
      }
   }

   void Uno::initialize(Statistics& statistics, Iterate& current_iterate, const Options& options) {
      statistics.start_new_line();
      statistics.set("iter", 0);
      statistics.set("status", "initial point");
      this->globalization_mechanism.initialize(statistics, current_iterate, options);
      options.print_used();
      if (Logger::level == INFO) statistics.print_current_line();
      current_iterate.status = TerminationStatus::NOT_OPTIMAL;
   }

   Statistics Uno::create_statistics(const Model& model, const Options& options) {
      Statistics statistics(options);
      statistics.add_column("iter", Statistics::int_width, options.get_int("statistics_major_column_order"));
      statistics.add_column("step norm", Statistics::double_width - 5, options.get_int("statistics_step_norm_column_order"));
      statistics.add_column("objective", Statistics::double_width - 5, options.get_int("statistics_objective_column_order"));
      if (model.is_constrained()) {
         statistics.add_column("primal feas", Statistics::double_width - 4, options.get_int("statistics_primal_feasibility_column_order"));
      }
      statistics.add_column("stationarity", Statistics::double_width - 3, options.get_int("statistics_stationarity_column_order"));
      statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
      statistics.add_column("status", Statistics::string_width - 9, options.get_int("statistics_status_column_order"));
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

   std::string Uno::current_version() {
      return "1.3.0";
   }

   void Uno::print_available_strategies() {
      std::cout << "Available strategies:\n";
      std::cout << "- Constraint relaxation strategies: " << join(ConstraintRelaxationStrategyFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Globalization mechanisms: " << join(GlobalizationMechanismFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Globalization strategies: " << join(GlobalizationStrategyFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Subproblems: " << join(SubproblemFactory::available_strategies(), ", ") << '\n';
      std::cout << "- QP solvers: " << join(QPSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- LP solvers: " << join(LPSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- Linear solvers: " << join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", ") << '\n';
   }

   std::string Uno::get_strategy_combination(const Options& options) {
      return options.get_string("globalization_mechanism") + " " + options.get_string("constraint_relaxation_strategy") + " " +
                                options.get_string("globalization_strategy") + " " + options.get_string("subproblem");

   }

   void Uno::print_optimization_summary(const Result& result) {
      DISCRETE << "\nUno " << Uno::current_version() << " (" << this->strategy_combination << ")\n";
      DISCRETE << Timer::get_current_date();
      DISCRETE << "────────────────────────────────────────\n";
      result.print(this->print_solution);
   }
} // namespace
