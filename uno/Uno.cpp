// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/BoundRelaxedModel.hpp"
#include "model/FixedBoundsConstraintsModel.hpp"
#include "model/HomogeneousEqualityConstrainedModel.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "optimization/OptimizationStatus.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   Level Logger::level = INFO;

   // solve without user callbacks
   Result Uno::solve(const Model& model, const Options& options) {
      // pass user callbacks that do nothing
      NoUserCallbacks user_callbacks{};
      return this->solve(model, options, user_callbacks);
   }

   // solve with user callbacks
   Result Uno::solve(const Model& model, const Options& options, UserCallbacks& user_callbacks) {
      DISCRETE << "Original model " << model.name << '\n' << model.number_variables << " variables, " <<
         model.number_constraints << " constraints (" << model.get_equality_constraints().size() <<
         " equality, " << model.get_inequality_constraints().size() << " inequality)\n";

      // reformulate the model if it is to be solved with an interior-point method
      if (options.get_string("inequality_handling_method") == "primal_dual_interior_point") {
         // move the fixed variables to the set of general constraints
         const FixedBoundsConstraintsModel fixed_bound_model(model);
         // if an equality-constrained problem is required (e.g. interior points or AL), reformulate the model with slacks
         const HomogeneousEqualityConstrainedModel homogeneous_model(fixed_bound_model);
         // slightly relax the bound constraints
         const BoundRelaxedModel bound_relaxed_model(homogeneous_model, options);

         DISCRETE << "Reformulated model " << bound_relaxed_model.name << '\n' << bound_relaxed_model.number_variables << " variables, " <<
            bound_relaxed_model.number_constraints << " constraints (" << bound_relaxed_model.get_equality_constraints().size() <<
            " equality, " << bound_relaxed_model.get_inequality_constraints().size() << " inequality)\n";
         return uno_solve(bound_relaxed_model, options, user_callbacks);
      }
      else {
         return uno_solve(model, options, user_callbacks);
      }
   }

   // protected solve function
   Result Uno::uno_solve(const Model& model, const Options& options, UserCallbacks& user_callbacks) {
      const Timer timer{};
      // pick the ingredients based on the user-defined options
      Uno::pick_ingredients(model, options);
      Statistics statistics = Uno::create_statistics(model, options);
      WarmstartInformation warmstart_information{};
      warmstart_information.whole_problem_changed();

      // initialize initial primal and dual points
      Iterate current_iterate(model.number_variables, model.number_constraints);
      model.initial_primal_point(current_iterate.primals);
      model.initial_dual_point(current_iterate.multipliers.constraints);

      size_t major_iterations = 0;
      OptimizationStatus optimization_status = OptimizationStatus::SUCCESS;
      const size_t max_iterations = options.get_unsigned_int("max_iterations"); // maximum number of iterations
      const double time_limit = options.get_double("time_limit"); // CPU time limit (can be inf)
      try {
         // use the initial primal-dual point to initialize the strategies and generate the initial iterate
         this->initialize(statistics, model, current_iterate, options);
         // allocate the trial iterate once and for all here
         Iterate trial_iterate(current_iterate);

         try {
            bool termination = false;
            // check for termination
            while (!termination) {
               ++major_iterations;
               statistics.start_new_line();
               statistics.set("iter", major_iterations);
               DEBUG << "\n### Outer iteration " << major_iterations << '\n';

               // compute an acceptable iterate by solving a subproblem at the current point
               warmstart_information.iterate_changed();
               this->globalization_mechanism->compute_next_iterate(statistics, *this->constraint_relaxation_strategy,
                  *this->globalization_strategy, model, current_iterate, trial_iterate, this->direction, warmstart_information,
                  user_callbacks);
               GlobalizationMechanism::set_dual_residuals_statistics(statistics, trial_iterate);
               termination = Uno::termination_criteria(trial_iterate.status, major_iterations, max_iterations,
                  timer.get_duration(), time_limit, optimization_status);
               user_callbacks.notify_new_primals(trial_iterate.primals);
               user_callbacks.notify_new_multipliers(trial_iterate.multipliers);

               // the trial iterate becomes the current iterate for the next iteration
               std::swap(current_iterate, trial_iterate);
            }
         }
         catch (std::exception& exception) {
            statistics.start_new_line();
            statistics.set("status", exception.what());
            if (Logger::level == INFO) statistics.print_current_line();
            DEBUG << exception.what() << '\n';
            optimization_status = OptimizationStatus::ALGORITHMIC_ERROR;
         }
         if (Logger::level == INFO) statistics.print_footer();

         Uno::postprocess_iterate(model, current_iterate);
      }
      catch (const std::exception& e) {
         DISCRETE  << "An error occurred at the initial iterate: " << e.what()  << '\n';
         optimization_status = OptimizationStatus::EVALUATION_ERROR;
      }
      Result result = this->create_result(model, optimization_status, current_iterate, major_iterations, timer);
      this->print_optimization_summary(result, options.get_bool("print_solution"));
      return result;
   }
   
   std::string Uno::current_version() {
      return "2.2.0";
   }
   
   void Uno::print_available_strategies() {
      std::cout << "Available Uno strategies:\n";
      std::cout << "- Constraint relaxation strategies: " << join(ConstraintRelaxationStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- Globalization mechanisms: " << join(GlobalizationMechanismFactory::available_strategies, ", ") << '\n';
      std::cout << "- Globalization strategies: " << join(GlobalizationStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- Inequality handling methods: " << join(InequalityHandlingMethodFactory::available_strategies(), ", ") << '\n';
      std::cout << "- QP solvers: " << join(QPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- LP solvers: " << join(LPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- Linear solvers: " << join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- Presets: filtersqp, ipopt\n";
   }

   void Uno::pick_ingredients(const Model& model, const Options& options) {
      const bool unconstrained_model = (model.number_constraints == 0);
      this->constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(unconstrained_model, options);
      this->globalization_strategy = GlobalizationStrategyFactory::create(unconstrained_model, options);
      this->globalization_mechanism = GlobalizationMechanismFactory::create(options);
   }

   void Uno::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, const Options& options) {
      statistics.start_new_line();
      statistics.set("iter", 0);
      statistics.set("status", "initial point");

      model.project_onto_variable_bounds(current_iterate.primals);
      // TODO here we don't know if there's a trust-region radius yet!
      this->constraint_relaxation_strategy->initialize(statistics, model, current_iterate, this->direction, INF<double>, options);
      GlobalizationMechanism::set_primal_statistics(statistics, model, current_iterate);
      GlobalizationMechanism::set_dual_residuals_statistics(statistics, current_iterate);
      this->globalization_strategy->initialize(statistics, current_iterate, options);
      this->globalization_mechanism->initialize(statistics, options);

      options.print_used_overwritten();
      if (Logger::level == INFO) {
         statistics.print_header();
         statistics.print_current_line();
      }
      current_iterate.status = SolutionStatus::NOT_OPTIMAL;
   }

   Statistics Uno::create_statistics(const Model& model, const Options& options) {
      Statistics statistics{};
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

   bool Uno::termination_criteria(SolutionStatus solution_status, size_t iteration, size_t max_iterations, double current_time,
         double time_limit, OptimizationStatus& optimization_status) {
      if (solution_status != SolutionStatus::NOT_OPTIMAL) {
         return true;
      }
      else if (max_iterations <= iteration) {
         optimization_status = OptimizationStatus::ITERATION_LIMIT;
         return true;
      }
      else if (time_limit <= current_time) {
         optimization_status = OptimizationStatus::TIME_LIMIT;
         return true;
      }
      return false;
   }

   void Uno::postprocess_iterate(const Model& model, Iterate& iterate) {
      // in case the objective was not yet evaluated, evaluate it
      iterate.evaluate_objective(model);
      model.postprocess_solution(iterate);
      DEBUG2 << "Final iterate:\n" << iterate;
   }

   Result Uno::create_result(const Model& model, OptimizationStatus optimization_status, Iterate& solution, size_t major_iterations,
         const Timer& timer) const {
      const size_t number_subproblems_solved = this->constraint_relaxation_strategy->get_number_subproblems_solved();
      const size_t number_hessian_evaluations = this->constraint_relaxation_strategy->get_hessian_evaluation_count();
      return {model.number_variables, model.number_constraints, optimization_status, solution.status,
         solution.evaluations.objective, solution.progress.infeasibility, solution.residuals.stationarity,
         solution.residuals.complementarity, solution.primals, solution.multipliers.constraints,
         solution.multipliers.lower_bounds, solution.multipliers.upper_bounds, major_iterations, timer.get_duration(),
         Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_objective_gradient,
         Iterate::number_eval_jacobian, number_hessian_evaluations, number_subproblems_solved};
   }

   std::string Uno::get_strategy_combination() const {
      return this->globalization_mechanism->get_name() + " " + this->globalization_strategy->get_name() + " " +
         this->constraint_relaxation_strategy->get_name();
   }

   void Uno::print_optimization_summary(const Result& result, bool print_solution) const {
      DISCRETE << "\nUno " << Uno::current_version() << " (" << this->get_strategy_combination() << ")\n";
      DISCRETE << Timer::get_current_date();
      DISCRETE << "────────────────────────────────────────\n";
      result.print(print_solution);
   }
} // namespace