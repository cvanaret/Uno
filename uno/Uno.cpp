// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianSubproblemSolverJointFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "../interfaces/C/Uno_C_API.h"
#include "linear_algebra/Vector.hpp"
#include "model/BoundRelaxedModel.hpp"
#include "model/Model.hpp"
#include "model/ScaledModel.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationStatus.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/InterruptHandler.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   // solve without user callbacks
   Result Uno::solve(const Model& model, const Options& options) {
      // pass user callbacks that do nothing
      NoUserCallbacks user_callbacks{};
      return this->solve(model, options, user_callbacks);
   }

   // solve with user callbacks
   Result Uno::solve(const Model& model, const Options& options, UserCallbacks& user_callbacks) {
      DISCRETE << "Model " << model.name << " has " << model.number_variables <<
         " variables, " << model.number_constraints << " constraints (" << model.get_equality_constraints().size() <<
         " equality, " << model.get_inequality_constraints().size() << " inequality)\n";
      DISCRETE << "Problem class: " << to_string(model.get_problem_type()) << '\n';

      // the functions are scaled at the initial point
      Vector<double> initial_primals(model.number_variables);
      model.initial_primal_point(initial_primals);

      // reformulate the model if:
      // - the user wants an interior-point method with log barrier function
      // - the model has bound constraints or inequality constraints
      if (options.get_string("inequality_handling_method") == "interior_point" && options.get_string("barrier_function") == "log" &&
            (model.has_bound_constraints() || model.has_inequality_constraints())) {
         // slightly relax the bound constraints. As in IPOPT, the original (unscaled) bounds are relaxed: relaxing the
         // scaled bounds would allow violations of the original constraints of up to relaxation/constraint scaling
         const BoundRelaxedModel bound_relaxed_model(model, options);
         const ScaledModel scaled_model(bound_relaxed_model, initial_primals, options);
         return uno_solve(model, scaled_model, options, user_callbacks);
      }
      const ScaledModel scaled_model(model, initial_primals, options);
      return uno_solve(model, scaled_model, options, user_callbacks);
   }

   // protected solve function
   // user_model is the original model provided by the user, model is the possibly reformulated model (function scaling,
   // relaxed bounds)
   Result Uno::uno_solve(const Model& user_model, const Model& model, const Options& options, UserCallbacks& user_callbacks) {
      InterruptGuard interrupt_guard;
      // drop the ingredients of a previous solve: they must not postprocess the iterate of this solve
      this->globalization_mechanism.reset();
      this->method_description = "strategy combination not initialized";
      Statistics statistics = create_statistics(model, options.get_bool("print_extended_statistics"));
      statistics.timers.wallclock.start();

      // initialize initial primal and dual points
      Iterate current_iterate(model.number_variables, model.number_constraints);
      model.initial_primal_point(current_iterate.primals);
      model.initial_dual_point(current_iterate.multipliers.constraints);
      model.reset_number_evaluations();
      EvaluationCache evaluation_cache{model};

      bool termination = false;
      OptimizationStatus optimization_status = OptimizationStatus::SUCCESS;

      // initialize the strategies and generate the initial iterate
      std::vector<OptionOverride> option_overrides;
      const size_t max_iterations = options.get_unsigned_int("max_iterations"); // maximum number of iterations
      if (max_iterations == 0) {
         termination = true;
         optimization_status = OptimizationStatus::ITERATION_LIMIT;
      }
      else {
         const std::optional<OptimizationStatus> initialization_failure = this->initialize(statistics, model,
            current_iterate, evaluation_cache, options, option_overrides);
         if (initialization_failure.has_value()) {
            termination = true;
            optimization_status = *initialization_failure;
         }
      }

      // outer loop: compute a sequence of accepted iterates
      size_t major_iterations = 0;
      // allocate the trial iterate once and for all
      Iterate trial_iterate(current_iterate);
      const double time_limit = options.get_double("time_limit"); // CPU time limit (can be inf)
      WarmstartInformation warmstart_information{};
      warmstart_information.whole_problem_changed();
      try {
         while (!termination) {
            ++major_iterations;
            warmstart_information.iteration = major_iterations;
            statistics.start_new_line();
            statistics.set("Iter", major_iterations);
            DEBUG << "\n### Outer iteration " << major_iterations << '\n';

            // compute an acceptable iterate by solving a subproblem at the current point
            warmstart_information.iterate_changed();
            this->globalization_mechanism->compute_next_iterate(statistics, model, current_iterate, trial_iterate,
               evaluation_cache, warmstart_information, user_callbacks);
            termination = check_termination(trial_iterate, major_iterations, max_iterations,
               statistics.timers.wallclock.get_elapsed_time(), time_limit, InterruptGuard::interrupt_requested(),
               optimization_status, user_callbacks);

            // the trial iterate becomes the current iterate for the next iteration
            std::swap(current_iterate, trial_iterate);
            std::swap(evaluation_cache.current_evaluations, evaluation_cache.trial_evaluations);
         }
      }
      // catch errors during the optimization process
      catch (std::exception& exception) {
         statistics.start_new_line();
         statistics.set("Status", exception.what());
         statistics.set("Time", Timer::format_to_seconds(statistics.timers.wallclock.get_elapsed_time()));
         if (Logger::level == INFO) statistics.print_current_line();
         DEBUG << exception.what() << '\n';
         optimization_status = OptimizationStatus::ALGORITHMIC_ERROR;
      }

      this->postprocess_solution(user_model, current_iterate, evaluation_cache.current_evaluations);
      statistics.timers.wallclock.stop();
      Result result = this->create_result(model, optimization_status, std::move(current_iterate),
         std::move(evaluation_cache.current_evaluations), major_iterations, statistics.timers, std::move(option_overrides));
      postprocess_multipliers_signs(model, result);
      this->print_optimization_summary(result, options.get_bool("print_solution"));
      return result;
   }
   
   std::string Uno::current_version() {
      return std::to_string(UNO_VERSION_MAJOR) + "." + std::to_string(UNO_VERSION_MINOR) + "." + std::to_string(UNO_VERSION_PATCH);
   }
   
   void Uno::print_available_strategies() {
      std::cout << "Available Uno strategies:\n";
      std::cout << "- Constraint relaxation strategies: " << join(ConstraintRelaxationStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- Globalization mechanisms: " << join(GlobalizationMechanismFactory::available_strategies, ", ") << '\n';
      std::cout << "- Globalization strategies: " << join(GlobalizationStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- Inequality handling methods: " << join(InequalityHandlingMethodFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Hessian models: " << join(HessianSubproblemSolverJointFactory::available_strategies, ", ") << '\n';
      std::cout << "- Inertia correction strategies: " << join(InertiaCorrectionStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- QP solvers: " << join(QPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- LP solvers: " << join(LPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- Linear solvers: " << join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- Presets: auto, filtersqp, filterslp, funnelsqp, ipopt\n";
   }

   const std::string& Uno::get_method_description() const {
      return this->method_description;
   }

   std::optional<OptimizationStatus> Uno::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) {
      try {
         statistics.start_new_line();
         statistics.set("Iter", 0);
         statistics.set("Status", "initial point");

         model.project_onto_variable_bounds(current_iterate.primals);

         // set the ingredients based on the user-defined options
         this->globalization_mechanism = GlobalizationMechanismFactory::create(model, options, option_overrides);
         if (this->globalization_mechanism == nullptr) {
            throw std::runtime_error("The ingredients could not be initialized");
         }
         this->globalization_mechanism->initialize(statistics, model, current_iterate, evaluation_cache, options,
            option_overrides);
         this->method_description = this->globalization_mechanism->get_name();

         statistics.set("Time", Timer::format_to_seconds(statistics.timers.wallclock.get_elapsed_time()));
         if (Logger::level == INFO) {
            statistics.print_header();
            statistics.print_current_line();
         }
         return std::nullopt;
      }
      // the functions could not be evaluated at the initial iterate
      catch (const EvaluationError& e) {
         DISCRETE << "An evaluation error occurred at the initial iterate: " << e.what() << '\n';
         this->discard_ingredients(model, current_iterate, evaluation_cache);
         return OptimizationStatus::EVALUATION_ERROR;
      }
      // the ingredients could not be created or initialized (invalid option, unavailable solver, ...)
      catch (const std::exception& e) {
         DISCRETE << "An error occurred during the initialization: " << e.what() << '\n';
         this->discard_ingredients(model, current_iterate, evaluation_cache);
         return OptimizationStatus::ALGORITHMIC_ERROR;
      }
   }

   // after a failed initialization, the ingredients may be partially constructed (e.g. a constraint relaxation strategy
   // without inequality handling method) and the iterate may have been reformulated (slacks, fixed variables, ...):
   // discard both, so that the result is the projected initial point of the model
   void Uno::discard_ingredients(const Model& model, Iterate& current_iterate, EvaluationCache& evaluation_cache) {
      this->globalization_mechanism.reset();
      this->method_description = "strategy combination not initialized";
      current_iterate = Iterate(model.number_variables, model.number_constraints);
      model.initial_primal_point(current_iterate.primals);
      model.project_onto_variable_bounds(current_iterate.primals);
      model.initial_dual_point(current_iterate.multipliers.constraints);
      evaluation_cache.current_evaluations.reset();
      evaluation_cache.trial_evaluations.reset();
   }

   Statistics Uno::create_statistics(const Model& model, bool print_extended_statistics) {
      Statistics statistics{print_extended_statistics};
      statistics.add_column("Iter", Statistics::int_width + 1, 3);
      statistics.add_column("||Step||", Statistics::double_width + 1, 2, /* is_extended = */ true);
      statistics.add_column("Objective", Statistics::double_width + 3, 4);
      if (model.is_constrained()) {
         statistics.add_column("Infeas", Statistics::double_width, 2);
      }
      statistics.add_column("Statio", Statistics::double_width + 1, 3);
      statistics.add_column("Compl", Statistics::double_width + 1, 3);
      statistics.add_column("Status", Statistics::string_width, 3);
      statistics.add_column("Time", Statistics::int_width + 2, 3);
      return statistics;
   }

   bool Uno::check_termination(const Iterate& trial_iterate, size_t iteration, size_t max_iterations, double current_time,
         double time_limit, bool user_interruption, OptimizationStatus& optimization_status, UserCallbacks& user_callbacks) {
      if (trial_iterate.status != SolutionStatus::NOT_OPTIMAL) {
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
      else if (user_callbacks.termination(trial_iterate.primals, trial_iterate.multipliers, trial_iterate.objective_multiplier,
            trial_iterate.scaled_residuals.primal_infeasibility, trial_iterate.scaled_residuals.stationarity,
            trial_iterate.scaled_residuals.complementarity)) {
         optimization_status = OptimizationStatus::USER_TERMINATION;
         return true;
      }
      else if (user_interruption) {
         optimization_status = OptimizationStatus::USER_INTERRUPTION;
         return true;
      }
      return false;
   }

   void Uno::postprocess_solution(const Model& user_model, Iterate& iterate, Evaluations& evaluations) const {
      // 1. primal-dual postprocessing only (unscaling of the multipliers, removal of slacks, fixed variables, ...)
      if (this->globalization_mechanism != nullptr) {
         this->globalization_mechanism->postprocess_solution(iterate, evaluations);
      }
      // 2. function values of the original model at the final primals. The postprocessing chain has either converted
      // the cached values to the units of the original model (ScaledModel) or invalidated them (e.g. snapped fixed
      // variables): only the invalidated values are evaluated, so that the user functions are not called again needlessly
      try {
         evaluations.evaluate_objective(user_model, iterate.primals);
      }
      catch (const std::exception& e) {
         DISCRETE  << "The objective could not be evaluated at the final iterate: " << e.what()  << '\n';
         evaluations.objective = Inf;
      }
      try {
         evaluations.evaluate_constraints(user_model, iterate.primals);
      }
      catch (const std::exception& e) {
         DISCRETE  << "The constraints could not be evaluated at the final iterate: " << e.what()  << '\n';
      }
      DEBUG2 << "Final iterate:\n" << iterate;
   }

   Result Uno::create_result(const Model& model, OptimizationStatus optimization_status, Iterate&& solution,
         Evaluations&& evaluations, size_t major_iterations, const Timers& timers, std::vector<OptionOverride>&& option_overrides) const {
      const size_t number_subproblems_solved = (this->globalization_mechanism != nullptr) ?
         this->globalization_mechanism->get_number_subproblems_solved() : 0;
      return {model.number_variables, model.number_constraints, model.base_indexing, optimization_status, solution.status,
         evaluations.objective, solution.scaled_residuals.primal_infeasibility, solution.scaled_residuals.stationarity,
         solution.scaled_residuals.complementarity, std::move(solution.primals), std::move(solution.multipliers.constraints),
         std::move(solution.multipliers.lower_bounds), std::move(solution.multipliers.upper_bounds),
         std::move(evaluations.constraints), major_iterations, timers, model.number_model_objective_evaluations(),
         model.number_model_constraints_evaluations(), model.number_model_objective_gradient_evaluations(),
         model.number_model_jacobian_evaluations(), model.number_model_hessian_evaluations(), number_subproblems_solved,
         std::move(option_overrides)};
   }

   // flip the signs of the multipliers, depending on what the sign convention of the Lagrangian is, and whether
   // we maximize
   void Uno::postprocess_multipliers_signs(const Model& model, Result& result) {
      if ((model.optimization_sense == 1. && model.lagrangian_sign_convention == -1.) ||
            (model.optimization_sense == -1. && model.lagrangian_sign_convention == 1.)) {
         // do nothing
      }
      else {
         // change the signs of the multipliers
         result.constraint_dual_solution.scale(-1.);
         result.lower_bound_dual_solution.scale(-1.);
         result.upper_bound_dual_solution.scale(-1.);
      }
      result.solution_objective *= model.optimization_sense;
   }

   void Uno::print_optimization_summary(const Result& result, bool print_solution) const {
      DISCRETE << "\nUno " << current_version() << " (" << this->method_description << ")\n";
      DISCRETE << Timer::get_current_date() << "\n\n";
      result.print(print_solution);
   }
} // namespace
