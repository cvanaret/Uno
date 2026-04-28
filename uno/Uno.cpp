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
#include "model/FixedBoundsConstraintsModel.hpp"
#include "model/HomogeneousEqualityConstrainedModel.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationCache.hpp"
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
   Result Uno::solve(const Model& model, Options& options) {
      // pass user callbacks that do nothing
      NoUserCallbacks user_callbacks{};
      return this->solve(model, options, user_callbacks);
   }

   // solve with user callbacks
   Result Uno::solve(const Model& model, Options& options, UserCallbacks& user_callbacks) {
      DISCRETE << "Original model " << model.name << '\n' << model.number_variables << " variables, " <<
         model.number_constraints << " constraints (" << model.get_equality_constraints().size() <<
         " equality, " << model.get_inequality_constraints().size() << " inequality)\n";
      INFO << "Problem type: " << to_string(model.get_problem_type()) << '\n';

      // reformulate the model if it is to be solved with an interior-point method with log barrier function
      if (options.get_string("inequality_handling_method") == "interior_point" && options.get_string("barrier_function") == "log") {
         // move the fixed variables to the set of general constraints
         const FixedBoundsConstraintsModel fixed_bound_model(model);
         // if an equality-constrained problem is required (e.g. interior points or AL), reformulate the model with slacks
         const HomogeneousEqualityConstrainedModel homogeneous_model(fixed_bound_model);
         // slightly relax the bound constraints
         const BoundRelaxedModel bound_relaxed_model(homogeneous_model, options);

         DISCRETE << "Reformulated model " << bound_relaxed_model.name << '\n' << bound_relaxed_model.number_variables << " variables, " <<
            bound_relaxed_model.number_constraints << " constraints (" << bound_relaxed_model.get_equality_constraints().size() <<
            " equality, " << bound_relaxed_model.get_inequality_constraints().size() << " inequality)\n";
         Result result = uno_solve(bound_relaxed_model, options, user_callbacks);
         // fix the dimensions
         result.number_variables = model.number_variables;
         result.number_constraints = model.number_constraints;
         return result;
      }
      else {
         return uno_solve(model, options, user_callbacks);
      }
   }

   // protected solve function
   Result Uno::uno_solve(const Model& model, Options& options, UserCallbacks& user_callbacks) {
      const Timer timer{};
      size_t major_iterations = 0;

      // initialize initial primal and dual points
      Iterate current_iterate(model.number_variables, model.number_constraints);
      model.initial_primal_point(current_iterate.primals);
      model.initial_dual_point(current_iterate.multipliers.constraints);
      model.reset_number_evaluations();
      EvaluationCache evaluation_cache{model};
      Statistics statistics = Uno::create_statistics(model);
      WarmstartInformation warmstart_information{};
      warmstart_information.whole_problem_changed();

      bool termination = false;
      OptimizationStatus optimization_status = OptimizationStatus::SUCCESS;

      // initialize the strategies and generate the initial iterate
      const bool initialization_success = this->initialize(statistics, model, current_iterate, options, evaluation_cache);
      if (!initialization_success) {
         termination = true;
         optimization_status = OptimizationStatus::EVALUATION_ERROR;
      }

      try {
         Iterate trial_iterate(current_iterate); // allocate the trial iterate once and for all here
         const size_t max_iterations = options.get_unsigned_int("max_iterations"); // maximum number of iterations
         const double time_limit = options.get_double("time_limit"); // CPU time limit (can be inf)
         if (max_iterations == 0) {
            termination = true;
            optimization_status = OptimizationStatus::ITERATION_LIMIT;
         }
         // outer loop: compute a sequence of accepted iterates
         while (!termination) {
            ++major_iterations;
            statistics.start_new_line();
            statistics.set("Major", major_iterations);
            DEBUG << "\n### Outer iteration " << major_iterations << '\n';

            // compute an acceptable iterate by solving a subproblem at the current point
            warmstart_information.iterate_changed();
            this->globalization_mechanism->compute_next_iterate(statistics, model, current_iterate, trial_iterate,
               evaluation_cache, warmstart_information, user_callbacks);
            const bool user_termination = user_callbacks.user_termination(trial_iterate.primals, trial_iterate.multipliers,
               trial_iterate.objective_multiplier, trial_iterate.progress.infeasibility, trial_iterate.residuals.stationarity,
               trial_iterate.residuals.complementarity);
            termination = Uno::check_termination(trial_iterate.status, major_iterations, max_iterations,
               timer.get_duration(), time_limit, user_termination, optimization_status);

            // the trial iterate becomes the current iterate for the next iteration
            std::swap(current_iterate, trial_iterate);
            std::swap(evaluation_cache.current_evaluations, evaluation_cache.trial_evaluations);
            evaluation_cache.trial_evaluations.reset();
         }
      }
      // catch errors during the optimization process
      catch (std::exception& exception) {
         statistics.start_new_line();
         statistics.set("Status", exception.what());
         if (Logger::level == INFO) statistics.print_current_line();
         DEBUG << exception.what() << '\n';
         optimization_status = OptimizationStatus::ALGORITHMIC_ERROR;
      }
      if (Logger::level == INFO) statistics.print_footer();

      Uno::postprocess_solution(model, current_iterate, evaluation_cache.current_evaluations);
      Result result = this->create_result(model, optimization_status, current_iterate, evaluation_cache.current_evaluations,
         major_iterations, timer);
      Uno::postprocess_multipliers_signs(model, result);
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
      std::cout << "- Presets: filtersqp, ipopt\n";
   }

   bool Uno::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, Options& options,
         EvaluationCache& evaluation_cache) {
      try {
         statistics.start_new_line();
         statistics.set("Major", 0);
         statistics.set("Status", "initial point");

         model.project_onto_variable_bounds(current_iterate.primals);

         // set the ingredients based on the user-defined options
         if (this->globalization_mechanism == nullptr) {
            // first call: create the mechanism (performs symbolic analysis of the linear system)
            this->globalization_mechanism = GlobalizationMechanismFactory::create(model, options);
            this->globalization_mechanism->initialize(statistics, model, current_iterate, evaluation_cache, options);
         }
         else {
            // subsequent call: reuse the mechanism (skips symbolic analysis)
            this->globalization_mechanism->reinitialize(statistics, model, current_iterate, evaluation_cache, options);
         }

         options.print_non_default();
         if (Logger::level == INFO) {
            statistics.print_header();
            statistics.print_current_line();
         }
         return true;
      }
      // catch errors at startup/initial iterate
      catch (const std::exception& e) {
         DISCRETE  << "An error occurred at the initial iterate: " << e.what()  << '\n';
         return false;
      }
   }

   Statistics Uno::create_statistics(const Model& model) {
      Statistics statistics{};
      statistics.add_column("Major", Statistics::int_width, 3, Statistics::column_order.at("Major"));
      statistics.add_column("||Step||", Statistics::double_width, 2, Statistics::column_order.at("||Step||"));
      statistics.add_column("Objective", Statistics::double_width + 1, 3, Statistics::column_order.at("Objective"));
      if (model.is_constrained()) {
         statistics.add_column("Infeas", Statistics::double_width, 2, Statistics::column_order.at("Infeas"));
      }
      statistics.add_column("Statio", Statistics::double_width, 2, Statistics::column_order.at("Statio"));
      statistics.add_column("Compl", Statistics::double_width, 2, Statistics::column_order.at("Compl"));
      statistics.add_column("Status", Statistics::string_width, 3, Statistics::column_order.at("Status"));
      return statistics;
   }

   bool Uno::check_termination(SolutionStatus solution_status, size_t iteration, size_t max_iterations, double current_time,
         double time_limit, bool user_termination, OptimizationStatus& optimization_status) {
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
      else if (user_termination) {
         optimization_status = OptimizationStatus::USER_TERMINATION;
         return true;
      }
      return false;
   }

   void Uno::postprocess_solution(const Model& model, Iterate& iterate, Evaluations& evaluations) {
      try {
         // in case the objective was not yet evaluated, evaluate it
         evaluations.evaluate_objective(model, iterate.primals);
      }
      catch (const std::exception& e) {
         DISCRETE  << "The objective could not be evaluated at the final iterate: " << e.what()  << '\n';
         evaluations.objective = INF<double>;
      }
      model.postprocess_solution(iterate);
      DEBUG2 << "Final iterate:\n" << iterate;
   }

   Result Uno::create_result(const Model& model, OptimizationStatus optimization_status, const Iterate& solution,
         const Evaluations& evaluations, size_t major_iterations, const Timer& timer) const {
      const size_t number_subproblems_solved = (this->globalization_mechanism != nullptr) ?
         this->globalization_mechanism->get_number_subproblems_solved() : 0;
      return {model.number_variables, model.number_constraints, model.base_indexing, optimization_status, solution.status,
         evaluations.objective, solution.progress.infeasibility, solution.residuals.stationarity,
         solution.residuals.complementarity, solution.primals, solution.multipliers.constraints,
         solution.multipliers.lower_bounds, solution.multipliers.upper_bounds, major_iterations, timer.get_duration(),
         model.number_model_objective_evaluations(), model.number_model_constraints_evaluations(),
         model.number_model_objective_gradient_evaluations(), model.number_model_jacobian_evaluations(),
         model.number_model_hessian_evaluations(), number_subproblems_solved};
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

   std::string Uno::get_strategy_combination() const {
      return (this->globalization_mechanism != nullptr) ? this->globalization_mechanism->get_name() :
         "strategy combination not initialized";
   }

   void Uno::print_optimization_summary(const Result& result, bool print_solution) const {
      DISCRETE << "\nUno " << Uno::current_version() << " (" << this->get_strategy_combination() << ")\n";
      DISCRETE << Timer::get_current_date();
      DISCRETE << "────────────────────────────────────────\n";
      result.print(print_solution);
   }
} // namespace