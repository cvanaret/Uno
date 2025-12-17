// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
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
         return uno_solve(bound_relaxed_model, options, user_callbacks);
      }
      else {
         return uno_solve(model, options, user_callbacks);
      }
   }

   // protected solve function
   Result Uno::uno_solve(const Model& model, const Options& options, UserCallbacks& user_callbacks) {
      const Timer timer{};
      model.reset_number_evaluations();
      // pick the ingredients based on the user-defined options
      Uno::pick_ingredients(model, options);
      Statistics statistics = Uno::create_statistics(model);
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
               statistics.set("Major", major_iterations);
               DEBUG << "\n### Outer iteration " << major_iterations << '\n';

               // compute an acceptable iterate by solving a subproblem at the current point
               warmstart_information.iterate_changed();
               this->globalization_mechanism->compute_next_iterate(statistics, model, current_iterate, trial_iterate,
                  this->direction, warmstart_information, user_callbacks);
               GlobalizationMechanism::set_dual_residuals_statistics(statistics, trial_iterate);
               const bool user_termination = user_callbacks.user_termination(trial_iterate.primals, trial_iterate.multipliers,
                  trial_iterate.objective_multiplier, trial_iterate.progress.infeasibility, trial_iterate.residuals.stationarity,
                  trial_iterate.residuals.complementarity);
               termination = Uno::termination_criteria(trial_iterate.status, major_iterations, max_iterations,
                  timer.get_duration(), time_limit, user_termination, optimization_status);
               
               // the trial iterate becomes the current iterate for the next iteration
               std::swap(current_iterate, trial_iterate);
            }
         }
         catch (std::exception& exception) {
            statistics.start_new_line();
            statistics.set("Status", exception.what());
            if (Logger::level == INFO) statistics.print_current_line();
            DEBUG << exception.what() << '\n';
            optimization_status = OptimizationStatus::ALGORITHMIC_ERROR;
         }
         if (Logger::level == INFO) statistics.print_footer();

         Uno::postprocess_solution(model, current_iterate);
      }
      catch (const std::exception& e) {
         DISCRETE  << "An error occurred at the initial iterate: " << e.what()  << '\n';
         optimization_status = OptimizationStatus::EVALUATION_ERROR;
      }
      Result result = this->create_result(model, optimization_status, current_iterate, major_iterations, timer);
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
      std::cout << "- Hessian models: " << join(HessianModelFactory::available_strategies, ", ") << '\n';
      std::cout << "- Inertia correction strategies: " << join(InertiaCorrectionStrategyFactory::available_strategies, ", ") << '\n';
      std::cout << "- QP solvers: " << join(QPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- LP solvers: " << join(LPSolverFactory::available_solvers, ", ") << '\n';
      std::cout << "- Linear solvers: " << join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- Presets: filtersqp, ipopt\n";
   }

   void Uno::pick_ingredients(const Model& model, const Options& options) {
      this->globalization_mechanism = GlobalizationMechanismFactory::create(model, options);
   }

   void Uno::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, const Options& options) {
      statistics.start_new_line();
      statistics.set("Major", 0);
      statistics.set("Status", "initial point");

      model.project_onto_variable_bounds(current_iterate.primals);
      this->globalization_mechanism->initialize(statistics, model, current_iterate, this->direction);
      GlobalizationMechanism::set_primal_statistics(statistics, model, current_iterate);
      GlobalizationMechanism::set_dual_residuals_statistics(statistics, current_iterate);

      options.print_non_default();
      if (Logger::level == INFO) {
         statistics.print_header();
         statistics.print_current_line();
      }
      current_iterate.status = SolutionStatus::NOT_OPTIMAL;
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

   bool Uno::termination_criteria(SolutionStatus solution_status, size_t iteration, size_t max_iterations, double current_time,
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

   void Uno::postprocess_solution(const Model& model, Iterate& iterate) {
      // in case the objective was not yet evaluated, evaluate it
      iterate.evaluate_objective(model);
      model.postprocess_solution(iterate);
      DEBUG2 << "Final iterate:\n" << iterate;
   }

   Result Uno::create_result(const Model& model, OptimizationStatus optimization_status, const Iterate& solution,
         size_t major_iterations, const Timer& timer) const {
      const size_t number_subproblems_solved = this->globalization_mechanism->get_number_subproblems_solved();
      //const size_t number_hessian_evaluations = this->constraint_relaxation_strategy->get_hessian_evaluation_count();
      return {model.number_variables, model.number_constraints, optimization_status, solution.status,
         solution.evaluations.objective, solution.progress.infeasibility, solution.residuals.stationarity,
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
      return this->globalization_mechanism->get_name();
   }

   void Uno::print_optimization_summary(const Result& result, bool print_solution) const {
      DISCRETE << "\nUno " << Uno::current_version() << " (" << this->get_strategy_combination() << ")\n";
      DISCRETE << Timer::get_current_date();
      DISCRETE << "────────────────────────────────────────\n";
      result.print(print_solution);
   }
} // namespace