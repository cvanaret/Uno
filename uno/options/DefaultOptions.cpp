// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DefaultOptions.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   void DefaultOptions::load(Options& options) {
      DefaultOptions::determine_subproblem_solvers(options);

      // preset
      options.set_string("preset", "auto");

      /** termination **/
      // primal tolerance (constraint violation)
      options.set_double("primal_tolerance", 1e-8);
      // dual tolerance (stationarity and complementarity)
      options.set_double("dual_tolerance", 1e-8);
      // bound on the unscaled stationarity residual required for optimal termination
      options.set_double("tight_unscaled_stationarity_tolerance", 1.);
      // bound on the unscaled constraint violation required for optimal termination
      options.set_double("tight_unscaled_primal_infeasibility_tolerance", 1e-4);
      // bound on the unscaled complementarity residual
      options.set_double("tight_unscaled_complementarity_tolerance", 1e-4);
      // tolerance on the scaled primal infeasibility, stationarity and complementarity for a solution to be acceptable
      options.set_double("acceptable_tolerance", 1e-6);
      // bound on the unscaled stationarity residual for a solution to be acceptable
      options.set_double("acceptable_unscaled_stationarity_tolerance", 1e10);
      // bound on the unscaled constraint violation for a solution to be acceptable
      options.set_double("acceptable_unscaled_primal_infeasibility_tolerance", 1e-2);
      // bound on the unscaled complementarity residual for a solution to be acceptable
      options.set_double("acceptable_unscaled_complementarity_tolerance", 1e-2);
      // bound on the relative change of the objective between consecutive iterates for a solution to be acceptable
      options.set_double("acceptable_objective_change_tolerance", 1e20);
      // number of consecutive acceptable iterates until acceptable termination
      options.set_integer("acceptable_iteration_threshold", 15);

      // maximum number of outer iterations
      options.set_integer("max_iterations", 2000);
      // CPU time limit (in seconds)
      options.set_double("time_limit", Inf);
      // print optimal solution (true|false)
      options.set_bool("print_solution", false);
      // threshold on norm of primal variables to declare diverging iterate
      options.set_double("diverging_iterate_threshold", 1e20);
      // threshold on objective to declare unbounded NLP
      options.set_double("unbounded_objective_threshold", -1e20);

      /** main options **/
      options.set_bool("write_solution_to_file", false);
      // logging level (SILENT|DISCRETE|WARNING|INFO|DEBUG|DEBUG2|DEBUG3)
      options.set_string("logger", "INFO");
      // norm of the progress measures (L1|L2|INF)
      options.set_string("progress_norm", "L1");
      // norm of the primal-dual residuals (L1|L2|INF)
      options.set_string("residual_norm", "INF");
      options.set_double("residual_scaling_threshold", 100.);
      options.set_bool("protect_actual_reduction_against_roundoff", false);
      options.set_double("protected_actual_reduction_macheps_coefficient", 10.);
      options.set_bool("print_subproblem", false);
      // use function scaling based on the initial gradients (true|false)
      options.set_bool("use_function_scaling", false);
      // threshold on the magnitudes of the functions' gradients
      options.set_double("function_scaling_threshold", 100.);
      // lower bound on the function scaling coefficients
      options.set_double("function_scaling_lower_bound", 1e-8);
      // whether minor iterations (line search and trust-region method) are printed in the statistics table
      options.set_bool("print_minor_iterations", false);
      // whether the statistics table is printed in extended form
      options.set_bool("print_extended_statistics", false);
      // dump KKT systems
      //options.set_string("dump_kkt_path", "");
      options.set_integer("dump_kkt_frequency", 1);

      /** globalization strategy options **/
      options.set_double("armijo_decrease_fraction", 1e-4);
      options.set_double("armijo_tolerance", 1e-9);

      /** switching method options **/
      options.set_double("switching_delta", 0.999);
      options.set_double("switching_merit_exponent", 1.);
      options.set_double("switching_infeasibility_exponent", 2.);

      /** merit function options **/
      options.set_double("sufficient_infeasibility_decrease_ratio", 0.9);

      /** filter method options **/
      options.set_double("filter_beta", 0.999);
      options.set_double("filter_gamma", 0.001);
      options.set_double("filter_ubd", 1e2);
      options.set_double("filter_fact", 1.25);
      options.set_integer("filter_capacity", 50); // 0: unlimited
      // used by Waechter filter method
      options.set_double("filter_sufficient_infeasibility_decrease_factor", 0.9);
      options.set_integer("filter_reset_iteration_threshold", 5); // >= 1
      options.set_integer("max_number_filter_resets", 5); // >= 0
      options.set_integer("filter_ubd_raise_trigger", 0); // 0: disabled
      options.set_double("filter_ubd_raise_factor", 100.); // > 1
      options.set_integer("filter_ubd_max_raises", 4);
      options.set_double("filter_merit_max_increase", 5.); // orders of magnitude (IPOPT obj_max_inc)
      options.set_bool("filter_margin_uses_entry_infeasibility", false);

      /** funnel options **/
      options.set_double("funnel_kappa", 0.5);
      options.set_double("funnel_beta", 0.9999);
      options.set_double("funnel_gamma", 0.001);
      options.set_double("funnel_ubd", 1.0);
      options.set_double("funnel_fact", 1.5);
      options.set_integer("funnel_update_strategy", 1);
      options.set_bool("funnel_require_acceptance_wrt_current_iterate", false);

      /** line search options */
      // backtracking ratio
      options.set_double("LS_backtracking_ratio", 0.5);
      // use the primal-dual and dual step lengths to scale the dual directions when assembling the trial iterate
      options.set_bool("LS_scale_duals_with_step_length", true);
      // safety factor for the minimum step length
      options.set_double("gamma_alpha", 0.05);
      // watchdog: number of consecutive shortened steps that triggers the watchdog (0 disables it)
      options.set_integer("LS_watchdog_shortened_iterations_trigger", 10);
      // watchdog: maximum number of iterations in which the full step is accepted without globalization test
      options.set_integer("LS_watchdog_max_trial_iterations", 3);

      /** Second-Order Corrections (SOC) options */
      // maximum number of SOC iterations
      options.set_integer("SOC_max_iterations", 4);
      // fraction of the infeasibility in the SOC termination criterion
      options.set_double("SOC_infeasibility_fraction", 0.99);

      /* quasi-Newton options */
      options.set_integer("quasi_newton_memory_size", 6);
      options.set_double("quasi_newton_delta_lower_bound", 1e-8);
      options.set_double("quasi_newton_delta_upper_bound", 1e8);
      options.set_integer("LBFGS_max_skips_before_reset", 2);
      options.set_double("LSR1_pivot_max_magnitude", 1e-7);

      /** regularization options **/
      // regularization failure threshold
      options.set_double("regularization_failure_threshold", 1e20); // https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_max_hessian_perturbation
      // Hessian regularization: initial value
      options.set_double("primal_regularization_initial_factor", 1e-4);
      options.set_double("regularization_increase_factor", 2);
      // regularization of augmented system
      options.set_double("dual_regularization_fraction", 1e-8);
      options.set_double("primal_regularization_lb", 1e-20);
      options.set_double("primal_regularization_decrease_factor", 3.);
      options.set_double("primal_regularization_fast_increase_factor", 100.);
      options.set_double("primal_regularization_slow_increase_factor", 8.);
      options.set_integer("threshold_unsuccessful_attempts", 8);
      options.set_bool("regularize_all_variables", true);

      /** iterative refinement options (on the uncondensed Newton system) **/
      // IPOPT: 1 (https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_min_refinement_steps). 0: no refinement step when the
      // solution is already accurate
      options.set_integer("iterative_refinement_min_steps", 0);
      // https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_max_refinement_steps (unlike IPOPT, this is the exact maximum
      // number of steps; 0 disables iterative refinement altogether)
      options.set_integer("iterative_refinement_max_steps", 10);
      options.set_double("iterative_refinement_residual_ratio_max", 1e-10);
      options.set_double("iterative_refinement_residual_ratio_singular", 1e-5);
      options.set_double("iterative_refinement_improvement_factor", 0.999999999);

      /** trust region options **/
      // initial trust region radius
      options.set_double("TR_radius", 10.);
      // TR radius increase factor
      options.set_double("TR_increase_factor", 2);
      // TR radius decrease factor
      options.set_double("TR_decrease_factor", 2);
      // TR aggressive radius decrease factor
      options.set_double("TR_aggressive_decrease_factor", 4);
      // tolerance in TR constraint activity
      options.set_double("TR_activity_tolerance", 1e-6);
      // minimum TR radius
      options.set_double("TR_min_radius", 1e-12);
      // threshold below which the TR radius is reset
      options.set_double("TR_radius_reset_threshold", 1e-4);

      /** feasibility restoration options **/
      // test linearized feasibility when switching back to the optimality phase
      options.set_bool("switch_to_optimality_requires_linearized_feasibility", true);
      options.set_double("l1_constraint_violation_coefficient", 1.);
      options.set_bool("use_proximal_term", false);

      /** barrier subproblem options **/
      options.set_string("barrier_function", "log");
      options.set_double("barrier_initial_parameter", 0.1);
      options.set_double("barrier_default_multiplier", 1);
      // Ipopt parameters
      options.set_double("barrier_tau_min", 0.99);
      options.set_double("barrier_k_sigma", 1e10);
      options.set_double("barrier_k_mu", 0.2);
      options.set_double("barrier_theta_mu", 1.5);
      options.set_double("barrier_k_epsilon", 10);
      options.set_double("barrier_update_fraction", 10);
      options.set_double("barrier_regularization_exponent", 0.25);
      options.set_double("barrier_small_direction_factor", 10.);
      options.set_double("barrier_push_variable_to_interior_k1", 1e-2);
      options.set_double("barrier_push_variable_to_interior_k2", 1e-2);
      options.set_double("barrier_damping_factor", 1e-5);
      options.set_double("barrier_small_infeasibility_factor", 1e-4);
      options.set_double("least_square_multiplier_max_norm", 1e3);
      options.set_double("bound_multiplier_max_norm", 1e3);
      options.set_double("constraint_violation_tolerance", 1e-4);
      options.set_double("bound_relaxation_factor", 1e-8);

      /** BQPD options **/
      options.set_string("BQPD_kmax_heuristic", "filtersqp");

      /** MA57 options **/
      options.set_bool("MA57_use_scaling", false);
   }

   // determine default subproblem solvers, based on the available external dependencies
   void DefaultOptions::determine_subproblem_solvers(Options& options) {
      // QP solver
      if (0 < QPSolverFactory::available_solvers.size()) {
         options.set_string("QP_solver", *QPSolverFactory::available_solvers.begin());
      }
      // LP solver
      if (0 < LPSolverFactory::available_solvers.size()) {
         options.set_string("LP_solver", *LPSolverFactory::available_solvers.begin());
      }
      // linear solver
      // HSL runtime loading: this runs before user options are parsed, so it can't see "libhsl_path" below.
      // Set UNO_HSL_LIBRARY to auto-select MA27/MA57, or choose linear_solver explicitly.
      const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      if (!linear_solvers.empty()) {
         options.set_string("linear_solver", linear_solvers[0]);
      }
      // path to the HSL shared library to dlopen for MA27/MA57 at runtime
      // (empty = platform default libhsl.{so,dylib,dll})
      // options.set_string("libhsl_path", "");
   }
} // namespace