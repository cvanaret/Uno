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

      /** termination **/
      // primal tolerance (constraint violation)
      options.set_double("primal_tolerance", 1e-8);
      // dual tolerance (stationarity and complementarity)
      options.set_double("dual_tolerance", 1e-8);
      // loose tolerance used if dual tolerance cannot be reached
      options.set_double("loose_dual_tolerance", 1e-6);
      // number of iterations during which the loose tolerance is monitored
      options.set_integer("loose_tolerance_consecutive_iteration_threshold", 15);
      // maximum outer iterations
      options.set_integer("max_iterations", 2000);
      // CPU time limit (in seconds)
      options.set_double("time_limit", INF<double>);
      // print optimal solution (yes|no)
      options.set_bool("print_solution", false);
      // threshold on objective to declare unbounded NLP
      options.set_double("unbounded_objective_threshold", -1e20);
      // enforce linear constraints at the initial point (yes|no)
      options.set_bool("enforce_linear_constraints", false);

      /** statistics table **/
      options.set_integer("statistics_major_column_order", 1);
      options.set_integer("statistics_minor_column_order", 2);
      options.set_integer("statistics_penalty_parameter_column_order", 5);
      options.set_integer("statistics_barrier_parameter_column_order", 8);
      options.set_integer("statistics_SOC_column_order", 9);
      options.set_integer("statistics_TR_radius_column_order", 10);
      options.set_integer("statistics_LS_step_length_column_order", 10);
      options.set_integer("statistics_restoration_phase_column_order", 20);
      options.set_integer("statistics_regularization_column_order", 21);
      options.set_integer("statistics_funnel_width_column_order", 25);
      options.set_integer("statistics_step_norm_column_order", 31);
      options.set_integer("statistics_objective_column_order", 100);
      options.set_integer("statistics_primal_feasibility_column_order", 101);
      options.set_integer("statistics_dual_feasibility_column_order", 102);
      options.set_integer("statistics_stationarity_column_order", 104);
      options.set_integer("statistics_complementarity_column_order", 105);
      options.set_integer("statistics_status_column_order", 200);

      /** main options **/
      // logging level (SILENT|DISCRETE|WARNING|INFO|DEBUG|DEBUG2|DEBUG3)
      options.set_string("logger", "INFO");
      // Hessian model (exact|zero)
      options.set_string("hessian_model", "exact");
      options.set_string("regularization_strategy", "primal");
      // scale the functions (yes|no)
      options.set_bool("scale_functions", false);
      options.set_double("function_scaling_threshold", 100.);
      // factor scaling
      options.set_double("function_scaling_factor", 100.);
      // scale the errors with respect to the current point (yes|no)
      options.set_bool("scale_residuals", true);
      // norm of the progress measures (L1|L2|INF)
      options.set_string("progress_norm", "L1");
      // norm of the primal-dual residuals (L1|L2|INF)
      options.set_string("residual_norm", "INF");
      options.set_double("residual_scaling_threshold", 100.);
      options.set_bool("protect_actual_reduction_against_roundoff", false);
      options.set_bool("print_subproblem", false);

      /** globalization strategy options **/
      options.set_double("armijo_decrease_fraction", 1e-4);
      options.set_double("armijo_tolerance", 1e-9);

      /** switching method options **/
      options.set_double("switching_delta", 0.999);
      options.set_double("switching_infeasibility_exponent", 2);

      /** filter method options **/
      // filter type (standard|nonmonotone)
      options.set_string("filter_type", "standard");
      options.set_double("filter_beta", 0.999);
      options.set_double("filter_gamma", 0.001);
      options.set_double("filter_ubd", 1e2);
      options.set_double("filter_fact", 1.25);
      options.set_integer("filter_capacity", 50);
      // used by Waechter filter method
      options.set_double("filter_sufficient_infeasibility_decrease_factor", 0.9);
      // nonmonotone filter strategy
      options.set_integer("nonmonotone_filter_number_dominated_entries", 3);

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
      // minimum step length
      options.set_double("LS_min_step_length", 1e-12);
      // use the primal-dual and dual step lengths to scale the dual directions when assembling the trial iterate
      options.set_bool("LS_scale_duals_with_step_length", true);

      /** regularization options **/
      // regularization failure threshold
      options.set_double("regularization_failure_threshold", 1e40);
      // Hessian regularization: initial value
      options.set_double("regularization_initial_value", 1e-4);
      options.set_double("regularization_increase_factor", 2);
      // regularization of augmented system
      options.set_double("primal_regularization_initial_factor", 1e-4);
      options.set_double("dual_regularization_fraction", 1e-8);
      options.set_double("primal_regularization_lb", 1e-20);
      options.set_double("primal_regularization_decrease_factor", 3.);
      options.set_double("primal_regularization_fast_increase_factor", 100.);
      options.set_double("primal_regularization_slow_increase_factor", 8.);
      options.set_integer("threshold_unsuccessful_attempts", 8);

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
      options.set_double("TR_min_radius", 1e-7);
      // threshold below which the TR radius is reset
      options.set_double("TR_radius_reset_threshold", 1e-4);

      /** feasibility restoration options **/
      // test linearized feasibility when switching back to the optimality phase
      options.set_bool("switch_to_optimality_requires_linearized_feasibility", true);
      options.set_double("l1_constraint_violation_coefficient", 1);

      /** barrier subproblem options **/
      options.set_double("barrier_initial_parameter", 0.1);
      options.set_double("barrier_default_multiplier", 1);
      // Ipopt parameters
      options.set_double("barrier_tau_min", 0.99);
      options.set_double("barrier_k_sigma", 1e10);
      options.set_double("barrier_smax", 100);
      options.set_double("barrier_k_mu", 0.2);
      options.set_double("barrier_theta_mu", 1.5);
      options.set_double("barrier_k_epsilon", 10);
      options.set_double("barrier_update_fraction", 10);
      options.set_double("barrier_regularization_exponent", 0.25);
      options.set_double("barrier_small_direction_factor", 10.);
      options.set_double("barrier_push_variable_to_interior_k1", 1e-2);
      options.set_double("barrier_push_variable_to_interior_k2", 1e-2);
      options.set_double("barrier_damping_factor", 1e-5);
      options.set_double("least_square_multiplier_max_norm", 1e3);

      /** BQPD options **/
      options.set_integer("BQPD_kmax", 500);
   }

   // determine default subproblem solvers, based on the available external dependencies
   void DefaultOptions::determine_subproblem_solvers(Options& options) {
      // QP solver
      if (0 < QPSolverFactory::available_solvers.size()) {
         options.set_string("QP_solver", *QPSolverFactory::available_solvers.begin(), true);
      }
      // LP solver
      if (0 < LPSolverFactory::available_solvers.size()) {
         options.set_string("LP_solver", *LPSolverFactory::available_solvers.begin(), true);
      }
      // linear solver
      const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      if (!linear_solvers.empty()) {
         options.set_string("linear_solver", linear_solvers[0], true);
      }
   }
} // namespace