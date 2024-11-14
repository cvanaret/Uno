// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DefaultOptions.hpp"
#include "solvers/QPSolverFactory.hpp"
#include "solvers/LPSolverFactory.hpp"
#include "solvers/SymmetricIndefiniteLinearSolverFactory.hpp"

namespace uno {
   Options DefaultOptions::load() {
      Options options(true);

      /** termination **/
      // (tight) tolerance
      options["tolerance"] = "1e-8";
      // loose tolerance used if tight tolerance cannot be reached
      options["loose_tolerance"] = "1e-6";
      // number of iterations during which the loose tolerance is monitored
      options["loose_tolerance_consecutive_iteration_threshold"] = "15";
      // maximum outer iterations
      options["max_iterations"] = "2000";
      // CPU time limit (in seconds)
      options["time_limit"] = "inf";
      // print optimal solution (yes|no)
      options["print_solution"] = "no";
      // threshold on objective to declare unbounded NLP
      options["unbounded_objective_threshold"] = "-1e20";
      // enforce linear constraints at the initial point (yes|no)
      options["enforce_linear_constraints"] = "no";

      /** statistics table **/
      options["statistics_print_header_frequency"] = "15";
      options["statistics_major_column_order"] = "1";
      options["statistics_minor_column_order"] = "2";
      options["statistics_penalty_parameter_column_order"] = "5";
      options["statistics_barrier_parameter_column_order"] = "8";
      options["statistics_SOC_column_order"] = "9";
      options["statistics_TR_radius_column_order"] = "10";
      options["statistics_LS_step_length_column_order"] = "10";
      options["statistics_restoration_phase_column_order"] = "20";
      options["statistics_regularization_column_order"] = "21";
      options["statistics_funnel_width_column_order"] = "25";
      options["statistics_step_norm_column_order"] = "31";
      options["statistics_objective_column_order"] = "100";
      options["statistics_primal_feasibility_column_order"] = "101";
      options["statistics_dual_feasibility_column_order"] = "102";
      options["statistics_stationarity_column_order"] = "104";
      options["statistics_complementarity_column_order"] = "105";
      options["statistics_status_column_order"] = "200";

      /** main options **/
      // logging level (SILENT|DISCRETE|WARNING|INFO|DEBUG|DEBUG2|DEBUG3)
      options["logger"] = "INFO";
      // Hessian model (exact|zero)
      options["hessian_model"] = "exact";
      // sparse matrix format (COO|CSC)
      options["sparse_format"] = "COO";
      // scale the functions (yes|no)
      options["scale_functions"] = "no";
      options["function_scaling_threshold"] = "100";
      // factor scaling
      options["function_scaling_factor"] = "100";
      // scale the errors with respect to the current point (yes|no)
      options["scale_residuals"] = "yes";
      // norm of the progress measures (L1|L2|INF)
      options["progress_norm"] = "L1";
      // norm of the primal-dual residuals (L1|L2|INF)
      options["residual_norm"] = "INF";
      options["residual_scaling_threshold"] = "100.";
      options["protect_actual_reduction_against_roundoff"] = "no";
      options["print_subproblem"] = "no";

      /** globalization strategy options **/
      options["armijo_decrease_fraction"] = "1e-4";
      options["armijo_tolerance"] = "1e-9";

      /** switching method options **/
      options["switching_delta"] = "0.999";
      options["switching_infeasibility_exponent"] = "2";

      /** filter method options **/
      // filter type (standard|nonmonotone)
      options["filter_type"] = "standard";
      options["filter_beta"] = "0.999";
      options["filter_gamma"] = "0.001";
      options["filter_ubd"] = "1e2";
      options["filter_fact"] = "1.25";
      options["filter_capacity"] = "50";
      // used by Waechter filter method
      options["filter_sufficient_infeasibility_decrease_factor"] = "0.9";
      // nonmonotone filter strategy
      options["nonmonotone_filter_number_dominated_entries"] = "3";

      /** funnel options **/
      options["funnel_kappa"] = "0.5";
      options["funnel_beta"] = "0.9999";
      options["funnel_gamma"] = "0.001";
      options["funnel_ubd"] = "1.0";
      options["funnel_fact"] = "1.5";
      options["funnel_update_strategy"] = "1";
      options["funnel_require_acceptance_wrt_current_iterate"] = "no";

      /** line search options */
      // backtracking ratio
      options["LS_backtracking_ratio"] = "0.5";
      // minimum step length
      options["LS_min_step_length"] = "1e-12";
      // use the primal-dual and dual step lengths to scale the dual directions when assembling the trial iterate
      options["LS_scale_duals_with_step_length"] = "yes";

      /** regularization options **/
      // regularization failure threshold
      options["regularization_failure_threshold"] = "1e40";
      // Hessian regularization: initial value
      options["regularization_initial_value"] = "1e-4";
      options["regularization_increase_factor"] = "2";
      // regularization of augmented system
      options["primal_regularization_initial_factor"] = "1e-4";
      options["dual_regularization_fraction"] = "1e-8";
      options["primal_regularization_lb"] = "1e-20";
      options["primal_regularization_decrease_factor"] = "3.";
      options["primal_regularization_fast_increase_factor"] = "100.";
      options["primal_regularization_slow_increase_factor"] = "8.";
      options["threshold_unsuccessful_attempts"] = "8";

      /** trust region options **/
      // initial trust region radius
      options["TR_radius"] = "10.";
      // TR radius increase factor
      options["TR_increase_factor"] = "2";
      // TR radius decrease factor
      options["TR_decrease_factor"] = "2";
      // TR aggressive radius decrease factor
      options["TR_aggressive_decrease_factor"] = "4";
      // tolerance in TR constraint activity
      options["TR_activity_tolerance"] = "1e-6";
      // minimum TR radius
      options["TR_min_radius"] = "1e-7";
      // threshold below which the TR radius is reset
      options["TR_radius_reset_threshold"] = "1e-4";
      // force QP convexification when in a trust-region setting
      options["convexify_QP"] = "false";

      /** constraint relaxation options **/
      // l1 relaxation options //
      // initial value of the penalty parameter
      options["l1_relaxation_initial_parameter"] = "1.";
      // use a fixed parameter (yes|no)
      options["l1_relaxation_fixed_parameter"] = "no";
      // decrease (multiplicative) factor of penalty parameter
      options["l1_relaxation_decrease_factor"] = "10.";
      // epsilon constants in Byrd's article
      options["l1_relaxation_epsilon1"] = "0.1";
      options["l1_relaxation_epsilon2"] = "0.1";
      options["l1_relaxation_residual_small_threshold"] = "1e-12";
      // coefficient of constraint violation
      options["l1_constraint_violation_coefficient"] = "1";
      // threshold for determining if duals have a zero norm
      options["l1_small_duals_threshold"] = "1e-10";

      /** feasibility restoration options **/
      // test linearized feasibility when switching back to the optimality phase
      options["switch_to_optimality_requires_linearized_feasibility"] = "yes";

      /** barrier subproblem options **/
      options["barrier_initial_parameter"] = "0.1";
      options["barrier_default_multiplier"] = "1";
      // Ipopt parameters
      options["barrier_tau_min"] = "0.99";
      options["barrier_k_sigma"] = "1e10";
      options["barrier_smax"] = "100";
      options["barrier_k_mu"] = "0.2";
      options["barrier_theta_mu"] = "1.5";
      options["barrier_k_epsilon"] = "10";
      options["barrier_update_fraction"] = "10";
      options["barrier_regularization_exponent"] = "0.25";
      options["barrier_small_direction_factor"] = "10.";
      options["barrier_push_variable_to_interior_k1"] = "1e-2";
      options["barrier_push_variable_to_interior_k2"] = "1e-2";
      options["barrier_damping_factor"] = "1e-5";
      options["least_square_multiplier_max_norm"] = "1e3";

      /** BQPD options **/
      options["BQPD_kmax"] = "500";

      /** AMPL options **/
      options["AMPL_write_solution_to_file"] = "yes";

      return options;
   }

   Options DefaultOptions::determine_solvers() {
      Options options(false);

      /** solvers: check the available solvers **/
      // QP solver
      const auto QP_solvers = QPSolverFactory::available_solvers();
      if (not QP_solvers.empty()) {
         options["QP_solver"] = QP_solvers[0];
      }
      // LP solver
      const auto LP_solvers = LPSolverFactory::available_solvers();
      if (not LP_solvers.empty()) {
         options["LP_solver"] = LP_solvers[0];
      }
      // linear solver
      const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      if (not linear_solvers.empty()) {
         options["linear_solver"] = linear_solvers[0];
      }
      return options;
   }
} // namespace
