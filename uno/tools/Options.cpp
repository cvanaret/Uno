// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <fstream>
#include <sstream>
#include "solvers/QPSolverFactory.hpp"
#include "solvers/LPSolverFactory.hpp"
#include "solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "Options.hpp"
#include "Logger.hpp"

namespace uno {
   Options::Options(bool are_default_options): are_default_options(are_default_options) { }

   size_t Options::size() const {
      return this->options.size();
   }

   // setter
   std::string& Options::operator[](const std::string& option_name) {
      this->is_default[option_name] = this->are_default_options;
      return this->options[option_name];
   }

   // getter
   const std::string& Options::at(const std::string& option_name) const {
      try {
         const std::string& option_value = this->options.at(option_name);
         this->used[option_name] = true;
         return option_value;
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option " + option_name + " was not found");
      }
   }

   const std::string& Options::get_string(const std::string& option_name) const {
      return this->at(option_name);
   }

   double Options::get_double(const std::string& option_name) const {
      const std::string& entry = this->at(option_name);
      return std::stod(entry);
   }

   int Options::get_int(const std::string& option_name) const {
      const std::string& entry = this->at(option_name);
      return std::stoi(entry);
   }

   size_t Options::get_unsigned_int(const std::string& option_name) const {
      const std::string& entry = this->at(option_name);
      return std::stoul(entry);
   }

   bool Options::get_bool(const std::string& option_name) const {
      const std::string& entry = this->at(option_name);
      return entry == "yes";
   }

   Options Options::get_default() {
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
      options["print_solution"] = "yes";
      // threshold on objective to declare unbounded NLP
      options["unbounded_objective_threshold"] = "-1e20";
      // enforce linear constraints at the initial point (yes|no)
      options["enforce_linear_constraints"] = "no";

      /** statistics table **/
      options["statistics_print_header_every_iterations"] = "15";
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
      options["BQPD_print_subproblem"] = "no";
      options["BQPD_kmax"] = "500";

      /** AMPL options **/
      options["AMPL_write_solution_to_file"] = "yes";

      /** solvers: check the available solvers **/
      // QP solver
      const auto QP_solvers = QPSolverFactory::available_solvers();
      if (not QP_solvers.empty()) {
         options["QP_solver"] = QP_solvers[0];
         options.is_default["QP_solver"] = false;
      }
      // LP solver
      const auto LP_solvers = LPSolverFactory::available_solvers();
      if (not LP_solvers.empty()) {
         options["LP_solver"] = LP_solvers[0];
         options.is_default["LP_solver"] = false;
      }
      // linear solver
      const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      if (not linear_solvers.empty()) {
         options["linear_solver"] = linear_solvers[0];
         options.is_default["linear_solver"] = false;
      }

      /** ingredients **/
      // default constraint relaxation strategy (feasibility_restoration|l1_relaxation)
      options["constraint_relaxation_strategy"] = "feasibility_restoration";
      // default subproblem (QP|LP|primal_dual_interior_point)
      options["subproblem"] = "QP";
      // default globalization strategy (l1_merit|fletcher_filter_method|waechter_filter_method)
      options["globalization_strategy"] = "fletcher_filter_method";
      // default globalization mechanism (TR|LS)
      options["globalization_mechanism"] = "TR";

      return options;
   }

   // argv[i] for i = 3..argc-1 are overwriting options
   Options Options::get_command_line_options(int argc, char* argv[]) {
      static const std::string delimiter = "=";
      Options overwriting_options(false);

      // build the (name, value) map
      for (size_t i = 3; i < static_cast<size_t>(argc); i++) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter + ".");
         }
         const std::string option_name = argument.substr(0, position);
         const std::string option_value = argument.substr(position + 1);
         if (option_name == "preset") {
            Options::set_preset(overwriting_options, option_value);
         }
         else if (option_name == "option_file") {
            Options::overwrite_with_option_file(overwriting_options, option_value);
         }
         else {
            overwriting_options[option_name] = option_value;
         }
      }
      return overwriting_options;
   }

   void Options::overwrite_with_option_file(Options& options, const std::string& file_name) {
      std::ifstream file;
      file.open(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      else {
         // register the default options
         std::string option_name, option_value;
         std::string line;
         while (std::getline(file, line)) {
            if (not line.empty() && line.find('#') != 0) {
               std::istringstream iss;
               iss.str(line);
               iss >> option_name >> option_value;
               options[option_name] = option_value;
            }
         }
         file.close();
      }
   }

   void Options::set_preset(Options& options, const std::string& preset_name) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["subproblem"] = "primal_dual_interior_point";
         options["globalization_mechanism"] = "LS";
         options["globalization_strategy"] = "waechter_filter_method";
         options["filter_type"] = "standard";
         options["filter_beta"] = "0.99999";
         options["filter_gamma"] = "1e-8";
         options["switching_delta"] = "1";
         options["filter_ubd"] = "1e4";
         options["filter_fact"] = "1e4";
         options["filter_switching_infeasibility_exponent"] = "1.1";
         options["armijo_decrease_fraction"] = "1e-8";
         options["LS_backtracking_ratio"] = "0.5";
         options["LS_min_step_length"] = "5e-7";
         options["barrier_tau_min"] = "0.99";
         options["barrier_damping_factor"] = "1e-5";
         options["l1_constraint_violation_coefficient"] = "1000.";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "INF";
         options["scale_functions"] = "yes";
         options["sparse_format"] = "COO";
         options["tolerance"] = "1e-8";
         options["loose_tolerance"] = "1e-6";
         options["loose_tolerance_consecutive_iteration_threshold"] = "15";
         options["switch_to_optimality_requires_linearized_feasibility"] = "no";
         options["LS_scale_duals_with_step_length"] = "yes";
         options["protect_actual_reduction_against_roundoff"] = "yes";
      }
      else if (preset_name == "filtersqp") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["subproblem"] = "QP";
         options["globalization_mechanism"] = "TR";
         options["globalization_strategy"] = "fletcher_filter_method";
         options["filter_type"] = "standard";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L2";
         options["sparse_format"] = "CSC";
         options["TR_radius"] = "10";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["enforce_linear_constraints"] = "yes";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["TR_min_radius"] = "1e-8";
         options["switch_to_optimality_requires_linearized_feasibility"] = "yes";
         options["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "byrd") {
         options["constraint_relaxation_strategy"] = "l1_relaxation";
         options["subproblem"] = "QP";
         options["globalization_mechanism"] = "LS";
         options["globalization_strategy"] = "l1_merit";
         options["l1_relaxation_initial_parameter"] = "1";
         options["LS_backtracking_ratio"] = "0.5";
         options["armijo_decrease_fraction"] = "1e-8";
         options["l1_relaxation_epsilon1"] = "0.1";
         options["l1_relaxation_epsilon2"] = "0.1";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L1";
         options["sparse_format"] = "CSC";
         options["LS_scale_duals_with_step_length"] = "no";
         options["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "funnelsqp") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["subproblem"] = "QP";
         options["globalization_mechanism"] = "TR";
         options["globalization_strategy"] = "funnel_method";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L2";
         options["sparse_format"] = "CSC";
         options["TR_radius"] = "10";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["enforce_linear_constraints"] = "yes";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["TR_min_radius"] = "1e-8";
         options["switch_to_optimality_requires_acceptance"] = "no";
         options["switch_to_optimality_requires_linearized_feasibility"] = "yes";

         options["funnel_beta"] = "0.9999";
         options["funnel_gamma"] = "0.001";
         options["switching_delta"] = "0.999";
         options["funnel_kappa"] = "0.5";
         options["funnel_ubd"] = "1.0";
         options["funnel_fact"] = "1.5";
         options["funnel_switching_infeasibility_exponent"] = "2";
         options["funnel_update_strategy"] = "2";
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known.");
      }
   }

   void Options::overwrite_with(const Options& overwriting_options) {
      for (const auto& [option_name, option_value]: overwriting_options) {
         (*this)[option_name] = option_value;
         this->is_default[option_name] = overwriting_options.is_default[option_name];
      }
   }

   void Options::print_used() const {
      size_t number_used_options = 0;
      std::string option_list{};
      for (const auto& [option_name, option_value]: this->options) {
         if (not this->is_default[option_name] && this->used[option_name]) {
            number_used_options++;
            option_list.append("- ").append(option_name).append(" = ").append(option_value).append("\n");
         }
      }
      // print the overwritten options
      if (number_used_options > 0) {
         DISCRETE << "Used overwritten options:\n" << option_list << '\n';
      }
   }

   std::map<std::string, std::string>::const_iterator Options::begin() const {
      return this->options.begin();
   }

   std::map<std::string, std::string>::const_iterator Options::end() const {
      return this->options.end();
   }
} // namespace
