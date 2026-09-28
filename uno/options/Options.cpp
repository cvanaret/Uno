// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include "Options.hpp"
#include "DefaultOptions.hpp"
#include "tools/Logger.hpp"

namespace uno {
   const std::unordered_map<std::string, OptionType> Options::option_types = {
      {"primal_tolerance", OptionType::DOUBLE},
      {"dual_tolerance", OptionType::DOUBLE},
      {"loose_primal_tolerance", OptionType::DOUBLE},
      {"loose_dual_tolerance", OptionType::DOUBLE},
      {"loose_tolerance_iteration_threshold", OptionType::INTEGER},
      {"max_iterations", OptionType::INTEGER},
      {"time_limit", OptionType::DOUBLE},
      {"print_solution", OptionType::BOOL},
      {"diverging_iterate_threshold", OptionType::DOUBLE},
      {"unbounded_objective_threshold", OptionType::DOUBLE},
      {"logger", OptionType::STRING},
      {"constraint_relaxation_strategy", OptionType::STRING},
      {"inequality_handling_method", OptionType::STRING},
      {"globalization_mechanism",OptionType::STRING},
      {"globalization_strategy", OptionType::STRING},
      {"hessian_model", OptionType::STRING},
      {"inertia_correction_strategy", OptionType::STRING},
      {"progress_norm", OptionType::STRING},
      {"residual_norm", OptionType::STRING},
      {"residual_scaling_threshold", OptionType::DOUBLE},
      {"protect_actual_reduction_against_roundoff", OptionType::BOOL},
      {"protected_actual_reduction_macheps_coefficient", OptionType::DOUBLE},
      {"print_subproblem", OptionType::BOOL},
      {"use_function_scaling", OptionType::BOOL},
      {"function_scaling_threshold", OptionType::DOUBLE},
      {"print_minor_iterations", OptionType::BOOL},
      {"print_extended_statistics", OptionType::BOOL},
      {"armijo_decrease_fraction", OptionType::DOUBLE},
      {"armijo_tolerance", OptionType::DOUBLE},
      {"switching_delta", OptionType::DOUBLE},
      {"switching_objective_exponent", OptionType::DOUBLE},
      {"switching_infeasibility_exponent", OptionType::DOUBLE},
      {"sufficient_infeasibility_decrease_ratio", OptionType::DOUBLE},
      {"filter_type", OptionType::STRING},
      {"filter_beta", OptionType::DOUBLE},
      {"filter_gamma", OptionType::DOUBLE},
      {"filter_ubd", OptionType::DOUBLE},
      {"filter_fact", OptionType::DOUBLE},
      {"filter_capacity", OptionType::INTEGER},
      {"filter_sufficient_infeasibility_decrease_factor", OptionType::DOUBLE},
      {"filter_reset_iteration_threshold", OptionType::INTEGER},
      {"max_number_filter_resets", OptionType::INTEGER},
      {"funnel_kappa", OptionType::DOUBLE},
      {"funnel_beta", OptionType::DOUBLE},
      {"funnel_gamma", OptionType::DOUBLE},
      {"funnel_ubd", OptionType::DOUBLE},
      {"funnel_fact", OptionType::DOUBLE},
      {"funnel_update_strategy", OptionType::INTEGER},
      {"funnel_require_acceptance_wrt_current_iterate", OptionType::BOOL},
      {"LS_backtracking_ratio", OptionType::DOUBLE},
      {"LS_scale_duals_with_step_length", OptionType::BOOL},
      {"gamma_alpha", OptionType::DOUBLE},
      {"quasi_newton_memory_size", OptionType::INTEGER},
      {"quasi_newton_delta_lower_bound", OptionType::DOUBLE},
      {"quasi_newton_delta_upper_bound", OptionType::DOUBLE},
      {"LBFGS_max_skips_before_reset", OptionType::INTEGER},
      {"LSR1_pivot_max_magnitude", OptionType::DOUBLE},
      {"regularization_failure_threshold", OptionType::DOUBLE},
      {"regularization_increase_factor", OptionType::DOUBLE},
      {"primal_regularization_initial_factor", OptionType::DOUBLE},
      {"dual_regularization_fraction", OptionType::DOUBLE},
      {"primal_regularization_lb", OptionType::DOUBLE},
      {"primal_regularization_decrease_factor", OptionType::DOUBLE},
      {"primal_regularization_fast_increase_factor", OptionType::DOUBLE},
      {"primal_regularization_slow_increase_factor", OptionType::DOUBLE},
      {"threshold_unsuccessful_attempts", OptionType::INTEGER},
      {"TR_radius", OptionType::DOUBLE},
      {"TR_increase_factor", OptionType::DOUBLE},
      {"TR_decrease_factor", OptionType::DOUBLE},
      {"TR_aggressive_decrease_factor", OptionType::DOUBLE},
      {"TR_activity_tolerance", OptionType::DOUBLE},
      {"TR_min_radius", OptionType::DOUBLE},
      {"TR_radius_reset_threshold", OptionType::DOUBLE},
      {"switch_to_optimality_requires_linearized_feasibility", OptionType::BOOL},
      {"l1_constraint_violation_coefficient", OptionType::DOUBLE},
      {"use_proximal_term", OptionType::BOOL},
      {"barrier_function", OptionType::STRING},
      {"barrier_initial_parameter", OptionType::DOUBLE},
      {"barrier_default_multiplier", OptionType::DOUBLE},
      {"barrier_tau_min", OptionType::DOUBLE},
      {"barrier_k_sigma", OptionType::DOUBLE},
      {"barrier_k_mu", OptionType::DOUBLE},
      {"barrier_theta_mu", OptionType::DOUBLE},
      {"barrier_k_epsilon", OptionType::DOUBLE},
      {"barrier_update_fraction", OptionType::DOUBLE},
      {"barrier_regularization_exponent", OptionType::DOUBLE},
      {"barrier_small_direction_factor", OptionType::DOUBLE},
      {"barrier_push_variable_to_interior_k1", OptionType::DOUBLE},
      {"barrier_push_variable_to_interior_k2", OptionType::DOUBLE},
      {"barrier_damping_factor", OptionType::DOUBLE},
      {"barrier_small_infeasibility_factor", OptionType::DOUBLE},
      {"least_square_multiplier_max_norm", OptionType::DOUBLE},
      {"bound_multiplier_max_norm", OptionType::DOUBLE},
      {"constraint_violation_tolerance", OptionType::DOUBLE},
      {"bound_relaxation_factor", OptionType::DOUBLE},
      {"BQPD_kmax_heuristic", OptionType::STRING},
      {"QP_solver", OptionType::STRING},
      {"LP_solver", OptionType::STRING},
      {"linear_solver", OptionType::STRING},
      {"preset", OptionType::STRING},
      {"option_file", OptionType::STRING},
      {"write_solution_to_file", OptionType::BOOL},
      {"SOC_max_iterations", OptionType::INTEGER},
      {"SOC_infeasibility_fraction", OptionType::DOUBLE},
      {"libhsl_path", OptionType::STRING},
      {"MA57_use_scaling", OptionType::BOOL},
   };

   // setters
   void Options::set_integer(const std::string& option_name, uno_int option_value, bool flag_as_overwritten) {
      this->values[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::set_double(const std::string& option_name, double option_value, bool flag_as_overwritten) {
      this->values[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::set_bool(const std::string& option_name, bool option_value, bool flag_as_overwritten) {
      this->values[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::set_string(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten) {
      this->values[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   // setter for option with unknown type
   void Options::set(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten) {
      try {
         const OptionType type = option_types.at(option_name);
         if (type == OptionType::INTEGER) {
            this->set_integer(option_name, std::stoi(option_value), flag_as_overwritten);
         }
         else if (type == OptionType::DOUBLE) {
            this->set_double(option_name, std::stod(option_value), flag_as_overwritten);
         }
         else if (type == OptionType::BOOL) {
            this->set_bool(option_name, option_value == "yes" || option_value == "true", flag_as_overwritten);
         }
         else if (type == OptionType::STRING) {
            this->set_string(option_name, option_value, flag_as_overwritten);
         }
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The type of the option with name " + option_name + " could not be found");
      }
   }

   void Options::overwrite(const Options& other) {
      for (const auto& [key, value]: other.values) {
         this->values[key] = value;
      }
   }

   // getters
   uno_int Options::get_int(const std::string& option_name) const {
      return this->get<uno_int>(option_name);
   }

   size_t Options::get_unsigned_int(const std::string& option_name) const {
      return static_cast<size_t>(this->get<uno_int>(option_name));
   }

   double Options::get_double(const std::string& option_name) const {
      return this->get<double>(option_name);
   }

   bool Options::get_bool(const std::string& option_name) const {
      return this->get<bool>(option_name);
   }

   const std::string& Options::get_string(const std::string& option_name) const {
      return this->get<std::string>(option_name);
   }

   std::optional<std::string> Options::get_string_optional(const std::string& option_name) const {
      return this->get_optional<std::string>(option_name);
   }

   OptionType Options::get_option_type(const std::string& option_name) const {
      try {
         return option_types.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The type of the option with name " + option_name + " could not be found");
      }
   }

   // argv[i] for i = offset..argc-1 are overwriting options
   std::vector<std::pair<std::string, std::string>> Options::get_command_line_options(int argc, char* argv[], size_t offset) {
      static const std::string delimiter = "=";
      std::vector<std::pair<std::string, std::string>> command_line_options;

      // build the (name, value) map
      for (size_t i = offset; i < static_cast<size_t>(argc); ++i) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter);
         }
         const std::string option_name = argument.substr(0, position);
         const std::string option_value = argument.substr(position + 1);
         command_line_options.emplace_back(option_name, option_value);
      }
      return command_line_options;
   }

   void trim(std::string& string) {
      string.erase(0, string.find_first_not_of(" \n\r\t"));
      string.erase(string.find_last_not_of(" \n\r\t")+1);
   }

   void Options::load_option_file(Options& options, const std::string& file_name) {
      std::ifstream file;
      file.open(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      else {
         std::string option_name, option_value;
         std::string line;
         while (std::getline(file, line)) {
            // Remove comments
            const size_t comment_position = line.find('#');
            if (comment_position != std::string::npos) {
               line.erase(comment_position);
            }
            trim(line);
            std::istringstream iss;
            iss.str(line);
            if (iss >> option_name >> option_value) {
               // set option (with unknown type)
               options.set(option_name, option_value);
            }
         }
         file.close();
      }
   }

   std::string Options::to_string(const OptionValue& value) {
      return std::visit([](const auto& typed_value) -> std::string {
         using T = std::decay_t<decltype(typed_value)>;
         if constexpr (std::is_same_v<T, bool>) {
            return typed_value ? "true" : "false";
         }
         else if constexpr (std::is_same_v<T, std::string>) {
            return typed_value;
         }
         else { // int, double
            std::ostringstream stream;
            stream << std::setprecision(std::numeric_limits<double>::digits10) << typed_value;
            return stream.str();
         }
      }, value);
   }

   void Options::dump_default_options() {
      std::cout << "preset\tauto\n";
      Options default_options;
      DefaultOptions::load(default_options);

      // sort the names for deterministic output
      std::vector<std::string> option_names;
      option_names.reserve(default_options.values.size());
      for (const auto& [option_name, option_value]: default_options.values) {
         option_names.push_back(option_name);
      }
      std::sort(option_names.begin(), option_names.end());

      for (const std::string& option_name: option_names) {
         std::cout << option_name << '\t' << to_string(default_options.values.at(option_name)) << '\n';
      }
   }

   void Options::print_non_default() const {
      size_t number_used_options = 0;
      std::string option_list{};
      for (const auto& [option_name, option_value]: this->values) {
         if (this->used[option_name] && this->overwritten_options[option_name]) {
            ++number_used_options;
            option_list.append(option_name).append(" = ").append(to_string(option_value)).append("\n");
         }
      }
      // print the overwritten options
      DISCRETE << '\n';
      if (number_used_options > 0) {
         DISCRETE << "Non-default options:\n" << option_list << '\n';
      }
   }
} // namespace