// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
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
      {"function_scaling_lower_bound", OptionType::DOUBLE},
      {"print_minor_iterations", OptionType::BOOL},
      {"print_extended_statistics", OptionType::BOOL},
      {"armijo_decrease_fraction", OptionType::DOUBLE},
      {"armijo_tolerance", OptionType::DOUBLE},
      {"switching_delta", OptionType::DOUBLE},
      {"switching_merit_exponent", OptionType::DOUBLE},
      {"switching_infeasibility_exponent", OptionType::DOUBLE},
      {"sufficient_infeasibility_decrease_ratio", OptionType::DOUBLE},
      {"filter_beta", OptionType::DOUBLE},
      {"filter_gamma", OptionType::DOUBLE},
      {"filter_ubd", OptionType::DOUBLE},
      {"filter_fact", OptionType::DOUBLE},
      {"filter_capacity", OptionType::INTEGER},
      {"filter_sufficient_infeasibility_decrease_factor", OptionType::DOUBLE},
      {"filter_reset_iteration_threshold", OptionType::INTEGER},
      {"max_number_filter_resets", OptionType::INTEGER},
      {"filter_margin_uses_entry_infeasibility", OptionType::BOOL},
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
      {"regularize_all_variables", OptionType::BOOL},
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
   void Options::set_integer(const std::string& option_name, uno_int option_value) {
      check_option_type(option_name, OptionType::INTEGER);
      this->values[option_name] = option_value;
   }

   void Options::set_double(const std::string& option_name, double option_value) {
      check_option_type(option_name, OptionType::DOUBLE);
      this->values[option_name] = option_value;
   }

   void Options::set_bool(const std::string& option_name, bool option_value) {
      check_option_type(option_name, OptionType::BOOL);
      this->values[option_name] = option_value;
   }

   void Options::set_string(const std::string& option_name, const std::string& option_value) {
      check_option_type(option_name, OptionType::STRING);
      this->values[option_name] = option_value;
   }

   namespace {
      uno_int parse_integer(const std::string& option_name, const std::string& option_value) {
         std::istringstream stream(option_value);
         stream.imbue(std::locale::classic());
         long long value{};
         if (!(stream >> value) || !(stream >> std::ws).eof()) {
            throw std::invalid_argument("The value " + option_value + " of option " + option_name + " is not an integer");
         }
         if (value < std::numeric_limits<uno_int>::min() || std::numeric_limits<uno_int>::max() < value) {
            throw std::out_of_range("The value " + option_value + " of option " + option_name + " is out of range");
         }
         return static_cast<uno_int>(value);
      }

      double parse_double(const std::string& option_name, const std::string& option_value) {
         std::istringstream stream(option_value);
         stream.imbue(std::locale::classic());
         double value{};
         if (!(stream >> value) || !(stream >> std::ws).eof()) {
            throw std::invalid_argument("The value " + option_value + " of option " + option_name + " is not a number");
         }
         return value;
      }

      bool parse_bool(const std::string& option_name, const std::string& option_value) {
         if (option_value == "yes" || option_value == "true") {
            return true;
         }
         if (option_value == "no" || option_value == "false") {
            return false;
         }
         throw std::invalid_argument("The value " + option_value + " of option " + option_name +
            " is not a boolean (expected yes/no/true/false)");
      }
   }

   // setter for option with unknown type
   void Options::set(const std::string& option_name, const std::string& option_value) {
      const auto type = option_types.find(option_name);
      if (type == option_types.end()) {
         throw std::out_of_range("The option with name " + option_name + " does not exist");
      }
      switch (type->second) {
         case OptionType::INTEGER:
            this->set_integer(option_name, parse_integer(option_name, option_value));
            break;
         case OptionType::DOUBLE:
            this->set_double(option_name, parse_double(option_name, option_value));
            break;
         case OptionType::BOOL:
            this->set_bool(option_name, parse_bool(option_name, option_value));
            break;
         case OptionType::STRING:
            this->set_string(option_name, option_value);
            break;
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
      const uno_int int_value = this->get<uno_int>(option_name);
      if (int_value < 0) {
         throw std::runtime_error("The unsigned int option with name " + option_name + " is negative");
      }
      return static_cast<size_t>(int_value);
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

   std::optional<uno_int> Options::get_int_optional(const std::string& option_name) const {
      return this->get_optional<uno_int>(option_name);
   }

   std::optional<size_t> Options::get_unsigned_int_optional(const std::string& option_name) const {
      const std::optional<uno_int> int_value = this->get_int_optional(option_name);
      if (!int_value.has_value()) {
         return std::nullopt;
      }
      if (*int_value < 0) {
         throw std::runtime_error("The unsigned int option with name " + option_name + " is negative");
      }
      return static_cast<size_t>(*int_value);
   }

   std::optional<double> Options::get_double_optional(const std::string& option_name) const {
      return this->get_optional<double>(option_name);
   }

   std::optional<bool> Options::get_bool_optional(const std::string& option_name) const {
      return this->get_optional<bool>(option_name);
   }

   std::optional<std::string> Options::get_string_optional(const std::string& option_name) const {
      return this->get_optional<std::string>(option_name);
   }

   OptionType Options::get_option_type(const std::string& option_name) const {
      try {
         return option_types.at(option_name);
      }
      catch (const std::out_of_range&) {
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
      std::ifstream file(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      std::string line;
      size_t line_number = 0;
      while (std::getline(file, line)) {
         ++line_number;
         // Remove comments
         const size_t comment_position = line.find('#');
         if (comment_position != std::string::npos) {
            line.erase(comment_position);
         }
         trim(line);
         if (line.empty()) {
            continue;
         }
         std::istringstream stream(line);
         std::string option_name, option_value, extra_token;
         if (!(stream >> option_name >> option_value) || (stream >> extra_token)) {
            throw std::invalid_argument("Malformed line " + std::to_string(line_number) + " in option file " + file_name + ": " + line);
         }
         options.set(option_name, option_value);
      }
      file.close();
   }

   std::string Options::to_string(const OptionValue& option_value) {
      if (const uno_int* value = std::get_if<uno_int>(&option_value)) {
         return std::to_string(*value);
      }
      if (const double* value = std::get_if<double>(&option_value)) {
         std::ostringstream stream;
         stream << std::setprecision(12) << *value;
         return stream.str();
      }
      if (const bool* value = std::get_if<bool>(&option_value)) {
         return *value ? "true" : "false";
      }
      if (const std::string* value = std::get_if<std::string>(&option_value)) {
         return *value;
      }
      return "<invalid>"; // valueless_by_exception
   }

   void Options::check_option_type(const std::string& option_name, OptionType expected_type) {
      const auto type = option_types.find(option_name);
      if (type == option_types.end()) {
         throw std::out_of_range("The option with name " + option_name + " does not exist");
      }
      if (type->second != expected_type) {
         throw std::invalid_argument("The option with name " + option_name + " has a different type");
      }
   }

   void Options::dump_default_options() {
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

   void Options::print(const std::string& header) const {
      if (!this->values.empty()) {
         // sort the names for deterministic output
         std::vector<std::string> option_names;
         option_names.reserve(this->values.size());
         for (const auto& [option_name, option_value]: this->values) {
            option_names.push_back(option_name);
         }
         std::sort(option_names.begin(), option_names.end());

         std::string option_list{};
         for (const std::string& option_name: option_names) {
            option_list.append(option_name).append(" = ").append(to_string(this->values.at(option_name))).append("\n");
         }
         DISCRETE << header << ":\n" << option_list << '\n';
      }
   }

   OptionOverride::OptionOverride(std::string option_name, std::optional<std::string> old_value, std::string new_value,
      std::string reason): option_name(std::move(option_name)), old_value(std::move(old_value)), new_value(std::move(new_value)),
      reason(std::move(reason)) {
   }
} // namespace