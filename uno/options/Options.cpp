// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <fstream>
#include <iostream>
#include <sstream>
#include "Options.hpp"
#include "DefaultOptions.hpp"
#include "Presets.hpp"
#include "tools/Logger.hpp"

namespace uno {
   const std::unordered_map<std::string, OptionType> Options::option_types = {
      {"primal_tolerance", OptionType::DOUBLE},
      {"dual_tolerance", OptionType::DOUBLE},
      {"loose_primal_tolerance", OptionType::DOUBLE},
      {"loose_dual_tolerance", OptionType::DOUBLE},
      {"loose_tolerance_consecutive_iteration_threshold", OptionType::INTEGER},
      {"max_iterations", OptionType::INTEGER},
      {"time_limit", OptionType::DOUBLE},
      {"print_solution", OptionType::BOOL},
      {"unbounded_objective_threshold", OptionType::DOUBLE},
      {"enforce_linear_constraints", OptionType::BOOL},
      {"logger", OptionType::STRING},
      {"constraint_relaxation_strategy", OptionType::STRING},
      {"inequality_handling_method", OptionType::STRING},
      {"globalization_mechanism",OptionType::STRING},
      {"globalization_strategy", OptionType::STRING},
      {"hessian_model", OptionType::STRING},
      {"inertia_correction_strategy", OptionType::STRING},
      {"scale_functions", OptionType::BOOL},
      {"function_scaling_threshold", OptionType::DOUBLE},
      {"function_scaling_factor", OptionType::DOUBLE},
      {"scale_residuals", OptionType::BOOL},
      {"progress_norm", OptionType::STRING},
      {"residual_norm", OptionType::STRING},
      {"residual_scaling_threshold", OptionType::DOUBLE},
      {"protect_actual_reduction_against_roundoff", OptionType::BOOL},
      {"print_subproblem", OptionType::BOOL},
      {"armijo_decrease_fraction", OptionType::DOUBLE},
      {"armijo_tolerance", OptionType::DOUBLE},
      {"switching_delta", OptionType::DOUBLE},
      {"switching_infeasibility_exponent", OptionType::DOUBLE},
      {"filter_type", OptionType::STRING},
      {"filter_beta", OptionType::DOUBLE},
      {"filter_gamma", OptionType::DOUBLE},
      {"filter_ubd", OptionType::DOUBLE},
      {"filter_fact", OptionType::DOUBLE},
      {"filter_capacity", OptionType::INTEGER},
      {"filter_sufficient_infeasibility_decrease_factor", OptionType::DOUBLE},
      {"nonmonotone_filter_number_dominated_entries", OptionType::INTEGER},
      {"funnel_kappa", OptionType::DOUBLE},
      {"funnel_beta", OptionType::DOUBLE},
      {"funnel_gamma", OptionType::DOUBLE},
      {"funnel_ubd", OptionType::DOUBLE},
      {"funnel_fact", OptionType::DOUBLE},
      {"funnel_update_strategy", OptionType::INTEGER},
      {"funnel_require_acceptance_wrt_current_iterate", OptionType::BOOL},
      {"LS_backtracking_ratio", OptionType::DOUBLE},
      {"LS_min_step_length", OptionType::DOUBLE},
      {"LS_scale_duals_with_step_length", OptionType::BOOL},
      {"regularization_failure_threshold", OptionType::DOUBLE},
      {"regularization_initial_value", OptionType::DOUBLE},
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
      {"barrier_initial_parameter", OptionType::DOUBLE},
      {"barrier_default_multiplier", OptionType::DOUBLE},
      {"barrier_tau_min", OptionType::DOUBLE},
      {"barrier_k_sigma", OptionType::DOUBLE},
      {"barrier_smax", OptionType::DOUBLE},
      {"barrier_k_mu", OptionType::DOUBLE},
      {"barrier_theta_mu", OptionType::DOUBLE},
      {"barrier_k_epsilon", OptionType::DOUBLE},
      {"barrier_update_fraction", OptionType::DOUBLE},
      {"barrier_regularization_exponent", OptionType::DOUBLE},
      {"barrier_small_direction_factor", OptionType::DOUBLE},
      {"barrier_push_variable_to_interior_k1", OptionType::DOUBLE},
      {"barrier_push_variable_to_interior_k2", OptionType::DOUBLE},
      {"barrier_damping_factor", OptionType::DOUBLE},
      {"least_square_multiplier_max_norm", OptionType::DOUBLE},
      {"BQPD_kmax_heuristic", OptionType::STRING},
      {"QP_solver", OptionType::STRING},
      {"LP_solver", OptionType::STRING},
      {"linear_solver", OptionType::STRING},
      {"preset", OptionType::STRING},
   };

   // setters
   void Options::set_integer(const std::string& option_name, uno_int option_value, bool flag_as_overwritten) {
      this->integer_options[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::dump_default_options() {
       Options defaults;
       DefaultOptions::load(defaults);
       // Header (tab-separated for easy parsing)
       std::cout << "name\ttype\tdefault\n";
       for (const auto &[option_name, option_type] : defaults.option_types) {
           std::cout << option_name << '\t';
           try {
               switch (option_type) {
                   case OptionType::INTEGER: 
                       std::cout << "INTEGER\t";
                       if (defaults.integer_options.find(option_name) != defaults.integer_options.end())
                           std::cout << defaults.integer_options.at(option_name);
                       else
                           std::cout << "<unset>";
                       break;
                   case OptionType::DOUBLE:
                       std::cout << "DOUBLE\t";
                       if (defaults.double_options.find(option_name) != defaults.double_options.end())
                           std::cout << defaults.double_options.at(option_name);
                       else
                           std::cout << "<unset>";
                       break;
                   case OptionType::BOOL:
                       std::cout << "BOOL\t";
                       if (defaults.bool_options.find(option_name) != defaults.bool_options.end())
                           std::cout << (defaults.bool_options.at(option_name) ? "true"
                                                                               : "false");
                       else
                           std::cout << "<unset>";
                       break;
                   case OptionType::STRING:
                       std::cout << "STRING\t";
                       if (defaults.string_options.find(option_name) != defaults.string_options.end())
                           std::cout << defaults.string_options.at(option_name);
                       else
                           std::cout << "<unset>";
                       break;
                   default:
                       std::cout << "UNKNOWN\t<unset>";
                       break;
               }
           } catch (const std::out_of_range &) {
               std::cout << "<error>";
           }
           std::cout << '\n';
       }
   }
   
   void Options::set_double(const std::string& option_name, double option_value, bool flag_as_overwritten) {
      this->double_options[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::set_bool(const std::string& option_name, bool option_value, bool flag_as_overwritten) {
      this->bool_options[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::set_string(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten) {
      this->string_options[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   // setter for option with unknown type
   void Options::set(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten) {
      try {
         const OptionType type = Options::option_types.at(option_name);
         if (type == OptionType::INTEGER) {
            this->set_integer(option_name, std::stoi(option_value), flag_as_overwritten);
         }
         else if (type == OptionType::DOUBLE) {
            this->set_double(option_name, std::stod(option_value), flag_as_overwritten);
         }
         else if (type == OptionType::BOOL) {
            this->set_bool(option_name, option_value == "yes", flag_as_overwritten);
         }
         else if (type == OptionType::STRING) {
            this->set_string(option_name, option_value, flag_as_overwritten);
         }
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The type of the option with name " + option_name + " could not be found");
      }
   }

   // getters
   uno_int Options::get_int(const std::string& option_name) const {
      this->used[option_name] = true;
      try {
         return this->integer_options.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   size_t Options::get_unsigned_int(const std::string& option_name) const {
      this->used[option_name] = true;
      try {
         return static_cast<size_t>(this->integer_options.at(option_name));
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   double Options::get_double(const std::string& option_name) const {
      this->used[option_name] = true;
      try {
         return this->double_options.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   bool Options::get_bool(const std::string& option_name) const {
      this->used[option_name] = true;
      try {
         return this->bool_options.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   const std::string& Options::get_string(const std::string& option_name) const {
      this->used[option_name] = true;
      try {
         return this->string_options.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   std::optional<std::string> Options::get_string_optional(const std::string& option_name) const {
      try {
         const std::string& option_value = this->string_options.at(option_name);
         this->used[option_name] = true;
         return option_value;
      }
      catch(const std::out_of_range&) {
         return std::nullopt;
      }
   }

   OptionType Options::get_option_type(const std::string& option_name) const {
      try {
         return this->option_types.at(option_name);
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The type of the option with name " + option_name + " could not be found");
      }
   }

   // argv[i] for i = offset..argc-1 are overwriting options
   Options Options::get_command_line_options(int argc, char* argv[], size_t offset) {
      static const std::string delimiter = "=";
      Options command_line_options;

      // build the (name, value) map
      for (size_t i = offset; i < static_cast<size_t>(argc); ++i) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter);
         }
         const std::string option_name = argument.substr(0, position);
         const std::string option_value = argument.substr(position + 1);
         // set option with unknown type
         command_line_options.set(option_name, option_value);
      }
      return command_line_options;
   }

   Options Options::load_option_file(const std::string& file_name) {
      Options options;
      std::ifstream file;
      file.open(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      else {
         std::string option_name, option_value;
         std::string line;
         while (std::getline(file, line)) {
            if (!line.empty() && line.find('#') != 0) {
               std::istringstream iss;
               iss.str(line);
               iss >> option_name >> option_value;
               // set option with unknown type
               options.set(option_name, option_value);
            }
         }
         file.close();
      }
      return options;
   }

   void Options::overwrite_with(const Options& overwriting_options) {
      for (const auto& [option_name, option_value]: overwriting_options.integer_options) {
         this->integer_options[option_name] = option_value;
      }
      for (const auto& [option_name, option_value]: overwriting_options.double_options) {
         this->double_options[option_name] = option_value;
      }
      for (const auto& [option_name, option_value]: overwriting_options.bool_options) {
         this->bool_options[option_name] = option_value;
      }
      for (const auto& [option_name, option_value]: overwriting_options.string_options) {
         this->string_options[option_name] = option_value;
      }
      // if the option already exists and is not the same, flag it as overwritten
      //const auto existing_option_value = this->at_optional(option_name);
      //bool flag_as_overwritten = (existing_option_value.has_value() && *existing_option_value != option_value);
      //this->set(option_name, option_value, flag_as_overwritten);
   }

   void Options::print_non_default() const {
      size_t number_used_options = 0;
      std::string option_list{};
      for (const auto& [option_name, option_value]: this->string_options) {
         if (this->used[option_name] && this->overwritten_options[option_name]) {
            ++number_used_options;
            option_list.append(option_name).append(" = ").append(option_value).append("\n");
         }
      }
      // print the overwritten options
      if (number_used_options > 0) {
         DISCRETE << "\nNon-default options:\n" << option_list << '\n';
      }
   }
} // namespace
