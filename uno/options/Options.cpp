// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <fstream>
#include <sstream>
#include "Options.hpp"
#include "tools/Logger.hpp"

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
         throw std::out_of_range("The option with name " + option_name + " was not found");
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

   // argv[i] for i = 3..argc-1 are overwriting options
   Options Options::get_command_line_options(int argc, char* argv[]) {
      static const std::string delimiter = "=";
      Options overwriting_options(false);

      // build the (name, value) map
      for (size_t i = 3; i < static_cast<size_t>(argc); i++) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter);
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
         throw std::runtime_error("The preset " + preset_name + " is not known");
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
         DISCRETE << "\nUsed overwritten options:\n" << option_list << '\n';
      }
   }

   std::map<std::string, std::string>::const_iterator Options::begin() const {
      return this->options.begin();
   }

   std::map<std::string, std::string>::const_iterator Options::end() const {
      return this->options.end();
   }
} // namespace
