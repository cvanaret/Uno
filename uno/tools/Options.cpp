// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <fstream>
#include <sstream>
#include "Options.hpp"

namespace uno {
   std::string& Options::operator[](const std::string& key) {
      return this->options[key];
   }

   const std::string& Options::at(const std::string& key) const {
      try {
         const std::string& value = this->options.at(key);
         this->is_used[key] = true;
         return value;
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option " + key + " was not found");
      }
   }

   const std::string& Options::get_string(const std::string& key) const {
      return this->at(key);
   }

   double Options::get_double(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stod(entry);
   }

   int Options::get_int(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stoi(entry);
   }

   size_t Options::get_unsigned_int(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stoul(entry);
   }

   bool Options::get_bool(const std::string& key) const {
      const std::string& entry = this->at(key);
      return entry == "yes";
   }

   void Options::print(bool only_used) const {
      std::cout << "Options:\n";
      for (const auto& [key, value]: this->options) {
         if (not only_used || this->is_used[key]) {
            std::cout << "- " << key << " = " << value << '\n';
         }
      }
   }

   Options Options::get_default_options(const std::string& file_name) {
      std::ifstream file;
      file.open(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      else {
         // register the default options
         Options options;
         std::string key, value;
         std::string line;
         while (std::getline(file, line)) {
            if (not line.empty() && line.find('#') != 0) {
               std::istringstream iss;
               iss.str(line);
               iss >> key >> value;
               options[key] = value;
            }
         }
         file.close();
         return options;
      }
   }

   void Options::find_preset(const std::string& preset_name) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         (*this)["constraint_relaxation_strategy"] = "feasibility_restoration";
         (*this)["subproblem"] = "primal_dual_interior_point";
         (*this)["globalization_mechanism"] = "LS";
         (*this)["globalization_strategy"] = "waechter_filter_method";
         (*this)["filter_type"] = "standard";
         (*this)["filter_beta"] = "0.99999";
         (*this)["filter_gamma"] = "1e-8";
         (*this)["switching_delta"] = "1";
         (*this)["filter_ubd"] = "1e4";
         (*this)["filter_fact"] = "1e4";
         (*this)["filter_switching_infeasibility_exponent"] = "1.1";
         (*this)["armijo_decrease_fraction"] = "1e-8";
         (*this)["LS_backtracking_ratio"] = "0.5";
         (*this)["LS_min_step_length"] = "5e-7";
         (*this)["barrier_tau_min"] = "0.99";
         (*this)["barrier_damping_factor"] = "1e-5";
         (*this)["l1_constraint_violation_coefficient"] = "1000.";
         (*this)["progress_norm"] = "L1";
         (*this)["residual_norm"] = "INF";
         (*this)["scale_functions"] = "yes";
         (*this)["sparse_format"] = "COO";
         (*this)["tolerance"] = "1e-8";
         (*this)["loose_tolerance"] = "1e-6";
         (*this)["loose_tolerance_consecutive_iteration_threshold"] = "15";
         (*this)["switch_to_optimality_requires_linearized_feasibility"] = "no";
         (*this)["LS_scale_duals_with_step_length"] = "yes";
         (*this)["protect_actual_reduction_against_roundoff"] = "yes";
      }
      else if (preset_name == "filtersqp") {
         (*this)["constraint_relaxation_strategy"] = "feasibility_restoration";
         (*this)["subproblem"] = "QP";
         (*this)["globalization_mechanism"] = "TR";
         (*this)["globalization_strategy"] = "fletcher_filter_method";
         (*this)["filter_type"] = "standard";
         (*this)["progress_norm"] = "L1";
         (*this)["residual_norm"] = "L2";
         (*this)["sparse_format"] = "CSC";
         (*this)["TR_radius"] = "10";
         (*this)["l1_constraint_violation_coefficient"] = "1.";
         (*this)["enforce_linear_constraints"] = "yes";
         (*this)["tolerance"] = "1e-6";
         (*this)["loose_tolerance"] = "1e-6";
         (*this)["TR_min_radius"] = "1e-8";
         (*this)["switch_to_optimality_requires_linearized_feasibility"] = "yes";
         (*this)["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "byrd") {
         (*this)["constraint_relaxation_strategy"] = "l1_relaxation";
         (*this)["subproblem"] = "QP";
         (*this)["globalization_mechanism"] = "LS";
         (*this)["globalization_strategy"] = "l1_merit";
         (*this)["l1_relaxation_initial_parameter"] = "1";
         (*this)["LS_backtracking_ratio"] = "0.5";
         (*this)["armijo_decrease_fraction"] = "1e-8";
         (*this)["l1_relaxation_epsilon1"] = "0.1";
         (*this)["l1_relaxation_epsilon2"] = "0.1";
         (*this)["l1_constraint_violation_coefficient"] = "1.";
         (*this)["tolerance"] = "1e-6";
         (*this)["loose_tolerance"] = "1e-6";
         (*this)["progress_norm"] = "L1";
         (*this)["residual_norm"] = "L1";
         (*this)["sparse_format"] = "CSC";
         (*this)["LS_scale_duals_with_step_length"] = "no";
         (*this)["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "funnelsqp") {
         (*this)["constraint_relaxation_strategy"] = "feasibility_restoration";
         (*this)["subproblem"] = "QP";
         (*this)["globalization_mechanism"] = "TR";
         (*this)["globalization_strategy"] = "funnel_method";
         (*this)["progress_norm"] = "L1";
         (*this)["residual_norm"] = "L2";
         (*this)["sparse_format"] = "CSC";
         (*this)["TR_radius"] = "10";
         (*this)["l1_constraint_violation_coefficient"] = "1.";
         (*this)["enforce_linear_constraints"] = "yes";
         (*this)["tolerance"] = "1e-6";
         (*this)["loose_tolerance"] = "1e-6";
         (*this)["TR_min_radius"] = "1e-8";
         (*this)["switch_to_optimality_requires_acceptance"] = "no";
         (*this)["switch_to_optimality_requires_linearized_feasibility"] = "yes";

         (*this)["funnel_beta"] = "0.9999";
         (*this)["funnel_gamma"] = "0.001";
         (*this)["switching_delta"] = "0.999";
         (*this)["funnel_kappa"] = "0.5";
         (*this)["funnel_ubd"] = "1.0";
         (*this)["funnel_fact"] = "1.5";
         (*this)["funnel_switching_infeasibility_exponent"] = "2";
         (*this)["funnel_update_strategy"] = "2";
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known.");
      }
   }


} // namespace
