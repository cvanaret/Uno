// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <fstream>
#include <sstream>
#include "Options.hpp"
#include "Logger.hpp"

std::string& Options::operator[](const std::string& key) {
   return this->options[key];
}

const std::string& Options::at(const std::string& key) const {
   try {
      return this->options.at(key);
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

void Options::print() const {
   std::cout << "Options:\n";
   for (const auto& [key, value]: this->options) {
      std::cout << "- " << key << " = " << value << '\n';
   }
}

Options get_default_options(const std::string& file_name) {
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

void find_preset(const std::string& preset_name, Options& options) {
   // shortcuts for state-of-the-art combinations
   if (preset_name == "ipopt") {
      options["constraint_relaxation_strategy"] = "feasibility_restoration";
      options["subproblem"] = "primal_dual_interior_point";
      options["globalization_mechanism"] = "LS";
      options["globalization_strategy"] = "waechter_filter_strategy";
      options["filter_type"] = "standard";
      options["filter_beta"] = "0.99999";
      options["filter_gamma"] = "1e-8";
      options["filter_delta"] = "1";
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
   }
   else if (preset_name == "filtersqp") {
      options["constraint_relaxation_strategy"] = "feasibility_restoration";
      options["subproblem"] = "QP";
      options["globalization_mechanism"] = "TR";
      options["globalization_strategy"] = "leyffer_filter_strategy";
      options["filter_type"] = "standard";
      options["progress_norm"] = "L1";
      options["residual_norm"] = "L2";
      options["sparse_format"] = "CSC";
      options["TR_radius"] = "10";
      options["l1_constraint_violation_coefficient"] = "1.";
      options["enforce_linear_constraints"] = "yes";
      options["tolerance"] = "1e-6";
      options["TR_min_radius"] = "1e-8";
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
      options["progress_norm"] = "L1";
      options["residual_norm"] = "L1";
      options["sparse_format"] = "CSC";
   }
}

void get_command_line_options(int argc, char* argv[], Options& options) {
   // build the (name, value) map
   int i = 1;
   while (i < argc - 1) {
      std::string argument = std::string(argv[i]);
      if (argument[0] == '-') {
         if (i < argc - 1) {
            // remove the '-'
            const std::string name = argument.substr(1);
            const std::string value = std::string(argv[i + 1]);
            if (name == "preset") {
               find_preset(value, options);
            }
            else {
               options[name] = value;
            }
            i += 2;
         }
      }
      else {
         WARNING << "Argument " << argument << " was ignored\n";
         i++;
      }
   }
}
