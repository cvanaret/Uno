// Copyright (c) 2022 Charlie Vanaret
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
      options["mechanism"] = "LS";
      options["constraint-relaxation"] = "feasibility-restoration";
      options["strategy"] = "waechter-filter-strategy";
      options["filter_type"] = "standard";
      options["subproblem"] = "barrier";
      options["filter_beta"] = "0.99999";
      options["filter_gamma"] = "1e-8";
      options["filter_delta"] = "1";
      options["filter_ubd"] = "1e4";
      options["filter_fact"] = "1e4";
      options["filter_switching_infeasibility_exponent"] = "1.1";
      options["armijo_decrease_fraction"] = "1e-8";
      options["LS_backtracking_ratio"] = "0.5";
      options["barrier_tau_min"] = "0.99";
      options["use_second_order_correction"] = "yes";
      options["residual_norm"] = "INF";
      options["scale_functions"] = "yes";
      options["sparse_format"] = "COO";
      options["tolerance"] = "1e-8";
   }
   else if (preset_name == "filtersqp") {
      options["mechanism"] = "TR";
      options["constraint-relaxation"] = "feasibility-restoration";
      options["strategy"] = "leyffer-filter-strategy";
      options["filter_type"] = "standard";
      options["subproblem"] = "QP";
      options["residual_norm"] = "L1";
      options["sparse_format"] = "CSC";
      options["TR_radius"] = "10";
      options["enforce_linear_constraints"] = "yes";
      options["tolerance"] = "1e-6";
   }
   else if (preset_name == "byrd") {
      options["mechanism"] = "LS";
      options["constraint-relaxation"] = "l1-relaxation";
      options["strategy"] = "merit";
      options["subproblem"] = "QP";
      options["l1_relaxation_initial_parameter"] = "1";
      options["LS_backtracking_ratio"] = "0.5";
      options["armijo_decrease_fraction"] = "1e-8";
      options["l1_relaxation_epsilon1"] = "0.1";
      options["l1_relaxation_epsilon2"] = "0.1";
      options["tolerance"] = "1e-6";
      options["residual_norm"] = "L1";
      options["sparse_format"] = "CSC";
   }
}

void get_command_line_options(int argc, char* argv[], Options& options) {
   // build the (argument, value) map
   int i = 1;
   while (i < argc - 1) {
      std::string argument_i = std::string(argv[i]);
      if (argument_i[0] == '-') {
         if (i < argc - 1) {
            // remove the '-'
            const std::string name = argument_i.substr(1);
            const std::string value_i = std::string(argv[i + 1]);
            if (name == "preset") {
               find_preset(value_i, options);
            }
            else {
               options[name] = value_i;
            }
            i += 2;
         }
      }
      else {
         WARNING << "Argument " << argument_i << " was ignored\n";
         i++;
      }
   }
}

void set_logger(const std::string& logger_level) {
   try {
      if (logger_level == "ERROR") {
         Logger::logger_level = ERROR;
      }
      else if (logger_level == "WARNING") {
         Logger::logger_level = WARNING;
      }
      else if (logger_level == "INFO") {
         Logger::logger_level = INFO;
      }
      else if (logger_level == "DEBUG") {
         Logger::logger_level = DEBUG;
      }
   }
   catch (const std::out_of_range&) {
      throw std::out_of_range("The logger level " + logger_level + " was not found");
   }
}
