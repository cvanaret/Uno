// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <fstream>
#include <sstream>
#include "Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // setter
   void Options::set(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten) {
      this->string_options[option_name] = option_value;
      this->overwritten_options[option_name] = flag_as_overwritten;
   }

   void Options::overwrite_with(const Options& overwriting_options) {
      for (const auto& [option_name, option_value]: overwriting_options) {
         // if the option already exists and is not the same, flag it as overwritten
         const auto existing_option_value = this->at_optional(option_name);
         bool flag_as_overwritten = (existing_option_value.has_value() && *existing_option_value != option_value);
         this->set(option_name, option_value, flag_as_overwritten);
      }
   }

   // getters
   const std::string& Options::at(const std::string& option_name) const {
      try {
         const std::string& option_value = this->string_options.at(option_name);
         this->used[option_name] = true;
         return option_value;
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option with name " + option_name + " was not found");
      }
   }

   std::optional<std::string> Options::at_optional(const std::string& option_name) const {
      try {
         const std::string& option_value = this->string_options.at(option_name);
         this->used[option_name] = true;
         return option_value;
      }
      catch(const std::out_of_range&) {
         return std::nullopt;
      }
   }

   const std::string& Options::get_string(const std::string& option_name) const {
      return this->at(option_name);
   }

   std::optional<std::string> Options::get_string_optional(const std::string& option_name) const {
      return this->at_optional(option_name);
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

   // argv[i] for i = offset..argc-1 are overwriting options
   Options Options::get_command_line_options(int argc, char* argv[], size_t offset) {
      static const std::string delimiter = "=";
      Options overwriting_options;

      // build the (name, value) map
      for (size_t i = offset; i < static_cast<size_t>(argc); ++i) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter);
         }
         const std::string option_name = argument.substr(0, position);
         const std::string option_value = argument.substr(position + 1);
         overwriting_options.set(option_name, option_value);
      }
      return overwriting_options;
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
               options.set(option_name, option_value);
            }
         }
         file.close();
      }
      return options;
   }

   void Options::print_used_overwritten() const {
      size_t number_used_options = 0;
      std::string option_list{};
      for (const auto& [option_name, option_value]: this->string_options) {
         if (this->used[option_name] && this->overwritten_options[option_name]) {
            ++number_used_options;
            option_list.append("- ").append(option_name).append(" = ").append(option_value).append("\n");
         }
      }
      // print the overwritten options
      if (number_used_options > 0) {
         DISCRETE << "\nUsed overwritten options:\n" << option_list << '\n';
      }
   }

   std::map<std::string, std::string>::const_iterator Options::begin() const {
      return this->string_options.begin();
   }

   std::map<std::string, std::string>::const_iterator Options::end() const {
      return this->string_options.end();
   }
} // namespace