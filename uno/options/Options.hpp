// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <unordered_map>
#include <optional>
#include <string>
#include <vector>
#include "../interfaces/C/uno_int.h"

namespace uno {
   enum class OptionType {INTEGER, DOUBLE, BOOL, STRING};

   class Options {
   public:
      Options() = default;

      void set_integer(const std::string& option_name, uno_int option_value, bool flag_as_overwritten = false);
      void set_double(const std::string& option_name, double option_value, bool flag_as_overwritten = false);
      void set_bool(const std::string& option_name, bool option_value, bool flag_as_overwritten = false);
      void set_string(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten = false);
      // setter for option with unknown type
      void set(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten = false);

      [[nodiscard]] uno_int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;
      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> get_string_optional(const std::string& option_name) const;
      [[nodiscard]] OptionType get_option_type(const std::string& option_name) const;
      [[nodiscard]] std::map<std::string, int32_t> get_integer_options() const;
      [[nodiscard]] std::map<std::string, double> get_double_options() const;
      [[nodiscard]] std::map<std::string, bool> get_bool_options() const;
      [[nodiscard]] std::map<std::string, std::string> get_string_options() const;

      [[nodiscard]] static std::vector<std::pair<std::string, std::string>> get_command_line_options(int argc, char* argv[],
         size_t offset);
      static void load_option_file(Options& options, const std::string& file_name);

      // Print all available options with their type and default value
      static void dump_default_options();

      void print_non_default() const;

      static const std::unordered_map<std::string, OptionType> option_types;

   private:
      std::unordered_map<std::string, uno_int> integer_options{};
      std::unordered_map<std::string, double> double_options{};
      std::unordered_map<std::string, bool> bool_options{};
      std::unordered_map<std::string, std::string> string_options{};

      mutable std::unordered_map<std::string, bool> used{};
      mutable std::unordered_map<std::string, bool> overwritten_options{};
   };
} // namespace

#endif // UNO_OPTIONS_H