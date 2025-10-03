// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <map>
#include <string>
#include <optional>

namespace uno {
   class Options {
   public:
      Options() = default;

      void set_integer(const std::string& option_name, int32_t option_value, bool flag_as_overwritten = false);
      void set_double(const std::string& option_name, double option_value, bool flag_as_overwritten = false);
      void set_bool(const std::string& option_name, bool option_value, bool flag_as_overwritten = false);
      void set_string(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten = false);
      void overwrite_with(const Options& overwriting_options);

      [[nodiscard]] int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;
      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> get_string_optional(const std::string& option_name) const;

      [[nodiscard]] static Options get_command_line_options(int argc, char* argv[], size_t offset);
      [[nodiscard]] static Options load_option_file(const std::string& file_name);
      
      void print_used_overwritten() const;

      [[nodiscard]] std::map<std::string, std::string>::const_iterator begin() const;
      [[nodiscard]] std::map<std::string, std::string>::const_iterator end() const;

   private:
      std::map<std::string, int32_t> integer_options{};
      std::map<std::string, double> double_options{};
      std::map<std::string, bool> bool_options{};
      std::map<std::string, std::string> string_options{};

      mutable std::map<std::string, bool> used{};
      mutable std::map<std::string, bool> overwritten_options{};

      [[nodiscard]] const std::string& at(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> at_optional(const std::string& option_name) const;
   };

   template <typename T>
   const T& get(const std::map<std::string, T>& options, const std::string& option_name) {
      return options.at(option_name);
   }

} // namespace

#endif // UNO_OPTIONS_H