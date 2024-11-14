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
      explicit Options(bool are_default_options);

      [[nodiscard]] size_t size() const;
      std::string& operator[](const std::string& option_name);

      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> get_string_optional(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;

      [[nodiscard]] static Options get_command_line_options(int argc, char* argv[], size_t offset);
      [[nodiscard]] static Options load_option_file(const std::string& file_name);
      void overwrite_with(const Options& overwriting_options);
      void print_used() const;

      [[nodiscard]] std::map<std::string, std::string>::const_iterator begin() const;
      [[nodiscard]] std::map<std::string, std::string>::const_iterator end() const;

   private:
      std::map<std::string, std::string> options{};
      mutable std::map<std::string, bool> used{};
      mutable std::map<std::string, bool> is_default{};
      const bool are_default_options;

      [[nodiscard]] const std::string& at(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> at_optional(const std::string& option_name) const;
   };
} // namespace

#endif // UNO_OPTIONS_H
