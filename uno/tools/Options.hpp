// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <map>
#include <string>

namespace uno {
   class Options {
   public:
      Options() = default;

      [[nodiscard]] size_t size() const;
      std::string& operator[](const std::string& option_name);

      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;

      static Options get_default();
      static Options get_command_line_options(int argc, char* argv[]);
      static void overwrite_with_option_file(Options& options, const std::string& file_name);
      static void set_preset(Options& options, const std::string& preset_name);
      void overwrite_with(const Options& overwriting_options);

      std::map<std::string, std::string>::const_iterator begin() const;
      std::map<std::string, std::string>::const_iterator end() const;

   private:
      std::map<std::string, std::string> options{};

      [[nodiscard]] const std::string& at(const std::string& option_name) const;
   };
} // namespace

#endif // UNO_OPTIONS_H
