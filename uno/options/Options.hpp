// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>
#include "../interfaces/C/uno_int.h"

namespace uno {
   enum class OptionType {INTEGER, DOUBLE, BOOL, STRING};

   using OptionValue = std::variant<uno_int, double, bool, std::string>;

   // errors raised when options are set. They derive from the standard exceptions previously thrown, so that existing
   // handlers keep working
   // no option has this name
   class OptionNotFound: public std::out_of_range {
   public:
      using std::out_of_range::out_of_range;
   };

   // the option exists, but has a different type
   class OptionWrongType: public std::invalid_argument {
   public:
      using std::invalid_argument::invalid_argument;
   };

   // the value is not valid for this option (e.g. not a number, or not one of the possible values of a string option)
   class InvalidOptionValue: public std::invalid_argument {
   public:
      using std::invalid_argument::invalid_argument;
   };

   // the option file cannot be opened or is malformed
   class InvalidOptionFile: public std::invalid_argument {
   public:
      using std::invalid_argument::invalid_argument;
   };

   class OptionOverride {
   public:
      std::string option_name;
      std::optional<std::string> old_value;
      std::string new_value;
      std::string reason;

      OptionOverride(std::string option_name, std::optional<std::string> old_value, std::string new_value, std::string reason);
   };

   class Options {
   public:
      // maximum number of characters of a string option value. Longer values (and values that contain a null character)
      // are rejected by the setters, so that all interfaces (C, Fortran, Julia, Python, AMPL, option files) store and
      // return strings of bounded length. Must match UNO_MAX_STRING_LENGTH in the C API
      static constexpr size_t max_string_length = 4095;

      Options() = default;

      void set_integer(const std::string& option_name, uno_int option_value);
      void set_double(const std::string& option_name, double option_value);
      void set_bool(const std::string& option_name, bool option_value);
      void set_string(const std::string& option_name, std::string option_value);
      // setter for option with unknown type
      void set(const std::string& option_name, std::string option_value);

      void overwrite(const Options& other);

      [[nodiscard]] uno_int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;
      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;

      [[nodiscard]] std::optional<uno_int> get_int_optional(const std::string& option_name) const;
      [[nodiscard]] std::optional<size_t> get_unsigned_int_optional(const std::string& option_name) const;
      [[nodiscard]] std::optional<double> get_double_optional(const std::string& option_name) const;
      [[nodiscard]] std::optional<bool> get_bool_optional(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> get_string_optional(const std::string& option_name) const;

      [[nodiscard]] static OptionType get_option_type(const std::string& option_name);
      // whether a value was set for the option
      [[nodiscard]] bool contains(const std::string& option_name) const;

      // possible values of a string option. Empty if the option accepts any string (e.g. a path)
      [[nodiscard]] static std::vector<std::string> get_possible_values(const std::string& option_name);

      [[nodiscard]] static std::vector<std::pair<std::string, std::string>> get_command_line_options(int argc, char* argv[],
         size_t offset);
      static void load_option_file(Options& options, const std::string& file_name);

      // print the default options
      static void dump_default_options();

      void print(const std::string& header) const;

      static const std::unordered_map<std::string, OptionType> option_types;

   protected:
      std::unordered_map<std::string, OptionValue> values;

      static std::string to_string(const OptionValue& value);
      static void check_option_type(const std::string& option_name, OptionType expected_type);
      static void check_string_value(const std::string& option_name, const std::string& option_value);
      static void check_string_option_value(const std::string& option_name, const std::string& option_value);

      template <typename T>
      static constexpr const char* option_type_name() {
         if constexpr (std::is_same_v<T, uno_int>) return "int";
         else if constexpr (std::is_same_v<T, double>) return "double";
         else if constexpr (std::is_same_v<T, bool>) return "bool";
         else if constexpr (std::is_same_v<T, std::string>) return "string";
         else {
            static_assert(sizeof(T) == 0, "Unsupported option type");
            return "";
         }
      }

      template <typename T>
      const T& get(const std::string& option_name) const {
         try {
            const OptionValue& value = this->values.at(option_name);
            if (const T* typed_value = std::get_if<T>(&value)) {
               return *typed_value;
            }
            throw std::runtime_error("The option " + option_name + " is not a " + option_type_name<T>() + " option");
         }
         catch (const std::out_of_range&) {
            throw std::out_of_range(std::string("The ") + option_type_name<T>() + " option " + option_name + " is not available");
         }
      }

      template <typename T>
      std::optional<T> get_optional(const std::string& option_name) const {
         const auto it = this->values.find(option_name);
         if (it == this->values.end()) {
            return std::nullopt;
         }
         if (const T* typed_value = std::get_if<T>(&it->second)) {
            return *typed_value;
         }
         throw std::runtime_error("The option " + option_name + " is not a " + option_type_name<T>() + " option");
      }
   };
} // namespace

#endif // UNO_OPTIONS_H