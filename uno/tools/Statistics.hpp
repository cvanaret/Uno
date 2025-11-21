// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <string_view>
#include <map>

namespace uno {
   class Statistics {
   public:
      Statistics() = default;

      static size_t int_width;
      static size_t double_width;
      static size_t string_width;
      static size_t numerical_format_size;

      void add_column(std::string_view name, size_t width, size_t precision, size_t order);
      void start_new_line();
      void set(std::string_view name, std::string value);
      void set(std::string_view name, int value);
      void set(std::string_view name, size_t value);
      void set(std::string_view name, double value);
      
      void print_horizontal_line();
      void print_header();
      void print_current_line();
      void print_footer();

      static const std::map<std::string_view, size_t> column_order;

   private:
      std::map<size_t, std::string> columns{};
      std::map<std::string_view, size_t> widths{};
      std::map<std::string_view, size_t> precisions{};
      std::map<std::string_view, std::string> current_line{};
      static std::string_view symbol(std::string_view value);
   };
} // namespace

#endif // UNO_STATISTICS_H
