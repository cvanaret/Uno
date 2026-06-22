// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <map>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>

namespace uno {
   class Statistics {
   public:
      Statistics() = default;
      static size_t int_width, double_width, string_width;

      void add_column(std::string_view name, size_t width, size_t precision);
      void start_new_line();

      void set(std::string_view name, std::string value);
      void set(std::string_view name, int value);
      void set(std::string_view name, size_t value);
      void set(std::string_view name, double value);

      void print_horizontal_line();
      void print_header();
      void print_current_line();
      void print_footer();

   protected:
      // name -> index mapping
      static const std::map<std::string_view, size_t> column_order;

      // one contiguous record per column, ordered by appearance (already sorted by index)
      struct Column {
         std::string_view name;
         size_t width;
         size_t precision;
         std::string value; // current-line cell, empty if !is_set
         bool is_set = false;
      };

      std::vector<Column> columns; // index = print order
      std::unordered_map<std::string_view, size_t> name_to_index; // name to position in `columns`
      bool finalized = false; // columns sorted and name_to_index built?

      // sort `columns` by index and build name_to_index, called after the last add_column()
      void finalize();
      [[nodiscard]] size_t index_of(std::string_view name);
      void set_value(size_t index, std::string_view value);
   };
}

#endif // UNO_STATISTICS_H