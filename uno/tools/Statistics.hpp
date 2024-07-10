// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <string_view>
#include <map>

// forward declaration
class Options;

class Statistics {
public:
   explicit Statistics(const Options& options);

   static int int_width;
   static int double_width;
   static int string_width;
   static int numerical_format_size;

   void add_column(std::string_view name, int width, int order);
   void set(std::string_view name, std::string value);
   void set(std::string_view name, int value);
   void set(std::string_view name, size_t value);
   void set(std::string_view name, double value);
   void print_header(bool first_occurrence);
   void print_current_line();
   void print_footer();
   void start_new_line();

private:
   size_t iteration{0};
   std::map<int, std::string> columns{};
   std::map<std::string_view, int> widths{};
   std::map<std::string_view, std::string> current_line{};

   size_t print_header_every_iterations{};
   static std::string_view symbol(std::string_view value);
};

#endif // UNO_STATISTICS_H