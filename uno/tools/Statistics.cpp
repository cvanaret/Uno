// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <iomanip>
#include "Statistics.hpp"
#include "Options.hpp"

// TODO move this to the option file
int Statistics::int_width = 7;
int Statistics::double_width = 17;
int Statistics::string_width = 26;
int Statistics::numerical_format_size = 6;

Statistics::Statistics(const Options& options): print_header_every_iterations(options.get_unsigned_int("statistics_print_header_every_iterations")) {
}

void Statistics::add_column(std::string_view name, int width, int order) {
   this->columns[order] = name;
   this->widths[name] = width;
}

void Statistics::set(std::string_view name, std::string value) {
   this->current_line[name] = std::move(value);
}

void Statistics::set(std::string_view name, int value) {
   this->set(name, std::to_string(value));
}

void Statistics::set(std::string_view name, size_t value) {
   this->set(name, std::to_string(value));
}

void Statistics::set(std::string_view name, double value) {
   std::ostringstream stream;
   stream << std::defaultfloat << std::setprecision(Statistics::numerical_format_size) << value;
   this->set(name, stream.str());
}

void Statistics::print_header(bool first_occurrence) {
   /* line above */
   std::cout << (first_occurrence ? Statistics::symbol("top-left") : Statistics::symbol("left-mid"));
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << (first_occurrence ? Statistics::symbol("top-mid") : Statistics::symbol("mid-mid"));
      }
      const std::string& header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("top");
      }
      k++;
   }
   std::cout << (first_occurrence ? Statistics::symbol("top-right") : Statistics::symbol("right-mid")) << '\n';
   /* headers */
   std::cout << Statistics::symbol("left");
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("middle");
      }
      const std::string& header = element.second;
      std::cout << " " << header;
      for (int j = 0; j < this->widths[header] - static_cast<int>(header.size()) - 1; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol("right") << '\n';
}

void Statistics::print_current_line() {
   if (this->iteration % this->print_header_every_iterations == 0) {
      this->print_header(this->iteration == 0);
   }
   std::cout << Statistics::symbol("left-mid");
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("mid-mid");
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("bottom");
      }
      k++;
   }
   std::cout << Statistics::symbol("right-mid") << '\n';
   /* headers */
   std::cout << Statistics::symbol("left");
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("middle");
      }
      const std::string& header = element.second;
      int size;
      try {
         std::string value = this->current_line.at(header);
         std::cout << " " << value;
         size = 1 + static_cast<int>(value.size());
      }
      catch (const std::out_of_range&) {
         std::cout << " -";
         size = 2;
      }
      int number_spaces = (size <= this->widths[header]) ? this->widths[header] - size : 0;
      for (int j = 0; j < number_spaces; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol("right") << '\n';
   this->iteration++;
}

void Statistics::print_footer() {
   std::cout << Statistics::symbol("bottom-left");
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("bottom-mid");
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("bottom");
      }
      k++;
   }
   std::cout << Statistics::symbol("bottom-right") << '\n';
}

void Statistics::start_new_line() {
   //this->current_line.clear();
   for (const auto& column: this->columns) {
      this->current_line[column.second] = "-";
   }
}

std::string_view Statistics::symbol(std::string_view value) {
   static std::map<std::string_view, std::string_view> symbols = {
         {"top", "─"},
         {"top-mid", "┬"},
         {"top-left", "┌"},
         {"top-right", "┐"},
         {"bottom", "─"},
         {"bottom-mid", "┴"},
         {"bottom-left", "└"},
         {"bottom-right", "┘"},
         {"left", "│"},
         {"left-mid", "├"},
         {"mid", "─"},
         {"mid-mid", "┼"},
         {"right", "│"},
         {"right-mid", "┤"},
         {"middle", "│"}
   };
   return symbols[value];
}