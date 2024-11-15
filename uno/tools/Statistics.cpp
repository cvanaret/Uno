// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <sstream>
#include <iostream>
#include <iomanip>
#include "Statistics.hpp"
#include "options/Options.hpp"

namespace uno {
   // TODO move this to the option file
   int Statistics::int_width = 7;
   int Statistics::double_width = 17;
   int Statistics::string_width = 26;
   int Statistics::numerical_format_size = 4;

   Statistics::Statistics(const Options& options): print_header_frequency(options.get_unsigned_int("statistics_print_header_frequency")) { }

   void Statistics::add_column(std::string_view name, int width, int order) {
      this->columns[order] = name;
      this->widths[name] = width;
   }
   
   void Statistics::start_new_line() {
      for (const auto& column: this->columns) {
         this->current_line[column.second] = "-";
      }
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
      stream << std::scientific << std::setprecision(Statistics::numerical_format_size) << value;
      this->set(name, stream.str());
   }
   
   void Statistics::print_horizontal_line() {
      for (const auto& element: this->columns) {
         std::string header = element.second;
         for (int j = 0; j < this->widths[header]; j++) {
            std::cout << Statistics::symbol("top");
         }
      }
      std::cout << '\n';
   }

   void Statistics::print_header() {
      /* line above */
      this->print_horizontal_line();
      /* headers */
      for (const auto& element: this->columns) {
         const std::string& header = element.second;
         std::cout << " " << header;
         for (int j = 0; j < this->widths[header] - static_cast<int>(header.size()) - 1; j++) {
            std::cout << " ";
         }
      }
      std::cout << '\n';
      /* line below */
      this->print_horizontal_line();
   }

   // https://cplusplus.com/forum/beginner/192031/
   std::size_t length_utf8(const std::string& str) {
      size_t length = 0;
      for (char c: str) {
         if ((c & 0xC0) != 0x80) {
            ++length;
         }
      }
      return length;
   }

   void Statistics::print_current_line() {
      if (this->line_number % this->print_header_frequency == 0) {
         this->print_header();
      }
      for (const auto& element: this->columns) {
         const auto& header = element.second;
         int length;
         try {
            const auto& value = this->current_line.at(header);
            std::cout << " " << value;
            length = 1 + static_cast<int>(length_utf8(value));
         }
         catch (const std::out_of_range&) {
            std::cout << " -";
            length = 2;
         }
         int number_spaces = (length <= this->widths[header]) ? this->widths[header] - length : 0;
         for (int j = 0; j < number_spaces; j++) {
            std::cout << " ";
         }
      }
      std::cout << '\n';
      this->line_number++;
   }

   void Statistics::print_footer() {
      for (const auto& element: this->columns) {
         const auto& header = element.second;
         for (int j = 0; j < this->widths[header]; j++) {
            std::cout << Statistics::symbol("bottom");
         }
      }
      std::cout << '\n';
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
} // namespace
