// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <sstream>
#include <iostream>
#include <iomanip>
#include "Statistics.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // TODO move this to the option file
   size_t Statistics::int_width = 7;
   size_t Statistics::double_width = 10;
   size_t Statistics::string_width = 15;

   const std::map<std::string_view, size_t> Statistics::column_order = {
      {"Major", 1},
      {"Minor", 2},
      {"Penalty", 5},
      {"Barrier", 8},
      {"Steplength", 10},
      {"Radius", 11},
      {"Phase", 20},
      {"Regulariz", 21},
      {"Funnel", 25},
      {"||Step||", 31},
      {"Objective", 100},
      {"Infeas", 101},
      {"Statio", 104},
      {"Compl", 105},
      {"Status", 200},
   };

   void Statistics::add_column(std::string_view name, size_t width, size_t precision, size_t order) {
      this->columns[order] = name;
      this->widths[name] = width;
      this->precisions[name] = precision;
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
      // get the associated precision
      const size_t precision = this->precisions.at(name);
      stream << std::scientific << std::setprecision(static_cast<int>(precision)) << value;
      this->set(name, stream.str());
   }
   
   void Statistics::print_horizontal_line() {
      for (const auto& element: this->columns) {
         std::string header = element.second;
         for (size_t j = 0; j < this->widths[header]; j++) {
            INFO << Statistics::symbol("top");
         }
      }
      INFO << '\n';
   }

   void Statistics::print_header() {
      // line above
      this->print_horizontal_line();
      // headers
      INFO << "  Iterations\n";
      for (const auto& element: this->columns) {
         const std::string& header = element.second;
         INFO << " " << header;
         const size_t number_spaces = (this->widths[header] >= header.size() + 1) ?
            this->widths[header] - (header.size() + 1) : static_cast<size_t>(0);
         for (size_t j = 0; j < number_spaces; j++) {
            INFO << " ";
         }
      }
      INFO << '\n';
      // line below
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
      for (const auto& element: this->columns) {
         const auto& header = element.second;
         size_t length;
         try {
            const auto& value = this->current_line.at(header);
            INFO << " " << value;
            length = 1 + length_utf8(value);
         }
         catch (const std::out_of_range&) {
            INFO << " -";
            length = 2;
         }
         const size_t number_spaces = (length <= this->widths[header]) ? this->widths[header] - length : static_cast<size_t>(0);
         for (size_t j = 0; j < number_spaces; j++) {
            INFO << " ";
         }
      }
      INFO << '\n';
   }

   void Statistics::print_footer() {
      /*
      for (const auto& element: this->columns) {
         const auto& header = element.second;
         for (size_t j = 0; j < this->widths[header]; j++) {
            INFO << Statistics::symbol("bottom");
         }
      }
      INFO << '\n';
      */
      Statistics::print_header();
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
