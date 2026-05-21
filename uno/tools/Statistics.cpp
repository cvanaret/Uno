// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <charconv>
#include <algorithm>
#include "Statistics.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   size_t Statistics::int_width = 7;
   size_t Statistics::double_width = 10;
   size_t Statistics::string_width = 15;

   const std::map<std::string_view, size_t> Statistics::column_order = {
      {"Major", 1}, {"Minor", 2}, {"Penalty", 5}, {"Barrier", 8},
      {"Steplength", 10}, {"Radius", 11}, {"Phase", 20}, {"Regulariz", 21},
      {"Funnel", 25}, {"|BFGS|", 26}, {"|SR1|", 26}, {"||Step||", 31},
      {"Objective", 100}, {"Infeas", 101}, {"Statio", 104}, {"Compl", 105},
      {"Status", 200},
   };

   // number of glyphs in a UTF-8 string (counts non-continuation bytes)
   static size_t length_utf8(std::string_view string) {
      size_t length = 0;
      for (unsigned char c: string) {
         if ((c & 0xC0) != 0x80) {
            ++length;
         }
      }
      return length;
   }

   // append leading space, `value`, then pad with spaces so the field occupies `width` glyphs
   static void append_cell(std::string& line, std::string_view value, size_t width) {
      line.push_back(' ');
      line.append(value);
      const size_t used = 1 + length_utf8(value);
      if (used < width) {
         line.append(width - used, ' ');
      }
   }

   size_t Statistics::index_of(std::string_view name) {
      if (!this->finalized) {
         this->finalize();
      }
      return this->name_to_index.at(name);
   }

   void Statistics::add_column(std::string_view name, size_t width, size_t precision) {
      // validate and append the column
      if (column_order.find(name) == column_order.end()) {
         throw std::invalid_argument("The column " + std::string(name) + " does not exist");
      }
      this->columns.push_back(Column{name, width, precision, std::string{}, false});
      // the name->index mapping is built in finalize()
      this->finalized = false;
   }

   void Statistics::finalize() {
      // order columns by their index
      std::sort(this->columns.begin(), this->columns.end(),
         [](const Column& a, const Column& b) {
            return column_order.at(a.name) < column_order.at(b.name);
         });
      // build name->index mapping
      this->name_to_index.clear();
      this->name_to_index.reserve(this->columns.size());
      for (size_t column_index: Range(this->columns.size())) {
         this->name_to_index.emplace(this->columns[column_index].name, column_index);
      }
      this->finalized = true;
   }

   void Statistics::start_new_line() {
      for (auto& column : this->columns) {
         column.is_set = false;
         column.value.clear(); // keeps capacity
      }
   }

   void Statistics::set_value(size_t index, std::string_view value) {
      Column& column = this->columns[index];
      column.value.assign(value);
      column.is_set = true;
   }

   void Statistics::set(std::string_view name, std::string value) {
      Column& column = this->columns[this->index_of(name)];
      column.value = std::move(value);
      column.is_set = true;
   }

   void Statistics::set(std::string_view name, int value) {
      char buffer[16];
      auto [pointer, _] = std::to_chars(buffer, buffer + sizeof(buffer), value);
      this->set_value(this->index_of(name), std::string_view(buffer, static_cast<size_t>(pointer - buffer)));
   }

   void Statistics::set(std::string_view name, size_t value) {
      char buffer[24];
      auto [pointer, _] = std::to_chars(buffer, buffer + sizeof(buffer), value);
      this->set_value(this->index_of(name), std::string_view(buffer, static_cast<size_t>(pointer - buffer)));
   }

   void Statistics::set(std::string_view name, double value) {
      const size_t index = this->index_of(name);
      char buffer[40];
      const int precision = static_cast<int>(this->columns[index].precision);
      auto [pointer, _] = std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::scientific, precision);
      this->set_value(index, std::string_view(buffer, static_cast<size_t>(pointer - buffer)));
   }

   void Statistics::print_horizontal_line() {
      if (!this->finalized) {
         this->finalize();
      }
      // count the glyphs (column-width units)
      size_t glyph_count = 0;
      for (const auto& column : this->columns) {
         glyph_count += column.width;
      }
      // create the horizontal line
      std::string line;
      line.reserve(glyph_count * top_symbol.size() + 1);
      for (size_t i = 0; i < glyph_count; ++i) {
         line.append(top_symbol);
      }
      line.push_back('\n');
      INFO << line; // single insertion
   }

   void Statistics::print_header() {
      if (!this->finalized) {
         this->finalize();
      }
      this->print_horizontal_line();
      INFO << "  Iterations\n";

      std::string line;
      line.reserve(128);
      for (const auto& column : this->columns) {
         append_cell(line, column.name, column.width);
      }
      line.push_back('\n');
      INFO << line; // single insertion

      this->print_horizontal_line();
   }

   void Statistics::print_current_line() {
      if (!this->finalized) {
         this->finalize();
      }
      std::string line;
      line.reserve(128);
      for (const auto& column : this->columns) {
         const std::string_view cell = column.is_set ? std::string_view(column.value) : "-";
         append_cell(line, cell, column.width);
      }
      line.push_back('\n');
      INFO << line; // single insertion
   }

   void Statistics::print_footer() {
      this->print_header();
   }
}