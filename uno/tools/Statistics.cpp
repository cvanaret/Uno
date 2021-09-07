#include <iostream>
#include <iomanip>
#include "Statistics.hpp"

std::map <std::string, std::string>Statistics::symbols = {
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

int Statistics::int_width = 7;
int Statistics::double_width = 18;
int Statistics::char_width = 7;

void Statistics::add_column(std::string name, int width, int order) {
   this->columns[order] = name;
   this->widths[std::move(name)] = width;
}

void Statistics::add_statistic(std::string name, std::string value) {
   this->current_line[std::move(name)] = std::move(value);
}

void Statistics::add_statistic(std::string name, int value) {
   add_statistic(std::move(name), std::move(std::to_string(value)));
}

void Statistics::add_statistic(std::string name, double value) {
   std::ostringstream stream;
   stream << std::scientific << std::setprecision(6) << value;
   add_statistic(std::move(name), std::move(stream.str()));
}

void Statistics::print_header(bool first_occurrence) {
   /* line above */
   std::cout << (first_occurrence ? Statistics::symbols["top-left"] : Statistics::symbols["left-mid"]);
   int k = 0;
   for (const std::pair<const int, std::string>& element: this->columns) {
      if (0 < k) {
         std::cout << (first_occurrence ? Statistics::symbols["top-mid"] : Statistics::symbols["mid-mid"]);
      }
      std::string header = element.second;
      for (size_t j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbols["top"];
      }
      k++;
   }
   std::cout << (first_occurrence ? Statistics::symbols["top-right"] : Statistics::symbols["right-mid"]) << "\n";
   /* headers */
   std::cout << Statistics::symbols["left"];
   k = 0;
   for (const std::pair<const int, std::string>& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbols["middle"];
      }
      std::string header = element.second;
      std::cout << " " << header;
      for (size_t j = 0; j < this->widths[header] - header.size() - 1; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbols["right"] << "\n";
}

void Statistics::print_current_line() {
   if (this->iteration % 15 == 0) {
      this->print_header(this->iteration == 0);
   }
   std::cout << Statistics::symbols["left-mid"];
   int k = 0;
   for (const std::pair<const int, std::string>& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbols["mid-mid"];
      }
      std::string header = element.second;
      for (size_t j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbols["bottom"];
      }
      k++;
   }
   std::cout << Statistics::symbols["right-mid"] << "\n";
   /* headers */
   std::cout << Statistics::symbols["left"];
   k = 0;
   for (const std::pair<const int, std::string>& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbols["middle"];
      }
      const std::string& header = element.second;
      size_t size;
      try {
         std::string value = this->current_line.at(header);
         std::cout << " " << value;
         size = 1 + value.size();
      }
      catch (const std::out_of_range&) {
         std::cout << " -";
         size = 2;
      }
      size_t number_spaces = std::max(0, (int) this->widths[header] - (int) size);
      for (size_t j = 0; j < number_spaces; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbols["right"] << "\n";
   this->iteration++;
}

void Statistics::print_footer() {
   std::cout << Statistics::symbols["bottom-left"];
   int k = 0;
   for (const std::pair<const int, std::string>& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbols["bottom-mid"];
      }
      std::string header = element.second;
      for (size_t j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbols["bottom"];
      }
      k++;
   }
   std::cout << Statistics::symbols["bottom-right"] << "\n";
}

void Statistics::new_line() {
   this->current_line.clear();
}
