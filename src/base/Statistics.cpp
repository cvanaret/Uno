#include <iostream>
#include "Statistics.hpp"

std::map<std::string, std::string> Statistics::symbols = {
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
int Statistics::double_width = 21;

void Statistics::add_column(std::string name, int width, int order) {
    this->columns_[order] = name;
    this->widths_[name] = width;
    return;
}

void Statistics::add_statistic(std::string name, std::string value) {
    this->current_line_[name] = value;
    return;
}

void Statistics::add_statistic(std::string name, int value) {
    return add_statistic(name, std::to_string(value));
}

void Statistics::add_statistic(std::string name, double value) {
    return add_statistic(name, std::to_string(value));
}

void Statistics::print_header(bool first_occurrence) {
    /* line above */
    std::cout << (first_occurrence ? Statistics::symbols["top-left"] : Statistics::symbols["left-mid"]);
    int k = 0;
    for (const std::pair<int, std::string>& element: this->columns_) {
        if (0 < k) {
            std::cout << (first_occurrence ? Statistics::symbols["top-mid"] : Statistics::symbols["mid-mid"]);
        }
        std::string header = element.second;
        for (int j = 0; j < this->widths_[header]; j++) {
            std::cout << Statistics::symbols["top"];
        }
        k++;
    }
    std::cout << (first_occurrence ? Statistics::symbols["top-right"] : Statistics::symbols["right-mid"]) << "\n";
    /* headers */
    std::cout << Statistics::symbols["left"];
    k = 0;
    for (const std::pair<int, std::string>& element: this->columns_) {
        if (0 < k) {
            std::cout << Statistics::symbols["middle"];
        }
        std::string header = element.second;
        std::cout << " " << header;
        for (unsigned int j = 0; j < this->widths_[header] - header.size() - 1; j++) {
            std::cout << " ";
        }
        k++;
    }
    std::cout << Statistics::symbols["right"] << "\n";
}

void Statistics::print_current_line() {
    std::cout << Statistics::symbols["left-mid"];
    int k = 0;
    for (const std::pair<int, std::string>& element: this->columns_) {
        if (0 < k) {
            std::cout << Statistics::symbols["mid-mid"];
        }
        std::string header = element.second;
        for (int j = 0; j < this->widths_[header]; j++) {
            std::cout << Statistics::symbols["bottom"];
        }
        k++;
    }
    std::cout << Statistics::symbols["right-mid"] << "\n";
    /* headers */
    std::cout << Statistics::symbols["left"];
    k = 0;
    for (const std::pair<int, std::string>& element: this->columns_) {
        if (0 < k) {
            std::cout << Statistics::symbols["middle"];
        }
        std::string header = element.second;
        unsigned int size;
        try {
            std::string value = this->current_line_.at(header);
            std::cout << " " << value;
            size = 1 + value.size();
        }
        catch (const std::out_of_range) {
            std::cout << " -";
            size = 2;
        }
        for (unsigned int j = 0; j < this->widths_[header] - size; j++) {
            std::cout << " ";
        }
        k++;
    }
    std::cout << Statistics::symbols["right"] << "\n";
}

void Statistics::print_footer() {
    std::cout << Statistics::symbols["bottom-left"];
    int k = 0;
    for (const std::pair<int, std::string>& element: this->columns_) {
        if (0 < k) {
            std::cout << Statistics::symbols["bottom-mid"];
        }
        std::string header = element.second;
        for (int j = 0; j < this->widths_[header]; j++) {
            std::cout << Statistics::symbols["bottom"];
        }
        k++;
    }
    std::cout << Statistics::symbols["bottom-right"] << "\n";
}

void Statistics::new_line() {
    this->current_line_.clear();
    return;
}