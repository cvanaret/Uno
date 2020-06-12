#ifndef STATISTICS_H
#define STATISTICS_H

#include <string>
#include <map>

class Statistics {
public:
    static std::map<std::string, std::string> symbols;
    static int int_width;
    static int double_width;
    
    void add_column(std::string name, int width, int order);
    void add_statistic(std::string name, std::string value);
    void add_statistic(std::string name, int value);
    void add_statistic(std::string name, double value);
    void print_header(bool first_occurrence);
    void print_current_line();
    void print_footer();
    void new_line();
    
private:
    std::map<int, std::string> columns_;
    std::map<std::string, int> widths_;
    std::map<std::string, std::string> current_line_;
};

#endif // STATISTICS_H