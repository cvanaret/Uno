#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <string>
#include <map>

class Statistics {
public:
   Statistics() = default;

   static std::map<std::string, std::string> symbols;
   static int int_width;
   static int double_width;
   static int char_width;

   void add_column(std::string name, int width, int order);
   void add_statistic(std::string name, std::string value);
   void add_statistic(std::string name, int value);
   void add_statistic(std::string name, size_t value);
   void add_statistic(std::string name, double value);
   void print_header(bool first_occurrence);
   void print_current_line();
   void print_footer();
   void new_line();

private:
   size_t iteration{0};
   std::map<int, std::string> columns{};
   std::map<std::string, size_t> widths{};
   std::map<std::string, std::string> current_line{};
};

#endif // UNO_STATISTICS_H