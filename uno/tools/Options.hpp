#ifndef OPTIONS_H
#define OPTIONS_H

#include <map>

class Options {
public:
   Options() = default;
   std::string& operator[](const std::string& key);
   [[nodiscard]] const std::string& at(const std::string& key) const;
   void print() const;

private:
   std::map<std::string, std::string> options{};
};

Options get_default_options(const std::string& file_name);
void get_command_line_options(int argc, char* argv[], Options& options);
void set_logger(const std::string& logger_level);

#endif // OPTIONS_H