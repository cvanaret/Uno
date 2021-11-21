#ifndef OPTIONS_H
#define OPTIONS_H

#include <map>

using Options = std::map<std::string, std::string>;

Options get_default_options(const std::string& file_name);
void get_command_line_options(int argc, char* argv[], Options& options);
void print_options(const Options& options);
void set_logger(const Options& options);

#endif // OPTIONS_H