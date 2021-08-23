#ifndef OPTIONS_H
#define OPTIONS_H

#include <unordered_map>

using Options = std::unordered_map<std::string, std::string>;

Options get_options(const std::string& file_name);
void get_command_options(int argc, char* argv[], Options& options);
void print_options(const Options& options);
void set_logger(const Options& options);

#endif // OPTIONS_H