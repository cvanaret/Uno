#include <fstream>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include "ConfigurationFile.hpp"

Value::Value(const std::string& value) : value(value) {
}

Value::operator int() const {
    std::stringstream stream(value);
    int int_value = 0;
    stream >> int_value;

    if (stream.fail()) {
        throw std::runtime_error("Cannot convert to int.");
    }
    return int_value;
}

Value::operator const char*() const {
    return this->value.c_str();
}

Value::operator std::string() const {
    return this->value;
}

Value::operator bool() const {
    return this->value == "true" || this->value == "True" || this->value == "TRUE";
}

ConfigurationFile::ConfigurationFile(const char* file_name, const char field_separator) {
    std::fstream file(file_name, std::ios::in);

    if (file) {
        std::string line;

        while (getline(file, line)) {
            size_t sep = line.find_first_of(field_separator);

            if (sep != std::string::npos) {
                std::string key = trim(line.substr(0, sep));
                std::string value = trim(line.substr(sep + 1));

                if (!key.empty() && !value.empty()) {
                    config_list[key] = Value(value);
                }
                else {
                    throw std::runtime_error("Error within configuration file at line: " + line);
                }
            }
            else {
                throw std::runtime_error("Error within configuration file at line: " + line);
            }
        }
    }
    else {
        throw std::runtime_error("Cannot open config file.");
    }
}

Value ConfigurationFile::get(const std::string& key) const {

    try {
        return this->config_list.at(key);
    }
    catch (std::out_of_range) {
        throw std::runtime_error("Cannot find config item.");
    }
}

std::ostream& operator<<(std::ostream& os, const ConfigurationFile& config) {
    for (auto element : config.config_list) {
        const std::string& key = element.first;
        const std::string val = element.second;

        std::cout << key << "\t" << val << "\n";
    }

    return os;
}

std::string ConfigurationFile::trim(std::string str) {
    size_t pos = str.find_first_not_of(" \t\n");

    if (pos != std::string::npos) {
        str.erase(0, pos);
    }

    pos = str.find_last_not_of(" \t\n");

    if (pos != std::string::npos) {
        str.erase(pos + 1);
    }

    return str;
}
