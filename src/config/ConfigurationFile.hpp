#ifndef CONFIGURATIONFILE_H
#define CONFIGURATIONFILE_H

#include <map>

class Value {
	public:
		Value(const std::string& value = "");

		operator int() const;

		operator const char*() const;

		operator std::string() const;

		operator bool() const;

	private:
		std::string value;
};


class ConfigurationFile {
	public:
		ConfigurationFile(const char* file_name, const char field_separator = ':');

		Value get(const std::string& key) const;

		friend std::ostream& operator<<(std::ostream& os, const ConfigurationFile& config);

	private:
		std::string trim(std::string str);

		std::map<std::string, Value> config_list;
};

#endif // CONFIGURATIONFILE_H
