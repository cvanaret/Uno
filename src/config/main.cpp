#include <iostream>
#include "ConfigurationFile.hpp"

int main(int argc, char** argv) {
	try {
		const char* file_name = argv[1];

		/* read the configuration file */
		ConfigurationFile config_file(file_name);

		int width = config_file.get("Screen Width");
		int height = config_file.get("Screen Height");
		std::string monitorType = config_file.get("Monitor Type");
		bool dpmsEnabled = config_file.get("DPMS Enabled");

		std::cout << "Screen Width : " << width << "\n"
			<< "Screen Height: " << height << "\n"
			<< "Monitor Type : " << monitorType << "\n"
			<< "DPMS Enabled : " << (dpmsEnabled ? "true" : "false")
			<< "\n";
	}
	catch (std::exception& e) {
		std::cerr << e.what() << "\n";
		return -1;
	}
}
