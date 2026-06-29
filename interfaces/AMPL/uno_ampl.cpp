// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "AMPLModel.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "Uno.hpp"

/*
size_t memory_allocation_amount = 0;

void* operator new(size_t size) {
   memory_allocation_amount += size;
   std::cout << "Memory: " << size << '\n';
   return malloc(size);
}
*/

namespace uno {
   void run_uno_ampl(const AMPLModel& model, Options& options) {
      Uno uno{};
      Result result = uno.solve(model, options);
      if (options.get_bool("write_solution_to_file")) {
         model.write_solution_to_file(result);
      }
      // std::cout << "memory_allocation_amount = " << memory_allocation_amount << '\n';
   }
} // namespace

int main(int argc, char* argv[]) {
   using namespace uno;

   if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--v")) {
      std::cout << "Uno " << Uno::current_version() << '\n';
   }
   else if (argc == 2 && std::string(argv[1]) == "--dump-options") {
      // Print all available options (type + default value) for automated tools
      Options::dump_default_options();
   }
   else if (argc == 2 && std::string(argv[1]) == "--strategies") {
      Uno::print_available_strategies();
   }
   else { // argc >= 2
      // AMPL expects: ./uno_ampl model.nl [-AMPL] [option_name=option_value, ...]
      // model name
      const char* model_name = argv[1];
      const AMPLModel model(model_name);

      // set the default options
      Options options;
      DefaultOptions::load(options);

      // the -AMPL flag indicates that the solution should be written to the AMPL solution file
      size_t offset = 2;
      if (argc > 2 && std::string(argv[2]) == "-AMPL") {
         options.set_bool("write_solution_to_file", true);
         ++offset;
      }
      else {
         options.set_bool("write_solution_to_file", false);
      }
      try {
         // get the command line arguments (options start at index offset)
         const auto command_line_options = Options::get_command_line_options(argc, argv, offset);

         // [optional] read an option file
         std::optional<std::string> optional_option_file{};
         std::string preset = "auto";
         for (const auto& [option_name, option_value]: command_line_options) {
            if (option_name == "option_file") {
               optional_option_file = option_value;
            }
            else if (option_name == "preset") {
               preset = option_value;
            }
         }
         if (optional_option_file.has_value()) {
            Options::load_option_file(options, *optional_option_file);
         }

         // set the preset (default: auto)
         const std::optional<std::string> optional_preset = options.get_string_optional("preset");
         if (optional_preset.has_value()) {
            preset = *optional_preset;
         }
         Presets::set(model, options, preset);

         // set the rest of the command line options
         for (const auto& [option_name, option_value]: command_line_options) {
            if (option_name != "option_file") {
               options.set(option_name, option_value);
            }
         }

         // solve the model
         Logger::set_logger(options.get_string("logger"));
         run_uno_ampl(model, options);
      }
      catch (const std::exception& e) {
         std::cout << "uno_ampl failed with the following error: " << e.what() << '\n';
         return EXIT_FAILURE;
      }
   }
   return EXIT_SUCCESS;
}
