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
   void run_uno_ampl(const std::string& model_name, Options& options) {
      const AMPLModel model(model_name);
      Uno uno{};
      Result result = uno.solve(model, options);
      if (result.optimization_status == OptimizationStatus::SUCCESS) {
         // check result.solution_status
      }
      else {
         // ...
      }
      if (options.get_bool("AMPL_write_solution_to_file")) {
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

      // gather the options
      Options options;
      DefaultOptions::load(options);
      // set default preset
      Presets::set_default(options);

      // the -AMPL flag indicates that the solution should be written to the AMPL solution file
      size_t offset = 2;
      if (argc > 2 && std::string(argv[2]) == "-AMPL") {
         options.set_bool("AMPL_write_solution_to_file", true);
         ++offset;
      }
      else {
         options.set_bool("AMPL_write_solution_to_file", false);
      }
      // get the command line arguments (options start at index offset)
      const auto command_line_options = Options::get_command_line_options(argc, argv, offset);
      std::optional<std::string> optional_option_file{};
      std::optional<std::string> optional_preset{};
      for (const auto& [option_name, option_value]: command_line_options) {
         if (option_name == "option_file") {
            optional_option_file = option_value;
         }
         else if (option_name == "preset") {
            optional_preset = option_value;
         }
      }

      // [optional] set options from an option file
      if (optional_option_file.has_value()) {
         Options::load_option_file(options, *optional_option_file);
      }

      // [optional] set a preset
      if (optional_preset.has_value()) {
         Presets::set(options, *optional_preset);
      }

      // set the rest of the command line options (note: we add preset to the options for debugging purposes)
      for (const auto& [option_name, option_value]: command_line_options) {
         if (option_name != "option_file" && option_name != "logger") {
            options.set(option_name, option_value);
         }
      }

      // solve the model
      Logger::set_logger(options.get_string("logger"));
      run_uno_ampl(model_name, options);
   }
   return EXIT_SUCCESS;
}
