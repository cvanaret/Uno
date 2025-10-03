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
   void run_uno_ampl(const std::string& model_name, const Options& options) {
      try {
         const AMPLModel model(model_name);
         Uno uno{};
         Result result = uno.solve(model, options);
         if (result.optimization_status == OptimizationStatus::SUCCESS) {
            // check result.solution.status
         }
         else {
            // ...
         }
         if (options.get_bool("AMPL_write_solution_to_file")) {
            model.write_solution_to_file(result);
         }
         // std::cout << "memory_allocation_amount = " << memory_allocation_amount << '\n';
      }
      catch (std::exception& exception) {
         DISCRETE << exception.what() << '\n';
      }
   }
} // namespace

int main(int argc, char* argv[]) {
   using namespace uno;

   try {
      if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--v")) {
         std::cout << "Uno " << Uno::current_version() << '\n';
      }
      else if (argc == 2 && std::string(argv[1]) == "--strategies") {
         Uno::print_available_strategies();
      }
      else {
         // AMPL expects: ./uno_ampl model.nl [-AMPL] [option_name=option_value, ...]
         // model name
         std::string model_name = std::string(argv[1]);

         // gather the options
         Options options;
         DefaultOptions::load(options);
         // the -AMPL flag indicates that the solution should be written to the AMPL solution file
         size_t offset = 2;
         if (argc > 2 && std::string(argv[2]) == "-AMPL") {
            options.set("AMPL_write_solution_to_file", "yes");
            ++offset;
         }
         else {
            options.set("AMPL_write_solution_to_file", "no");
         }
         // get the command line arguments (options start at index offset)
         const Options command_line_options = Options::get_command_line_options(argc, argv, offset);
         // possibly set options from an option file
         const auto optional_option_file = command_line_options.get_string_optional("option_file");
         if (optional_option_file.has_value()) {
            Options file_options = Options::load_option_file(*optional_option_file);
            options.overwrite_with(file_options);
         }
         // possibly set a preset
         const auto optional_preset = command_line_options.get_string_optional("preset");
         const Options preset_options = Presets::get_preset_options(optional_preset);
         options.overwrite_with(preset_options);
         // overwrite the options with the command line arguments
         options.overwrite_with(command_line_options);

         // solve the model
         Logger::set_logger(options.get_string("logger"));
         run_uno_ampl(model_name, options);
      }
   }
   catch (std::exception& exception) {
      DISCRETE << exception.what() << '\n';
   }
   return EXIT_SUCCESS;
}