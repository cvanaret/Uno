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
      // the -AMPL flag indicates that the solution should be written to the AMPL solution file
      const bool write_solution_to_file = (argc > 2 && std::string(argv[2]) == "-AMPL");
      const size_t offset = write_solution_to_file ? 3 : 2;
      try {
         // get the command line arguments (options start at index offset)
         const auto command_line_options = Options::get_command_line_options(argc, argv, offset);

         // user options: option file first, then command line (takes precedence)
         Options user_options;
         for (const auto& [option_name, option_value]: command_line_options) {
            if (option_name == "option_file") {
               Options::load_option_file(user_options, option_value);
            }
         }
         for (const auto& [option_name, option_value]: command_line_options) {
            if (option_name != "option_file") {
               user_options.set(option_name, option_value);
            }
         }
         user_options.set_bool("write_solution_to_file", write_solution_to_file);

         // defaults -> preset -> user options
         Options options;
         DefaultOptions::load(options);
         Logger::set_logger(user_options.get_string_optional("logger").value_or(options.get_string("logger")));

         user_options.print("User options");

         // create the model
         const char* model_name = argv[1];
         const AMPLModel model(model_name);

         const std::string preset = user_options.get_string_optional("preset").value_or(options.get_string("preset"));
         Presets::set(model, options, preset);
         options.overwrite(user_options);

         // solve the model
         run_uno_ampl(model, options);
      }
      catch (const std::exception& e) {
         std::cout << "uno_ampl failed with the following error: " << e.what() << '\n';
         return EXIT_FAILURE;
      }
   }
   return EXIT_SUCCESS;
}
