// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "AMPLModel.hpp"
#include "model/ModelFactory.hpp"
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
         // AMPL model
         std::unique_ptr<Model> ampl_model = std::make_unique<AMPLModel>(model_name, options);
         DISCRETE << "Original model " << ampl_model->name << '\n' << ampl_model->number_variables << " variables, " <<
            ampl_model->number_constraints << " constraints (" << ampl_model->get_equality_constraints().size() <<
            " equality, " << ampl_model->get_inequality_constraints().size() << " inequality)\n";

         // reformulate (scale, add slacks, relax the bounds, ...) if necessary
         std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(ampl_model), options);
         DISCRETE << "Reformulated model " << model->name << '\n' << model->number_variables << " variables, " <<
            model->number_constraints << " constraints (" << model->get_equality_constraints().size() <<
            " equality, " << model->get_inequality_constraints().size() << " inequality)\n";

         // initialize initial primal and dual points
         Iterate initial_iterate(model->number_variables, model->number_constraints);
         model->initial_primal_point(initial_iterate.primals);
         model->initial_dual_point(initial_iterate.multipliers.constraints);

         // solve the instance
         Uno uno{};
         const Result result = uno.solve(*model, initial_iterate, options);
         if (result.optimization_status == OptimizationStatus::SUCCESS) {
            // check result.solution.status
         }
         else {
            // ...
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
         Uno::print_instructions();
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
         // determine the default solvers based on the available libraries
         const Options subproblem_solvers = DefaultOptions::determine_subproblem_solvers();
         options.set(subproblem_solvers);
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
         Options command_line_options = Options::get_command_line_options(argc, argv, offset);
         // possibly set options from an option file
         const auto optional_option_file = command_line_options.get_string_optional("option_file");
         if (optional_option_file.has_value()) {
            Options file_options = Options::load_option_file(*optional_option_file);
            options.overwrite_with(file_options);
         }
         // possibly set a preset
         const auto optional_preset = command_line_options.get_string_optional("preset");
         Options preset_options = Presets::get_preset_options(optional_preset);
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