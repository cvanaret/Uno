// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <stdexcept>
#include "AMPLModel.hpp"
#include "AMPLUserCallbacks.hpp"
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
         model->project_onto_variable_bounds(initial_iterate.primals);
         model->initial_dual_point(initial_iterate.multipliers.constraints);
         initial_iterate.feasibility_multipliers.reset();

         // solve the instance
         Uno uno{model->number_constraints, options};
         AMPLUserCallbacks user_callbacks{};
         const Result result = uno.solve(*model, initial_iterate, options, user_callbacks);
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

   void print_uno_instructions() {
      std::cout << "Welcome in Uno " << Uno::current_version() << '\n';
      std::cout << "To solve an AMPL model, type ./uno_ampl model.nl -AMPL [option_name=option_value ...]\n";
      std::cout << "To choose a constraint relaxation strategy, use the argument constraint_relaxation_strategy="
                   "[feasibility_restoration|l1_relaxation]\n";
      std::cout << "To choose a subproblem method, use the argument subproblem=[QP|LP|primal_dual_interior_point]\n";
      std::cout << "To choose a globalization mechanism, use the argument globalization_mechanism=[LS|TR]\n";
      std::cout << "To choose a globalization strategy, use the argument globalization_strategy="
                   "[l1_merit|fletcher_filter_method|waechter_filter_method]\n";
      std::cout << "To choose a preset, use the argument preset=[filtersqp|ipopt|byrd]\n";
      std::cout << "The options can be combined in the same command line.\n";
   }
} // namespace

int main(int argc, char* argv[]) {
   using namespace uno;

   try {
      if (argc == 1) {
         print_uno_instructions();
      }
      else if (argc == 2 && std::string(argv[1]) == "--v") {
         print_uno_instructions();
      }
      else if (argc == 2 && std::string(argv[1]) == "--strategies") {
         Uno::print_available_strategies();
      }
      else {
         // AMPL expects: ./uno_ampl model.nl [-AMPL] [option_name=option_value, ...]
         // model name
         std::string model_name = std::string(argv[1]);

         // gather the options
         Options options = DefaultOptions::load();
         // determine the default solvers based on the available libraries
         Options solvers_options = DefaultOptions::determine_solvers();
         options.overwrite_with(solvers_options);
         // the -AMPL flag indicates that the solution should be written to the AMPL solution file
         size_t offset = 2;
         if (std::string(argv[2]) == "-AMPL") {
            options["AMPL_write_solution_to_file"] = "yes";
            ++offset;
         }
         else {
            options["AMPL_write_solution_to_file"] = "no";
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