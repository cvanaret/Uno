// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <stdexcept>
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "AMPLModel.hpp"
#include "Uno.hpp"
#include "model/ModelFactory.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"

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

         // initialize initial primal and dual points
         Iterate initial_iterate(ampl_model->number_variables, ampl_model->number_constraints);
         ampl_model->initial_primal_point(initial_iterate.primals);
         ampl_model->project_onto_variable_bounds(initial_iterate.primals);
         ampl_model->initial_dual_point(initial_iterate.multipliers.constraints);
         initial_iterate.feasibility_multipliers.reset();

         // reformulate (scale, add slacks, relax the bounds, ...) if necessary
         std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(ampl_model), initial_iterate, options);

         // create the constraint relaxation strategy, the globalization mechanism and the Uno solver
         auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(*model, options);
         auto globalization_mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
         Uno uno = Uno(*globalization_mechanism, options);

         // solve the instance
         Result result = uno.solve(*model, initial_iterate, options);
         Uno::print_optimization_summary(options, result);
         // std::cout << "memory_allocation_amount = " << memory_allocation_amount << '\n';
      }
      catch (std::exception& exception) {
         ERROR << exception.what() << '\n';
      }
   }

   // argv[0] is ./uno_ampl
   // argv[1] is the model name
   // argv[2] should be -AMPL
   // argv[i] for i = 3..argc-1 are options
   void get_command_line_options(Options& options, int argc, char* argv[]) {
      static const std::string delimiter = "=";

      // build the (name, value) map
      for (int i = 3; i < argc; i++) {
         const std::string argument = std::string(argv[i]);
         size_t position = argument.find_first_of(delimiter);
         if (position == std::string::npos) {
            throw std::runtime_error("The option " + argument + " does not contain the delimiter " + delimiter + ".");
         }
         const std::string key = argument.substr(0, position);
         const std::string value = argument.substr(position + 1);
         if (key == "preset") {
            options.find_preset(value);
         }
         else {
            options[key] = value;
         }
      }
      options.print(false);
   }

   void print_uno_instructions() {
      std::cout << "Welcome in Uno 1.1.0\n";
      std::cout << "To solve an AMPL model, type ./uno_ampl model.nl -AMPL [key=value ...]\n";
      std::cout << "To choose a constraint relaxation strategy, use the argument constraint_relaxation_strategy="
                   "[feasibility_restoration|l1_relaxation]\n";
      std::cout << "To choose a subproblem method, use the argument subproblem=[QP|LP|primal_dual_interior_point]\n";
      std::cout << "To choose a globalization mechanism, use the argument globalization_mechanism=[LS|TR]\n";
      std::cout << "To choose a globalization strategy, use the argument globalization_strategy="
                   "[l1_merit|fletcher_filter_method|waechter_filter_method]\n";
      std::cout << "To choose a preset, use the argument preset=[filtersqp|ipopt|byrd]\n";
      std::cout << "The options can be combined in the same command line.\n";
      //std::cout << "Autocompletion is possible (see README).\n";
   }
} // namespace

int main(int argc, char* argv[]) {
   using namespace uno;

   if (argc == 1) {
      print_uno_instructions();
   }
   else if (argc == 2) {
      if (std::string(argv[1]) == "--v") {
         print_uno_instructions();
      }
      else if (std::string(argv[1]) == "--strategies") {
         Uno::print_available_strategies();
      }
      else {
         throw std::runtime_error("The second command line argument should be -AMPL.");
      }
   }
   else if (argc >= 3) {
      // get the default options
      Options options = Options::get_default_options("uno.options");

      // AMPL expects: ./uno_ampl model.nl -AMPL [key=value, ...]
      // model name
      std::string model_name = std::string(argv[1]);
      std::cout << "Model " << model_name << '\n';

      // -AMPL
      if (std::string(argv[2]) != "-AMPL") {
         throw std::runtime_error("The second command line argument should be -AMPL.");
      }

      // override them with the command line arguments
      get_command_line_options(options, argc, argv);

      // solve the model
      Logger::set_logger(options.get_string("logger"));
      run_uno_ampl(model_name, options);
   }
   return EXIT_SUCCESS;
}
