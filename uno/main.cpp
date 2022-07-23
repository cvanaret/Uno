// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "Uno.hpp"
#include "optimization/ModelFactory.hpp"
#include "optimization/ScaledModel.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Timer.hpp"

void run_uno_ampl(const std::string& model_name, const Options& options) {
   // AMPL model
   auto original_model = ModelFactory::create(model_name, options);

   // initial primal and dual points
   Iterate first_iterate(original_model->number_variables, original_model->number_constraints);
   original_model->get_initial_primal_point(first_iterate.primals);
   original_model->get_initial_dual_point(first_iterate.multipliers.constraints);
   // project primal duals into the bounds
   original_model->project_point_onto_bounds(first_iterate.primals);

   // initialize the function scaling
   Scaling scaling(original_model->number_constraints, stod(options.at("function_scaling_threshold")));
   // function scaling
   const bool scale_functions = (options.at("scale_functions") == "yes");
   if (scale_functions) {
      // evaluate the gradients at the current point
      first_iterate.evaluate_objective_gradient(*original_model);
      first_iterate.evaluate_constraint_jacobian(*original_model);
      scaling.compute(first_iterate.original_evaluations.objective_gradient, first_iterate.original_evaluations.constraint_jacobian);
      // forget about these evaluations
      first_iterate.reset_evaluations();
   }
   const Model& model_to_solve = ScaledModel(*original_model, scaling);

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(model_to_solve, options);

   // create the globalization mechanism
   auto mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);

   // instantiate the combination of ingredients and solve the problem
   Uno uno = Uno(*mechanism, options);
   Result result = uno.solve(model_to_solve, first_iterate, options);
   Uno::postsolve_solution(*original_model, scaling, result.solution, result.status);

   std::string combination = options.at("mechanism") + " " + options.at("constraint-relaxation") + " " + options.at("strategy") + " " + options.at("subproblem");
   std::cout << "\nUno (" << combination << "): optimization summary\n";
   std::cout << Timer::get_current_date();
   std::cout << "==============================\n";

   const bool print_solution = (options.at("print_solution") == "yes");
   result.print(print_solution);
}

Level Logger::logger_level = INFO;

int main(int argc, char* argv[]) {
   if (1 < argc) {
      // get the default options
      Options options = get_default_options("uno.options");
      // get the command line options
      get_command_line_options(argc, argv, options);
      set_logger(options.at("logger"));

      options.print();

      if (std::string(argv[1]) == "-v") {
         std::cout << "Welcome in Uno\n";
         std::cout << "To solve an AMPL model, type ./uno_ampl path_to_file/file.nl\n";
         std::cout << "To choose a globalization mechanism, use the argument -mechanism [LS|TR]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint-relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a globalization strategy, use the argument -strategy [penalty|filter|nonmonotone-filter]\n";
         std::cout << "To choose a subproblem method, use the argument -subproblem [QP|LP|barrier]\n";
         std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|byrd]\n";
         std::cout << "The options can be combined in the same command line. Autocompletion is active.\n";
      }
      else {
         // run Uno on the .nl file (last command line argument)
         std::string model_name = std::string(argv[argc - 1]);
         run_uno_ampl(model_name, options);
      }
   }
   return EXIT_SUCCESS;
}
