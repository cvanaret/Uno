#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategyFactory.hpp"
#include "Uno.hpp"
#include "optimization/ModelFactory.hpp"
#include "optimization/ScaledModel.hpp"
#include "optimization/EqualityConstrainedModel.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Timer.hpp"

// new() overload to track heap allocations
size_t total_allocations = 0;

/*
void* operator new(size_t size) {
   std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}
 */

void run_uno_ampl(const std::string& model_name, const Options& options) {
   // AMPL model
   auto original_model = ModelFactory::create(model_name);
   INFO << "Heap allocations after AMPL: " << total_allocations << '\n';

   //auto reformulated_model = std::make_unique<EqualityConstrainedModel>(*original_model);
   Model* reformulated_model = original_model.get();

   // initial primal and dual points
   Iterate first_iterate(reformulated_model->number_variables, reformulated_model->number_constraints);
   reformulated_model->get_initial_primal_point(first_iterate.primals);
   reformulated_model->get_initial_dual_point(first_iterate.multipliers.constraints);
   // project x into the bounds
   reformulated_model->project_point_in_bounds(first_iterate.primals);

   // initialize the function scaling
   Scaling scaling(reformulated_model->number_constraints, stod(options.at("function_scaling_threshold")));
   // function scaling
   const bool scale_functions = (options.at("scale_functions") == "yes");
   if (scale_functions) {
      // evaluate the gradients at the current point
      first_iterate.evaluate_objective_gradient(*reformulated_model);
      first_iterate.evaluate_constraint_jacobian(*reformulated_model);
      scaling.compute(first_iterate.original_evaluations.objective_gradient, first_iterate.original_evaluations.constraint_jacobian);
      // forget about these evaluations
      first_iterate.reset_evaluations();
   }
   const Model& model_to_solve = ScaledModel(*reformulated_model, scaling);

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(model_to_solve, options);
   INFO << "Heap allocations after ConstraintRelaxation, Subproblem and Solver: " << total_allocations << '\n';

   // create the globalization mechanism
   auto mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   INFO << "Heap allocations after Mechanism: " << total_allocations << '\n';

   Uno uno = Uno(*mechanism, options);

   INFO << "Heap allocations before solving: " << total_allocations << '\n';
   const bool enforce_linear_constraints = (options.at("enforce_linear_constraints") == "yes");
   Result result = uno.solve(model_to_solve, first_iterate, enforce_linear_constraints);
   Uno::postsolve_solution(*reformulated_model, scaling, result.solution, result.status);

   std::string combination = options.at("mechanism") + " " + options.at("constraint-relaxation") + " " + options.at("strategy") + " " + options.at("subproblem");
   std::cout << "\nUno (" << combination << "): optimization summary\n";
   std::cout << Timer::get_current_date();
   std::cout << "==============================\n";

   const bool print_solution = (options.at("print_solution") == "yes");
   result.print(print_solution);
   INFO << "Heap allocations: " << total_allocations << '\n';
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
         std::cout << "To choose a globalization strategy, use the argument -strategy [penalty|filter|nonmonotone-filter]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint-relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a subproblem method, use the argument -subproblem [QP|LP|barrier]\n";
         std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|squid]\n";
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
