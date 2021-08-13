#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cassert>
#include <vector>
#include <map>
#include "interfaces/AMPL/AMPLModel.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategyFactory.hpp"
#include "Uno.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"

// new overload to track heap allocations
size_t total_allocations = 0;

void* operator new(size_t size) {
   // std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}

void run_uno(const std::string& problem_name, const Options& options) {
   // TODO: use a factory
   // AMPL model
   auto problem = std::make_unique<AMPLModel>(problem_name);
   std::cout << "Heap allocations after AMPL: " << total_allocations << "\n";

   /* create the constraint relaxation strategy and the subproblem */
   bool use_trust_region = (options.at("mechanism") == "TR");
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(options.at("constraint-relaxation"), *problem, options,
         use_trust_region);
   std::cout << "Heap allocations after ConstraintRelax, Subproblem and Solver: " << total_allocations << "\n";

   /* create the globalization mechanism */
   auto mechanism = GlobalizationMechanismFactory::create(options.at("mechanism"), *constraint_relaxation_strategy, options);
   std::cout << "Heap allocations after Mechanism: " << total_allocations << "\n";

   double tolerance = std::stod(options.at("tolerance"));
   int max_iterations = std::stoi(options.at("max_iterations"));
   bool use_preprocessing = (options.at("preprocessing") == "yes");
   Uno uno = Uno(*mechanism, tolerance, max_iterations);

   /* initial primal and dual points */
   std::vector<double> x(problem->number_variables);
   Multipliers multipliers(problem->number_variables, problem->number_constraints);
   problem->set_initial_primal_point(x);
   problem->set_initial_dual_point(multipliers.constraints);

   std::cout << "Heap allocations before solving: " << total_allocations << "\n";
   Result result = uno.solve(*problem, x, multipliers, use_preprocessing);

   /* remove auxiliary variables */
   result.solution.x.resize(problem->number_variables);
   result.solution.multipliers.lower_bounds.resize(problem->number_variables);
   result.solution.multipliers.upper_bounds.resize(problem->number_variables);
   bool print_solution = (options.at("print_solution") == "yes");
   result.display(print_solution);
   std::cout << "Heap allocations: " << total_allocations << "\n";
}

Level Logger::logger_level = DEBUG;

int main(int argc, char* argv[]) {
   if (1 < argc) {
      /* get the options */
      Options options = get_options("uno.cfg");
      /* get the command line options */
      get_command_options(argc, argv, options);
      print_options(options);
      set_logger(options);

      if (std::string(argv[1]) == "-v") {
         std::cout << "Welcome in UNO\n";
         std::cout << "To solve an AMPL problem, type ./uno_ampl path_to_file/file.nl\n";
         std::cout << "To choose a globalization mechanism, use the argument -mechanism [LS|TR]\n";
         std::cout << "To choose a globalization strategy, use the argument -strategy [penalty|filter|nonmonotone-filter]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint-relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a subproblem, use the argument -subproblem [SQP|SLP|IPM]\n";
         std::cout << "The four options can be combined in the same command line. Autocompletion is active.\n";
      }
      else {
         std::string problem_name = std::string(argv[argc - 1]);
         /* run Uno */
         run_uno(problem_name, options);
      }
   }
   return EXIT_SUCCESS;
}
