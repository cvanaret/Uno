#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <functional>
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

/*
void* operator new(size_t size) {
   //std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}
 */

struct PredictedReduction {
   const double full_step_value;
   const std::function<double (double step_length)> partial_step_value;

   PredictedReduction(double full_step_value, const std::function<double ()>& partial_step_value);
   double evaluate(double step_length);
};

void run_uno(const std::string& problem_name, const Options& options) {
   const std::string_view mechanism_type = options.at("mechanism");
   const std::string_view constraint_relaxation_type = options.at("constraint-relaxation");
   const std::string_view subproblem_type = options.at("subproblem");

   // TODO: use a factory
   // AMPL model
   auto problem = std::make_unique<AMPLModel>(problem_name);
   INFO << "Heap allocations after AMPL: " << total_allocations << "\n";

   /* create the subproblem */
   const bool use_trust_region = (mechanism_type == "TR");
   const size_t number_variables = ConstraintRelaxationStrategyFactory::get_number_variables(constraint_relaxation_type, *problem);
   auto subproblem = SubproblemFactory::create(*problem, number_variables, subproblem_type, options, use_trust_region);

   /* create the constraint relaxation strategy */
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(constraint_relaxation_type, *problem, *subproblem, options);
   INFO << "Heap allocations after ConstraintRelax, Subproblem and Solver: " << total_allocations << "\n";

   /* create the globalization mechanism */
   auto mechanism = GlobalizationMechanismFactory::create(mechanism_type, *constraint_relaxation_strategy, options);
   INFO << "Heap allocations after Mechanism: " << total_allocations << "\n";

   const double tolerance = std::stod(options.at("tolerance"));
   const int max_iterations = std::stoi(options.at("max_iterations"));
   Uno uno = Uno(*mechanism, tolerance, max_iterations);

   /* initial primal and dual points */
   Iterate first_iterate(problem->number_variables, problem->number_constraints);
   problem->set_initial_primal_point(first_iterate.x);
   problem->set_initial_dual_point(first_iterate.multipliers.constraints);

   INFO << "Heap allocations before solving: " << total_allocations << "\n";
   bool use_preprocessing = (options.at("preprocessing") == "yes");
   Result result = uno.solve(*problem, first_iterate, use_preprocessing);

   /* remove auxiliary variables */
   result.solution.x.resize(problem->number_variables);
   result.solution.multipliers.lower_bounds.resize(problem->number_variables);
   result.solution.multipliers.upper_bounds.resize(problem->number_variables);
   bool print_solution = (options.at("print_solution") == "yes");
   result.display(print_solution);
   INFO << "Heap allocations: " << total_allocations << "\n";
}

Level Logger::logger_level = INFO;

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
