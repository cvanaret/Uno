#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cassert>
#include <vector>
#include <map>
#include "AMPLModel.hpp"
#include "SubproblemFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "GlobalizationMechanismFactory.hpp"
#include <ConstraintRelaxationStrategyFactory.hpp>
#include "FeasibilityRestoration.hpp"
#include "Uno.hpp"
#include "Logger.hpp"

// new overload to track heap allocations
size_t total_allocations = 0;

void* operator new(size_t size) {
   // std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}

void run_uno(const std::string& problem_name, const std::map<std::string, std::string>& options) {
   // TODO: use a factory
   // AMPL model (Hessian with a Fortran indexing starting at 1)
   auto problem = std::make_unique<AMPLModel>(problem_name, 1);
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

std::map<std::string, std::string> get_command_options(int argc, char* argv[], std::map<std::string, std::string>& options) {
   /* build the (argument, value) map */
   int i = 1;
   while (i < argc - 1) {
      std::string argument_i = std::string(argv[i]);

      if (argument_i[0] == '-' && i < argc - 1) {
         std::string value_i = std::string(argv[i + 1]);
         std::cout << "(" << argument_i << ", " << value_i << ")\n";
         options[argument_i.substr(1)] = value_i;
         i += 2;
      }
      else {
         std::cout << "Argument " << argument_i << " was ignored\n";
         i++;
      }
   }
   return options;
}

Level Logger::logger_level = DEBUG;

void set_logger(std::map<std::string, std::string> options) {
   try {
      const std::string logger_level = options.at("logger");

      if (logger_level == "ERROR") {
         Logger::logger_level = ERROR;
      }
      else if (logger_level == "WARNING") {
         Logger::logger_level = WARNING;
      }
      else if (logger_level == "INFO") {
         Logger::logger_level = INFO;
      }
      else if (logger_level == "DEBUG") {
         Logger::logger_level = DEBUG;
      }
   }
   catch (std::out_of_range&) {
   }
}

std::map<std::string, std::string> get_options(const std::string& file_name) {
   std::ifstream file;
   file.open(file_name);
   if (!file) {
      throw std::invalid_argument("The configuration file was not found");
   }
   else {
      /* register the default values */
      std::map<std::string, std::string> options;
      std::string key, value;
      std::string line;
      while (std::getline(file, line)) {
         if (!line.empty() && line.find("#") != 0) {
            std::istringstream iss;
            iss.str(line);
            iss >> key >> value;
            std::cout << "Option " << key << " = " << value << "\n";
            options[key] = value;
         }
      }
      file.close();
      return options;
   }
}

int main(int argc, char* argv[]) {
   if (1 < argc) {
      /* get the options */
      std::map<std::string, std::string> options = get_options("uno.cfg");
      /* get the command line options */
      get_command_options(argc, argv, options);
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
