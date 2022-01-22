#include "interfaces/AMPL/AMPLModel.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategyFactory.hpp"
#include "Uno.hpp"
#include "optimization/ScaledReformulation.hpp"
#include "optimization/SlackReformulation.hpp"
#include "optimization/ElasticReformulation.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "linear_algebra/CSCSymmetricMatrix.hpp"

// new() overload to track heap allocations
size_t total_allocations = 0;

/*
void* operator new(size_t size) {
   std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}
*/

void run_uno_ampl(const std::string& problem_name, const Options& options) {
   // TODO: use a factory
   // AMPL model
   auto original_problem = std::make_unique<AMPLModel>(problem_name);
   INFO << "Heap allocations after AMPL: " << total_allocations << "\n";

   auto reformulated_problem = std::make_unique<SlackReformulation>(*original_problem);
   //const Problem& reformulated_problem = (options.at("subproblem") == "barrier") ? SlackReformulation(*original_problem) : *original_problem;

   // initial primal and dual points
   Iterate first_iterate(reformulated_problem->number_variables, reformulated_problem->number_constraints);
   reformulated_problem->get_initial_primal_point(first_iterate.x);
   reformulated_problem->get_initial_dual_point(first_iterate.multipliers.constraints);
   // project x into the bounds
   reformulated_problem->project_point_in_bounds(first_iterate.x);

   // initialize the function scaling TODO put this constant in option file
   Scaling scaling(reformulated_problem->number_constraints, 100.);
   // function scaling
   const bool scale_functions = (options.at("scale_functions") == "yes");
   if (scale_functions) {
      // evaluate the gradients at the current point
      first_iterate.evaluate_objective_gradient(*reformulated_problem);
      first_iterate.evaluate_constraint_jacobian(*reformulated_problem);
      scaling.compute(first_iterate.objective_gradient, first_iterate.constraint_jacobian);
      // forget about these evaluations
      first_iterate.reset_evaluations();
   }
   const Problem& problem_to_solve = ScaledReformulation(*reformulated_problem, scaling);

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(problem_to_solve, options);
   INFO << "Heap allocations after ConstraintRelax, Subproblem and Solver: " << total_allocations << "\n";

   // create the globalization mechanism
   auto mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   INFO << "Heap allocations after Mechanism: " << total_allocations << "\n";

   Uno uno = Uno(*mechanism, options);

   INFO << "Heap allocations before solving: " << total_allocations << "\n";
   const bool enforce_linear_constraints = (options.at("enforce_linear_constraints") == "yes");
   Result result = uno.solve(problem_to_solve, first_iterate, enforce_linear_constraints);
   Uno::postsolve_solution(*reformulated_problem, scaling, result.solution, result.status);

   const bool print_solution = (options.at("print_solution") == "yes");
   result.print(print_solution);
   INFO << "Heap allocations: " << total_allocations << "\n";
}

Level Logger::logger_level = INFO;

void test_problem_with_slacks(const std::string& problem_name) {
   auto original_problem = std::make_unique<AMPLModel>(problem_name);
   // add slacks
   const Problem& problem_to_solve = SlackReformulation(*original_problem);

   // create the first iterate (slacks set to 0)
   Iterate first_iterate(problem_to_solve.number_variables, problem_to_solve.number_constraints);
   problem_to_solve.get_initial_primal_point(first_iterate.x);
   problem_to_solve.get_initial_dual_point(first_iterate.multipliers.constraints);

   // set slack values
   for (size_t i = original_problem->number_variables; i < problem_to_solve.number_variables; i++) {
      const Range bounds = {problem_to_solve.get_variable_lower_bound(i), problem_to_solve.get_variable_upper_bound(i)};
      first_iterate.x[i] = Subproblem::push_variable_to_interior(0., bounds);
   }

   std::cout << "x = "; print_vector(std::cout, first_iterate.x);
   std::cout << "multipliers = "; print_vector(std::cout, first_iterate.multipliers.constraints);

   for (size_t i = 0; i < problem_to_solve.number_variables; i++) {
      std::cout << "Bounds of x" << i << ": [" << problem_to_solve.get_variable_lower_bound(i) << ", " <<
         problem_to_solve.get_variable_upper_bound(i) << "]\n";
   }
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      std::cout << "Bounds of c" << j << ": [" << problem_to_solve.get_constraint_lower_bound(j) << ", " <<
                problem_to_solve.get_constraint_upper_bound(j) << "]\n";
   }

   // evaluations
   std::vector<SparseVector<double>> jacobian(problem_to_solve.number_constraints);
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      jacobian[j].reserve(problem_to_solve.number_variables);
   }
   problem_to_solve.evaluate_constraint_jacobian(first_iterate.x, jacobian);
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      std::cout << jacobian[j];
   }

   CSCSymmetricMatrix hessian(problem_to_solve.number_variables, problem_to_solve.get_hessian_maximum_number_nonzeros());
   problem_to_solve.evaluate_lagrangian_hessian(first_iterate.x, 1., first_iterate.multipliers.constraints, hessian);
   std::cout << hessian;
}

void test_problem_with_elastics(const std::string& problem_name) {
   auto original_problem = std::make_unique<AMPLModel>(problem_name);
   // add slacks
   const Problem& problem_to_solve = ElasticReformulation(*original_problem, 1.);

   // create the first iterate (slacks set to 0)
   Iterate first_iterate(problem_to_solve.number_variables, problem_to_solve.number_constraints);
   problem_to_solve.get_initial_primal_point(first_iterate.x);
   problem_to_solve.get_initial_dual_point(first_iterate.multipliers.constraints);

   std::cout << "x = "; print_vector(std::cout, first_iterate.x);
   std::cout << "multipliers = "; print_vector(std::cout, first_iterate.multipliers.constraints);

   std::cout << "Variables\n";
   for (size_t i = 0; i < problem_to_solve.number_variables; i++) {
      std::cout << "Bounds of x" << i << ": [" << problem_to_solve.get_variable_lower_bound(i) << ", " <<
                problem_to_solve.get_variable_upper_bound(i) << "]\n";
   }
   std::cout << "Constraints\n";
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      std::cout << "Bounds of c" << j << ": [" << problem_to_solve.get_constraint_lower_bound(j) << ", " <<
                problem_to_solve.get_constraint_upper_bound(j) << "]\n";
   }

   // evaluations
   std::vector<SparseVector<double>> jacobian(problem_to_solve.number_constraints);
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      jacobian[j].reserve(problem_to_solve.number_variables);
   }
   problem_to_solve.evaluate_constraint_jacobian(first_iterate.x, jacobian);
   for (size_t j = 0; j < problem_to_solve.number_constraints; j++) {
      std::cout << jacobian[j];
   }

   CSCSymmetricMatrix hessian(problem_to_solve.number_variables, problem_to_solve.get_hessian_maximum_number_nonzeros());
   problem_to_solve.evaluate_lagrangian_hessian(first_iterate.x, 1., first_iterate.multipliers.constraints, hessian);
   std::cout << hessian;
}

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
         std::cout << "To solve an AMPL problem, type ./uno_ampl path_to_file/file.nl\n";
         std::cout << "To choose a globalization mechanism, use the argument -mechanism [LS|TR]\n";
         std::cout << "To choose a globalization strategy, use the argument -strategy [penalty|filter|nonmonotone-filter]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint-relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a subproblem, use the argument -subproblem [QP|LP|barrier]\n";
         std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|byrd]\n";
         std::cout << "The options can be combined in the same command line. Autocompletion is active.\n";
      }
      else {
         // run Uno on the .nl file (last command line argument)
         std::string problem_name = std::string(argv[argc - 1]);
         run_uno_ampl(problem_name, options);
      }
   }
   return EXIT_SUCCESS;
}
