// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "model/CustomModel.hpp"
#include "Uno.hpp"
#include "model/ModelFactory.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"

double objective(const Vector<double>& x) {
   return 100.*std::pow(x[1] - std::pow(x[0], 2.), 2.) + std::pow(1. - x[0], 2.);
}

void objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) {
   gradient.insert(0, -400.*x[0]*(x[1] - std::pow(x[0], 2.)) - 2.*(1. - x[0]));
   gradient.insert(0, 200.*(x[1] - std::pow(x[0], 2.)));
}

void constraints(const Vector<double>& x, std::vector<double>& constraints) {
   constraints[0] = x[0]*x[1];
   constraints[1] = x[0] + std::pow(x[1], 2.);
}

void constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) {
   // first constraint
   constraint_jacobian[0].insert(0, x[1]);
   constraint_jacobian[0].insert(1, x[0]);
   // second constraint
   constraint_jacobian[1].insert(0, 1.);
   constraint_jacobian[1].insert(1, 2.*x[1]);
}

void lagrangian_hessian(const Vector<double>& /*x*/, double /*objective_multiplier*/, const Vector<double>& /*multipliers*/,
      SymmetricMatrix<double>& /*hessian*/) {
   // TODO
}

void run_uno(const Options& options) {
   const size_t number_variables = 2;
   const size_t number_constraints = 2;
   const double objective_sign = 1.;
   const CustomModel custom_model = CustomModel{"custom", number_variables, number_constraints, objective_sign, objective, objective_gradient,
         constraints, constraint_jacobian, lagrangian_hessian};

   // initialize initial primal and dual points
   Iterate initial_iterate(custom_model.number_variables, custom_model.number_constraints);
   custom_model.initial_primal_point(initial_iterate.primals);
   custom_model.project_onto_variable_bounds(initial_iterate.primals);
   custom_model.initial_dual_point(initial_iterate.multipliers.constraints);
   initial_iterate.feasibility_multipliers.reset();

   // reformulate (scale, add slacks, relax the bounds, ...) if necessary
   //std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(custom_model), initial_iterate, options);
   
   // create the constraint relaxation strategy, the globalization mechanism and the Uno solver
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(custom_model, options);
   auto globalization_mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   Uno uno = Uno(*globalization_mechanism, options);

   // solve the instance
   Result result = uno.solve(custom_model, initial_iterate, options);
   Uno::print_optimization_summary(options, result);
}

Level Logger::level = INFO;

int main(int argc, char* argv[]) {
   // get the default options
   Options options = get_default_options("uno.options");
   // override them with the command line arguments
   get_command_line_arguments(argc, argv, options);
   Logger::set_logger(options.get_string("logger"));

   if (1 < argc && std::string(argv[1]) == "-v") {
      Uno::print_uno_version();
   }
   else if (1 < argc && std::string(argv[1]) == "--strategies") {
      Uno::print_available_strategies();
   }
   else {
      run_uno(options);
   }
   return EXIT_SUCCESS;
}
