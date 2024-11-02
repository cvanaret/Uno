#include <iostream>
#include <string>
#include <stdexcept>
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "hs71.hpp"
#include "Uno.hpp"
#include "model/ModelFactory.hpp"
#include "options/Options.hpp"
#include "options/DefaultOptions.hpp"
#include "tools/Logger.hpp"



int main() {

   uno::Options options = uno::DefaultOptions::load();
   uno::Options solvers_options = uno::DefaultOptions::determine_solvers_and_preset();
   // uno::Options::set_preset(solvers_options, "ipopt");
   options.overwrite_with(solvers_options);

   std::unique_ptr<uno::Model> hs_model = std::make_unique<local::HS71>();
   
   uno::Iterate initial_iterate(hs_model->number_variables, hs_model->number_constraints);
   hs_model->initial_primal_point(initial_iterate.primals);
   hs_model->project_onto_variable_bounds(initial_iterate.primals);
   hs_model->initial_dual_point(initial_iterate.multipliers.constraints);
   initial_iterate.feasibility_multipliers.reset();

   std::unique_ptr<uno::Model> model = uno::ModelFactory::reformulate(std::move(hs_model), initial_iterate, options);

   auto constraint_relaxation_strategy = uno::ConstraintRelaxationStrategyFactory::create(*model, options);
   auto globalization_mechanism = uno::GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   uno::Uno uno = uno::Uno(*globalization_mechanism, options);

   // solve the instance
   uno::Result result = uno.solve(*model, initial_iterate, options);
   uno.print_optimization_summary(result);

   return EXIT_SUCCESS;
}