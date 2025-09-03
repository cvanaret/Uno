// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnoSolver.hpp"
#include "PythonModel.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "tools/Logger.hpp"

namespace uno {
   UnoSolver::UnoSolver() {
      DefaultOptions::load(this->options);
   }

   void UnoSolver::optimize(const Model& model) {
      // TODO reformulate (scale, add slacks, relax the bounds, ...) if necessary
      DISCRETE << "Reformulated model " << model.name << '\n' << model.number_variables << " variables, " <<
         model.number_constraints << " constraints (" << model.get_equality_constraints().size() <<
         " equality, " << model.get_inequality_constraints().size() << " inequality)\n";

      Logger::set_logger(this->options.get_string("logger"));
      DISCRETE << "Original model " << model.name << '\n' << model.number_variables << " variables, " <<
         model.number_constraints << " constraints (" << model.get_equality_constraints().size() <<
         " equality, " << model.get_inequality_constraints().size() << " inequality)\n";

      // initialize initial primal and dual points
      Iterate initial_iterate(model.number_variables, model.number_constraints);
      model.initial_primal_point(initial_iterate.primals);
      model.project_onto_variable_bounds(initial_iterate.primals);
      model.initial_dual_point(initial_iterate.multipliers.constraints);

      // solve the instance
      Result result = this->uno_solver.solve(model, initial_iterate, this->options);
   }
} // namespace