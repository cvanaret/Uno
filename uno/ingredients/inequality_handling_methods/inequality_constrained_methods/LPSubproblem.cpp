// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSubproblem.hpp"
#include "ingredients/local_models/LagrangeNewtonSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/LPSolver.hpp"
#include "solvers/LPSolverFactory.hpp"

namespace uno {
   LPSubproblem::LPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, const Options& options) :
         InequalityConstrainedMethod("zero", number_variables, 0, false, options),
         solver(LPSolverFactory::create(number_variables, number_constraints,
               number_objective_gradient_nonzeros, number_jacobian_nonzeros, options)) {
   }

   LPSubproblem::~LPSubproblem() { }

   void LPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   }

   void LPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, const WarmstartInformation& warmstart_information) {
      const LagrangeNewtonSubproblem subproblem(problem, current_iterate, current_multipliers.constraints, *this->hessian_model, this->trust_region_radius);
      this->solver->solve_LP(subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }
} // namespace