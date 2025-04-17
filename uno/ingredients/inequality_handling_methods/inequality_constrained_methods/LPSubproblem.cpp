// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "options/Options.hpp"

namespace uno {
   LPSubproblem::LPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, const Options& options) :
         InequalityConstrainedMethod("zero", number_variables, number_constraints, 0, false, options),
         solver(LPSolverFactory::create(number_variables, number_constraints,
               number_objective_gradient_nonzeros, number_jacobian_nonzeros, options)) {
   }

   LPSubproblem::~LPSubproblem() { }

   void LPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   }

   void LPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, WarmstartInformation& warmstart_information) {
      this->solver->solve_LP(problem, current_iterate, this->initial_point, direction, this->trust_region_radius, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   double LPSubproblem::hessian_quadratic_product(const Vector<double>& /*primal_direction*/) const {
      return 0.;
   }
} // namespace