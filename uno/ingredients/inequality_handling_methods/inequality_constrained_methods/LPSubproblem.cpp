// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"

namespace uno {
   LPSubproblem::LPSubproblem(const Options& options):
         InequalityConstrainedMethod(),
         solver(LPSolverFactory::create(options)) {
   }

   LPSubproblem::~LPSubproblem() { }

   void LPSubproblem::initialize(const OptimizationProblem& first_reformulation, const HessianModel& hessian_model) {
      InequalityConstrainedMethod::initialize(first_reformulation, hessian_model);
      this->solver->initialize_memory(first_reformulation, hessian_model);
   }

   void LPSubproblem::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   }

   void LPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, HessianModel& /*hessian_model*/, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      this->solver->solve_LP(problem, current_iterate, this->initial_point, direction, trust_region_radius, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   double LPSubproblem::hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.;
   }

   std::string LPSubproblem::get_strategy_combination() const {
      return "LP method";
   }
} // namespace