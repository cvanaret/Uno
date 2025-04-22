// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "options/Options.hpp"

namespace uno {
   QPSubproblem::QPSubproblem(size_t number_variables, size_t number_constraints, size_t number_hessian_nonzeros, const Options& options) :
         InequalityConstrainedMethod(number_variables, number_constraints),
         enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
         // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
         solver(QPSolverFactory::create(number_variables,
               // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
               std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, number_hessian_nonzeros),
               options)) {
   }

   QPSubproblem::~QPSubproblem() { }

   void QPSubproblem::initialize(const OptimizationProblem& first_reformulation) {
      this->solver->initialize_memory(first_reformulation);
   }

   void QPSubproblem::generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) {
      if (this->enforce_linear_constraints_at_initial_iterate) {
         Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, *this->solver);
      }
   }

   void QPSubproblem::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         Direction& direction, HessianModel& hessian_model, WarmstartInformation& warmstart_information) {
      this->solver->solve_QP(statistics, problem, current_iterate, current_multipliers, this->initial_point, direction,
         hessian_model, this->trust_region_radius, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   double QPSubproblem::hessian_quadratic_product(const Vector<double>& primal_direction) const {
      return this->solver->hessian_quadratic_product(primal_direction);
   }
} // namespace
