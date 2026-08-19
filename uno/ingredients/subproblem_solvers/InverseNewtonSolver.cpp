// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InverseNewtonSolver.hpp"
#include "ingredients/hessian_models/quasi_newton/inverse/InverseLBFGSHessian.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   InverseNewtonSolver::InverseNewtonSolver(InverseLBFGSHessian& hessian_model): hessian_model(hessian_model) {
   }

   void InverseNewtonSolver::initialize_memory(const Subproblem& subproblem) {
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->rhs.resize(subproblem.number_variables);
   }

   void InverseNewtonSolver::generate_initial_iterate(const Subproblem& /*subproblem*/, Iterate& /*initial_iterate*/,
         Evaluations& /*evaluations*/) {
      // do nothing
   }

   void InverseNewtonSolver::compute_least_squares_multipliers(const Subproblem& /*subproblem*/, Iterate& /*iterate*/,
         Evaluations& /*evaluations*/, double /*multipliers_threshold*/) {
      DEBUG << "The inverse Newton solver does not compute least-squares multipliers, keeping existing multipliers";
   }

   const Direction& InverseNewtonSolver::solve(Statistics& /*statistics*/, const Subproblem& subproblem, const Iterate& current_iterate,
         [[maybe_unused]] double trust_region_radius, const Vector<double>& /*initial_point*/, Evaluations& current_evaluations,
         const WarmstartInformation& /*warmstart_information*/) {
      assert(is_infinite(trust_region_radius));
      this->direction.reset();

      // store -gradient in this->rhs
      current_evaluations.evaluate_objective_gradient(subproblem.problem.model, current_iterate.primals);
      this->rhs = -current_evaluations.objective_gradient;

      // compute the Newton step d = -H⁻¹ g
      this->hessian_model.compute_inverse_hessian_vector_product(current_iterate.primals.data(),
         this->rhs.data(), this->direction.primals.data());
      direction.norm = norm_inf(view(direction.primals, 0, subproblem.problem.get_number_original_variables()));
      return this->direction;
   }

   bool InverseNewtonSolver::has_second_order_corrections() const {
      return false;
   }

   const Direction& InverseNewtonSolver::compute_second_order_correction(const Subproblem& /*subproblem*/,
         const Iterate& /*current_iterate*/, const Vector<double>& /*constraints_SOC*/) {
      throw std::runtime_error("No SOC implemented in InverseNewtonSolver");
   }

   const SolverWorkspace& InverseNewtonSolver::get_workspace() const {
      return this->workspace;
   }
} // namespace