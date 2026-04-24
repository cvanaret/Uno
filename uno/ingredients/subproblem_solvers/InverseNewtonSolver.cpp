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
      this->rhs.resize(subproblem.number_variables);
   }

   void InverseNewtonSolver::solve(Statistics& /*statistics*/, const Subproblem& subproblem, [[maybe_unused]] double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& /*warmstart_information*/) {
      assert(is_infinite(trust_region_radius));

      // store -gradient in this->rhs
      current_evaluations.evaluate_objective_gradient(subproblem.problem.model, subproblem.current_iterate.primals);
      this->rhs = -current_evaluations.objective_gradient;

      // compute the Newton step d = -H⁻¹ g
      this->hessian_model.compute_inverse_hessian_vector_product(subproblem.current_iterate.primals.data(),
         this->rhs.data(), direction.primals.data());
   }

   SolverWorkspace& InverseNewtonSolver::get_workspace() {
      return this->workspace;
   }
} // namespace