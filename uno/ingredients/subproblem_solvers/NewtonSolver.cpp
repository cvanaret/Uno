// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "NewtonSolver.hpp"
#include "ingredients/hessian_models/quasi_newton/InverseLBFGSHessian.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   NewtonSolver::NewtonSolver(InverseLBFGSHessian& hessian_model): hessian_model(hessian_model) {
   }

   void NewtonSolver::initialize_memory(const Subproblem& subproblem) {
      this->rhs.resize(subproblem.number_variables);
   }

   void NewtonSolver::solve(Statistics& /*statistics*/, const Subproblem& subproblem, double trust_region_radius,
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

   SolverWorkspace& NewtonSolver::get_workspace() {
      return this->workspace;
   }
} // namespace