// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cstddef>
#include "LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"

namespace uno {
   LagrangeNewtonSubproblem::LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate,
         const Vector<double>& current_multipliers, const HessianModel& hessian_model, double trust_region_radius):
         number_variables(problem.number_variables), number_constraints(problem.number_constraints),
         problem(problem), current_iterate(current_iterate), current_multipliers(current_multipliers), hessian_model(hessian_model),
         trust_region_radius(trust_region_radius) { }

   void LagrangeNewtonSubproblem::evaluate_objective_gradient(SparseVector<double>& gradient) const {
      gradient.clear();
      this->problem.evaluate_objective_gradient(this->current_iterate, gradient);
   }

   void LagrangeNewtonSubproblem::evaluate_constraints(Vector<double>& constraints) const {
      this->problem.evaluate_constraints(this->current_iterate, constraints);
   }

   void LagrangeNewtonSubproblem::evaluate_constraint_jacobian(RectangularMatrix<double>& jacobian) const {
      this->problem.evaluate_constraint_jacobian(this->current_iterate, jacobian);
   }

   void LagrangeNewtonSubproblem::evaluate_lagrangian_hessian(const Vector<double>& x, SymmetricMatrix<size_t, double>& hessian) const {
      this->problem.evaluate_lagrangian_hessian(x, this->current_multipliers, hessian);
   }

   void LagrangeNewtonSubproblem::compute_hessian_vector_product(const Vector<double>& x, Vector<double>& result) const {
      this->problem.compute_hessian_vector_product(x, this->current_multipliers, result);
   }
} // namespace