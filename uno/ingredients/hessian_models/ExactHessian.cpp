// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   // exact Hessian
   ExactHessian::ExactHessian(): HessianModel() {
   }

   void ExactHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const { }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& primal_variables,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(problem.number_variables);
      problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, hessian);
      this->evaluation_count++;
   }

   void ExactHessian::compute_hessian_vector_product(const OptimizationProblem& problem,
         const Vector<double>& vector, const Vector<double>& constraint_multipliers, Vector<double>& result) {
      problem.compute_hessian_vector_product(vector, constraint_multipliers, result);
      this->evaluation_count++;
   }
} // namespace