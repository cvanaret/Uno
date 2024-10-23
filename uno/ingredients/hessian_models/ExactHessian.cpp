// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "options/Options.hpp"

namespace uno {
   // exact Hessian
   ExactHessian::ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options) :
         HessianModel(dimension, maximum_number_nonzeros, options.get_string("sparse_format"), /* use_regularization = */false) {
   }

   void ExactHessian::evaluate(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& primal_variables,
         const Vector<double>& constraint_multipliers) {
      // evaluate Lagrangian Hessian
      this->hessian.set_dimension(problem.number_variables);
      problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, this->hessian);
      this->evaluation_count++;
   }
} // namespace