// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "options/Options.hpp"

namespace uno {
   ZeroHessian::ZeroHessian(size_t dimension, const Options& options) :
         HessianModel(dimension, 0, options.get_string("sparse_format"), /* use_regularization = */false) { }

   void ZeroHessian::evaluate(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& /*primal_variables*/,
         const Vector<double>& /*constraint_multipliers*/) {
      this->hessian.set_dimension(problem.number_variables);
   }
}