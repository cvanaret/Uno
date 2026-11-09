// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "options/Options.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options):
         problem(problem),
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))) {
   }
} // namespace