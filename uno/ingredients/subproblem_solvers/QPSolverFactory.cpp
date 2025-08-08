// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "QPSolverFactory.hpp"
#include "QPSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"

#ifdef HAS_BQPD
#include "ingredients/subproblem_solvers/BQPD/BQPDSolver.hpp"
#endif
#ifdef HAS_HIGHS
#include "ingredients/subproblem_solvers/HiGHS/HiGHSSolver.hpp"
#endif

namespace uno {
   std::unique_ptr<QPSolver> QPSolverFactory::create([[maybe_unused]] const Options& options) {
      try {
         [[maybe_unused]] const std::string& QP_solver_name = options.get_string("QP_solver");
#ifdef HAS_BQPD
         if (QP_solver_name == "BQPD") {
            return std::make_unique<BQPDSolver>(options);
         }
#endif
#ifdef HAS_HIGHS
         if (QP_solver_name == "HiGHS") {
            return std::make_unique<HiGHSSolver>(options);
         }
#endif
         std::string message = "The QP solver ";
         message.append(QP_solver_name).append(" is unknown").append("\n").append("The following values are available: ")
            .append(join(QPSolverFactory::available_solvers, ", "));
         throw std::invalid_argument(message);
      }
      catch (const std::out_of_range& exception) {
         std::string message = exception.what();
         message.append("\n").append("The following values are available: ").append(join(QPSolverFactory::available_solvers, ", "));
         throw std::out_of_range(message);
      }
   }
} // namespace
