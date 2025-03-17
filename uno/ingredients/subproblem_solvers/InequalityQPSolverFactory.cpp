// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "InequalityQPSolverFactory.hpp"
#include "InequalityQPSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"

#ifdef HAS_BQPD
#include "ingredients/subproblem_solvers/BQPD/BQPDSolver.hpp"
#endif

namespace uno {
   std::unique_ptr<InequalityQPSolver> InequalityQPSolverFactory::create([[maybe_unused]] size_t number_variables, [[maybe_unused]] size_t number_constraints,
         [[maybe_unused]] size_t number_objective_gradient_nonzeros, [[maybe_unused]] size_t number_jacobian_nonzeros,
         [[maybe_unused]] size_t number_hessian_nonzeros, [[maybe_unused]] const Options& options) {
      try {
         [[maybe_unused]] const std::string& QP_solver_name = options.get_string("inequality_QP_solver");
#ifdef HAS_BQPD
         if (QP_solver_name == "BQPD") {
            return std::make_unique<BQPDSolver>(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
                  number_hessian_nonzeros, options);
         }
#endif
         std::string message = "The inequality QP solver ";
         message.append(QP_solver_name).append(" is unknown").append("\n").append("The following values are available: ")
               .append(join(InequalityQPSolverFactory::available_solvers(), ", "));
         throw std::invalid_argument(message);
      }
      catch (const std::out_of_range& exception) {
         std::string message = exception.what();
         message.append("\n").append("The following values are available: ").append(join(InequalityQPSolverFactory::available_solvers(), ", "));
         throw std::out_of_range(message);
      }
   }

   // return the list of available QP solvers
   std::vector<std::string> InequalityQPSolverFactory::available_solvers() {
      std::vector<std::string> solvers{};
#ifdef HAS_BQPD
      solvers.emplace_back("BQPD");
#endif
      return solvers;
   }
} // namespace
