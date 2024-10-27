// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVERFACTORY_H
#define UNO_LPSOLVERFACTORY_H

#include <memory>
#include <stdexcept>
#include "options/Options.hpp"
#include "solvers/LPSolver.hpp"

#ifdef HAS_BQPD
#include "solvers/BQPD/BQPDSolver.hpp"
#endif

namespace uno {
   class LPSolverFactory {
   public:
      static std::unique_ptr<LPSolver> create([[maybe_unused]] size_t number_variables, [[maybe_unused]] size_t number_constraints,
            [[maybe_unused]] size_t number_objective_gradient_nonzeros, [[maybe_unused]] size_t number_jacobian_nonzeros,
            [[maybe_unused]] const Options& options) {
         try {
            [[maybe_unused]] const std::string& LP_solver_name = options.get_string("LP_solver");
#ifdef HAS_BQPD
            if (LP_solver_name == "BQPD") {
               return std::make_unique<BQPDSolver>(number_variables, number_constraints, number_objective_gradient_nonzeros,
                     number_jacobian_nonzeros, 0, BQPDProblemType::LP, options);
            }
#endif
            std::string message = "The LP solver ";
            message.append(LP_solver_name).append(" is unknown").append("\n").append("The following values are available: ")
                  .append(join(LPSolverFactory::available_solvers(), ", "));
            throw std::invalid_argument(message);
         }
         catch (const std::out_of_range& exception) {
            std::string message = exception.what();
            message.append("\n").append("The following values are available: ").append(join(LPSolverFactory::available_solvers(), ", "));
            throw std::out_of_range(message);
         }
      }

      // return the list of available LP solvers
      static std::vector<std::string> available_solvers() {
         std::vector<std::string> solvers{};
#ifdef HAS_BQPD
         solvers.emplace_back("BQPD");
#endif
         return solvers;
      }
   };
} // namespace

#endif // UNO_LPSOLVERFACTORY_H
