// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include <stdexcept>
#include <string>
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

#ifdef HAS_MA57
#include "solvers/MA57/MA57Solver.hpp"
#endif

#ifdef HAS_MUMPS
#include "solvers/MUMPS/MUMPSSolver.hpp"
#endif

namespace uno {
   class SymmetricIndefiniteLinearSolverFactory {
   public:
      static std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> create([[maybe_unused]] size_t dimension,
            [[maybe_unused]] size_t number_nonzeros, const Options& options) {
         try {
            [[maybe_unused]] const std::string& linear_solver_name = options.get_string("linear_solver");
#ifdef HAS_MA57
            if (linear_solver_name == "MA57") {
               return std::make_unique<MA57Solver>(dimension, number_nonzeros);
            }
#endif
#ifdef HAS_MUMPS
            if (linear_solver_name == "MUMPS") {
               return std::make_unique<MUMPSSolver>(dimension, number_nonzeros);
            }
#endif
            std::string message = "The linear solver ";
            message.append(linear_solver_name).append(" is unknown").append("\n").append("The following values are available: ")
                  .append(join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", "));
            throw std::invalid_argument(message);
         }
         catch (const std::out_of_range& exception) {
            std::string message = exception.what();
            message.append("\n").append("The following values are available: ")
                  .append(join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", "));
            throw std::out_of_range(message);
         }
      }

      // return the list of available solvers
      static std::vector<std::string> available_solvers() {
         std::vector<std::string> solvers{};
#ifdef HAS_MA57
         solvers.emplace_back("MA57");
#endif
#ifdef HAS_MUMPS
         solvers.emplace_back("MUMPS");
#endif
         return solvers;
      }
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H
