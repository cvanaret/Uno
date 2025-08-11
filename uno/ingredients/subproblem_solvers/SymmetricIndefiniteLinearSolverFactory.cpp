// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/Vector.hpp"

#if defined(HAS_HSL) || defined(HAS_MA57)
#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#endif

#if defined(HAS_HSL) || defined(HAS_MA27)
#include "ingredients/subproblem_solvers/MA27/MA27Solver.hpp"
#endif

#ifdef HAS_HSL
namespace uno {
   extern "C" {
      bool LIBHSL_isfunctional();
   }
}
#endif

#ifdef HAS_MUMPS
#include "ingredients/subproblem_solvers/MUMPS/MUMPSSolver.hpp"
#endif

namespace uno {
   std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> SymmetricIndefiniteLinearSolverFactory::create(const std::string& linear_solver) {
#if defined(HAS_HSL) || defined(HAS_MA57)
      if (linear_solver == "MA57"
   #ifdef HAS_HSL
         && LIBHSL_isfunctional()
   #endif
            ) {
         return std::make_unique<MA57Solver>();
      }
#endif

#if defined(HAS_HSL) || defined(HAS_MA27)
      if (linear_solver == "MA27"
   # ifdef HAS_HSL
         && LIBHSL_isfunctional()
   # endif
      ) {
         return std::make_unique<MA27Solver>();
      }
#endif // HAS_HSL || HAS_MA27

#ifdef HAS_MUMPS
      if (linear_solver == "MUMPS") {
         return std::make_unique<MUMPSSolver>();
      }
#endif
      std::string message = "The linear solver ";
      message.append(linear_solver).append(" is unknown").append("\n").append("The following values are available: ")
            .append(join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", "));
      throw std::invalid_argument(message);
   }

   // return the list of available solvers
   std::vector<std::string> SymmetricIndefiniteLinearSolverFactory::available_solvers() {
      std::vector<std::string> solvers{};
#ifdef HAS_HSL
      if (LIBHSL_isfunctional()) {
            solvers.emplace_back("MA57");
            solvers.emplace_back("MA27");
         }
#else
   #ifdef HAS_MA57
      solvers.emplace_back("MA57");
   #endif
   #ifdef HAS_MA27
      solvers.emplace_back("MA27");
   #endif
#endif

#ifdef HAS_MUMPS
      solvers.emplace_back("MUMPS");
#endif
      return solvers;
   }
} // namespace
