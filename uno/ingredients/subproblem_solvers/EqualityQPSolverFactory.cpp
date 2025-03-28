// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "EqualityQPSolverFactory.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"

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
   std::unique_ptr<EqualityQPSolver<size_t, double>> EqualityQPSolverFactory::create(size_t number_variables,
         size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) {
      try {
         [[maybe_unused]] const std::string& solver_name = options.get_string("equality_QP_solver");
#if defined(HAS_HSL) || defined(HAS_MA57)
         if (solver_name == "MA57"
   #ifdef HAS_HSL
            && LIBHSL_isfunctional()
   #endif
               ) {
            return std::make_unique<MA57Solver>(number_variables, number_constraints, number_jacobian_nonzeros, number_hessian_nonzeros);
         }
#endif

#if defined(HAS_HSL) || defined(HAS_MA27)
         if (solver_name == "MA27"
   # ifdef HAS_HSL
            && LIBHSL_isfunctional()         
   # endif
         ) {
            return std::make_unique<MA27Solver>(number_variables, number_constraints, number_jacobian_nonzeros, number_hessian_nonzeros);
         }
#endif // HAS_HSL || HAS_MA27

#ifdef HAS_MUMPS
         if (solver_name == "MUMPS") {
            return std::make_unique<MUMPSSolver>(number_variables, number_constraints, number_jacobian_nonzeros, number_hessian_nonzeros);
         }
#endif
         std::string message = "The equality QP solver ";
         message.append(solver_name).append(" is unknown").append("\n").append("The following values are available: ")
               .append(join(EqualityQPSolverFactory::available_solvers(), ", "));
         throw std::invalid_argument(message);
      }
      catch (const std::out_of_range& exception) {
         std::string message = exception.what();
         message.append("\n").append("The following values are available: ")
               .append(join(EqualityQPSolverFactory::available_solvers(), ", "));
         throw std::out_of_range(message);
      }
   }

   // return the list of available solvers
   std::vector<std::string> EqualityQPSolverFactory::available_solvers() {
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
