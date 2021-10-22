#ifndef LINEARSOLVERFACTORY_H
#define LINEARSOLVERFACTORY_H

#include <memory>
#include "LinearSolver.hpp"

#ifdef HAS_MA57
#include "MA57Solver.hpp"
#endif
#ifdef HAS_PARDISO
#include "PardisoSolver.hpp"
#endif

class LinearSolverFactory {
public:
   static std::unique_ptr<LinearSolver> create(const std::string& linear_solver_name, size_t max_dimension, size_t max_number_nonzeros) {
#ifdef HAS_MA57
      if (linear_solver_name == "MA57") {
         return std::make_unique<MA57Solver>(max_dimension, max_number_nonzeros);
      }
#endif
/*
#ifdef HAS_PARDISO
      if (linear_solver_name == "PARDISO") {
         return std::make_unique<PardisoSolver>(max_dimension, max_number_nonzeros);
      }
#endif
 */
      throw std::invalid_argument("Linear solver name is unknown");
   }
};

#endif // LINEARSOLVERFACTORY_H
