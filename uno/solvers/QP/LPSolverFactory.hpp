#ifndef LPSOLVERFACTORY_H
#define LPSOLVERFACTORY_H

#include <memory>
#include "LPSolver.hpp"

template<class LPSolverType>
class LPSolverFactory;

// specialize the template factory with concrete solver types
#ifdef HAS_BQPD
#include "BQPDSolver.hpp"

template<>
class LPSolverFactory<BQPDSolver> {
public:
   static std::unique_ptr<BQPDSolver> create(size_t number_variables, size_t number_constraints) {
      return std::make_unique<BQPDSolver>(number_variables, number_constraints, 0, false);
   }
};
#endif

#endif // LPSOLVERFACTORY_H
