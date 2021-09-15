#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"

template<class QPSolverType>
class QPSolverFactory;

// specialize the template factory with concrete solver types
#ifdef HAS_BQPD
#include "BQPDSolver.hpp"

template<>
class QPSolverFactory<BQPDSolver> {
public:
   static std::unique_ptr<QPSolver<typename BQPDSolver::matrix_type>> create(size_t number_variables, size_t number_constraints,
         size_t maximum_number_nonzeros, bool quadratic_programming) {
      return std::make_unique<BQPDSolver>(number_variables, number_constraints, maximum_number_nonzeros, quadratic_programming);
   }
};
#endif

#endif // QPSOLVERFACTORY_H
