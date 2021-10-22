#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"

#ifdef HAS_BQPD
#include "BQPDSolver.hpp"
#endif

class QPSolverFactory {
public:
   static std::unique_ptr<QPSolver> create(const std::string& QP_solver_name, size_t number_variables, size_t number_constraints,
         size_t maximum_number_nonzeros, bool quadratic_programming) {
#ifdef HAS_BQPD
      if (QP_solver_name == "BQPD") {
         return std::make_unique<BQPDSolver>(number_variables, number_constraints, maximum_number_nonzeros, quadratic_programming);
      }
#endif
      throw std::invalid_argument("QP solver name is unknown");
   }
};

#endif // QPSOLVERFACTORY_H
