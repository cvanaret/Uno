#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"
#include "Problem.hpp"

#ifdef HAS_BQPD
#include "BQPDSolver.hpp"
#endif

template<class SparseSymmetricMatrix>
class QPSolverFactory;

// specialize the template factory with concrete matrix types (CSC and COO)
template<>
class QPSolverFactory<CSCSymmetricMatrix> {
public:
   static std::unique_ptr<QPSolver<CSCSymmetricMatrix> > create(const std::string& QP_solver_name, size_t number_variables, size_t number_constraints,
         size_t maximum_number_nonzeros, bool quadratic_programming) {
#ifdef HAS_BQPD
      if (QP_solver_name == "BQPD") {
         return std::make_unique<BQPDSolver>(number_variables, number_constraints, maximum_number_nonzeros, quadratic_programming);
      }
#endif
      throw std::invalid_argument("CSC QP solver " + QP_solver_name + " does not exist.");
   }
};

template<>
class QPSolverFactory<COOSymmetricMatrix> {
public:
   static std::unique_ptr<QPSolver<COOSymmetricMatrix> > create(const std::string& QP_solver_name, size_t /*number_variables*/, size_t
   /*number_constraints*/, size_t /*maximum_number_nonzeros*/, bool /*quadratic_programming*/) {
      throw std::invalid_argument("COO QP solver " + QP_solver_name + " does not exist.");
   }
};

#endif // QPSOLVERFACTORY_H
