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

template <typename SparseSymmetricMatrix>
class LinearSolverFactory;

// specialize the template factory with concrete matrix types (CSC and COO)
template<>
class LinearSolverFactory<CSCSymmetricMatrix> {
public:
   static std::unique_ptr<LinearSolver<CSCSymmetricMatrix> > create(const std::string& linear_solver_name) {
#ifdef HAS_PARDISO
      if (linear_solver_name == "PARDISO") {
         return std::make_unique<PardisoSolver>();
      }
#endif
      throw std::invalid_argument("CSC linear solver " + linear_solver_name + " does not exist.");
   }
};

template<>
class LinearSolverFactory<COOSymmetricMatrix> {
public:
   static std::unique_ptr<LinearSolver<COOSymmetricMatrix> > create(const std::string& linear_solver_name) {
#ifdef HAS_MA57
      if (linear_solver_name == "MA57") {
         return std::make_unique<MA57Solver>();
      }
#endif
      throw std::invalid_argument("COO linear solver " + linear_solver_name + " does not exist.");
   }
};

#endif // LINEARSOLVERFACTORY_H
