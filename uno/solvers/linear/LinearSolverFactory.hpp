#ifndef LINEARSOLVERFACTORY_H
#define LINEARSOLVERFACTORY_H

#include <memory>
#include "LinearSolver.hpp"

template<typename LinearSolverType>
class LinearSolverFactory;

// specialize the template factory with concrete matrix types (CSC and COO)
/*
#ifdef HAS_PARDISO
#include "PardisoSolver.hpp"

template<>
class LinearSolverFactory<PardisoSolver> {
public:
   static std::unique_ptr<LinearSolver<typename PardisoSolver::matrix_type>> create(size_t max_dimension) {
      return std::make_unique<PardisoSolver>(max_dimension);
   }
};
#endif
 */

#ifdef HAS_MA57
#include "MA57Solver.hpp"

template<>
class LinearSolverFactory<MA57Solver> {
public:
   static std::unique_ptr<LinearSolver<typename MA57Solver::matrix_type>> create(size_t max_dimension) {
      return std::make_unique<MA57Solver>(max_dimension);
   }
};
#endif

#endif // LINEARSOLVERFACTORY_H
