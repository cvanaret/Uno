#ifndef LINEARSOLVERFACTORY_H
#define LINEARSOLVERFACTORY_H

#include <memory>
#include "LinearSolver.hpp"
#include "Problem.hpp"

class LinearSolverFactory {
public:
    static std::shared_ptr<LinearSolver> create(const std::string& linear_solver, int number_variables, int number_constraints, int maximum_number_nonzeros);
};

#endif // LINEARSOLVERFACTORY_H
