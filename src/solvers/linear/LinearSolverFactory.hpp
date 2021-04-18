#ifndef LINEARSOLVERFACTORY_H
#define LINEARSOLVERFACTORY_H

#include <memory>
#include "LinearSolver.hpp"
#include "Problem.hpp"

class LinearSolverFactory {
public:
    static std::unique_ptr<LinearSolver> create(const std::string& linear_solver);
};

#endif // LINEARSOLVERFACTORY_H
