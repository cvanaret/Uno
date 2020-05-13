#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include "Matrix.hpp"

/*! \class LinearSolver
 * \brief QP solver
 *
 */
class LinearSolver {
public:

    virtual ~LinearSolver() {
    };
    virtual void solve(COOMatrix& matrix, std::vector<double>& rhs) = 0;
};

#endif // LINEARSOLVER_H
