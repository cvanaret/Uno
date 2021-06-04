#ifndef LPSOLVER_H
#define LPSOLVER_H

#include <vector>
#include "Direction.hpp"

/*! \class LPSolver
 * \brief LP solver
 *
 */
class LPSolver {
public:

    virtual ~LPSolver() {
    };
    virtual Direction solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, SparseVector& linear_objective, std::vector<SparseVector>& constraints_jacobian, std::vector<double>& x) = 0;
};

#endif // LPSOLVER_H
