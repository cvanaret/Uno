#ifndef LPSOLVER_H
#define LPSOLVER_H

#include <vector>
#include "SubproblemSolution.hpp"

/*! \class LPSolver
 * \brief LP solver
 *
 */
class LPSolver {
public:

    virtual ~LPSolver() {
    };
    virtual SubproblemSolution solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x) = 0;
};

#endif // LPSOLVER_H
