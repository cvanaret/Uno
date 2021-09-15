#ifndef LPSOLVER_H
#define LPSOLVER_H

#include <vector>
#include "ingredients/subproblem/Direction.hpp"
#include "linear_algebra/SparseVector.hpp"

/*! \class LPSolver
 * \brief LP solver
 *
 */
class LPSolver {
public:
   LPSolver() = default;
   virtual ~LPSolver() = default;
   virtual Direction solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
   linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const std::vector<double>& initial_point) = 0;
};

#endif // LPSOLVER_H
