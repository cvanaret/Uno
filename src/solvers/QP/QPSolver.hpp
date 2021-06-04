#ifndef QPSOLVER_H
#define QPSOLVER_H

#include <vector>
#include "LPSolver.hpp"
#include "Direction.hpp"
#include "Matrix.hpp"

/*! \class QPSolver
 * \brief QP solver
 *
 */
class QPSolver : public LPSolver {
public:

    virtual ~QPSolver() {
    };
    virtual Direction solve_QP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, SparseVector& linear_objective, std::vector<SparseVector>& constraints_jacobian, CSCMatrix& hessian, std::vector<double>& x0) = 0;
    virtual Direction solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, SparseVector& linear_objective, std::vector<SparseVector>& constraints_jacobian, std::vector<double>& x0) = 0;
};

#endif // QPSOLVER_H
