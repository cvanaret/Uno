#ifndef QPSOLVER_H
#define QPSOLVER_H

#include <vector>
#include "LPSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

/*! \class QPSolver
 * \brief QP solver
 *
 */
class QPSolver : public LPSolver {
public:
   QPSolver() = default;
   ~QPSolver() override = default;
   virtual Direction solve_QP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds,
         const SparseVector<double>& linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const SymmetricMatrix&
         hessian, const std::vector<double>& initial_point) = 0;
   Direction solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
         linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const std::vector<double>& initial_point) override = 0;
};

#endif // QPSOLVER_H
