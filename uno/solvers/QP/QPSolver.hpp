// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include <vector>
#include "solvers/LP/LPSolver.hpp"

// forward declarations
class Direction;
template <typename ElementType>
class RectangularMatrix;
template <typename ElementType>
class SparseVector;
template <typename ElementType>
class SymmetricMatrix;
struct WarmstartInformation;

/*! \class QPSolver
 * \brief QP solver
 *
 */
class QPSolver : public LPSolver {
public:
   QPSolver();
   ~QPSolver() override = default;

   void solve_LP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) override = 0;

   virtual void solve_QP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) = 0;
};

inline QPSolver::QPSolver(): LPSolver() {
}

#endif // UNO_QPSOLVER_H
