// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include <vector>
#include "solvers/LP/LPSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/RectangularMatrix.hpp"

/*! \class QPSolver
 * \brief QP solver
 *
 */
class QPSolver : public LPSolver {
public:
   QPSolver() = default;
   ~QPSolver() override = default;
   virtual Direction solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) = 0;
   Direction solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) override = 0;
};

#endif // UNO_QPSOLVER_H
