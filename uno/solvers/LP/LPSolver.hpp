// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

#include <vector>
#include "ingredients/subproblem/Direction.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "optimization/Model.hpp"
#include "optimization/WarmstartInformation.hpp"

/*! \class LPSolver
 * \brief LP solver
 *
 */
class LPSolver {
public:
   LPSolver() = default;
   virtual ~LPSolver() = default;
   virtual Direction solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) = 0;
};

#endif // UNO_LPSOLVER_H
