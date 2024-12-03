// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

#include <vector>

namespace uno {
   // forward declarations
   class Direction;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   template <typename ElementType>
   class Vector;
   struct WarmstartInformation;

   /*! \class LPSolver
    * \brief LP solver
    *
    */
   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver() = default;

      virtual void solve_LP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
            const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
            const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H
