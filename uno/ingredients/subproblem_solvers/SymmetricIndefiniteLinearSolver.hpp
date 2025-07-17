// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>

namespace uno {
   // forward declarations
   class Direction;
   class Statistics;
   class Subproblem;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   template <typename IndexType, typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      SymmetricIndefiniteLinearSolver() = default;
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void initialize_memory(size_t number_variables, size_t number_constraints, size_t number_hessian_nonzeros,
         size_t regularization_size) = 0;

      virtual void solve_indefinite_system(const SymmetricMatrix<IndexType, ElementType>& matrix, const Vector<ElementType>& rhs,
         Vector<ElementType>& result) = 0;
      virtual void solve_equality_constrained_subproblem(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H