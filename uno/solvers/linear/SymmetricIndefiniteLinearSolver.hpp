// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include <tuple>

namespace uno {
   // forward declarations
   template <typename IndexType, typename NumericalType>
   class SymmetricMatrix;

   template <typename NumericalType>
   class Vector;

   template <typename IndexType, typename NumericalType, typename LinearOperator>
   class SymmetricIndefiniteLinearSolver {
   public:
      explicit SymmetricIndefiniteLinearSolver(size_t dimension) : dimension(dimension) { };
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void solve_indefinite_system(const SymmetricMatrix<IndexType, NumericalType>& matrix, const Vector<NumericalType>& rhs,
            Vector<NumericalType>& result) = 0;
      virtual void solve_indefinite_system(const LinearOperator& linear_operator, const Vector<NumericalType>& rhs, Vector<NumericalType>& result) = 0;

   protected:
      const size_t dimension;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H