// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <vector>
#include <linear_algebra/SymmetricMatrix.hpp>

template <typename IndexType, typename ElementType>
class SymmetricIndefiniteLinearSolver {
public:
   explicit SymmetricIndefiniteLinearSolver(size_t max_dimension): max_dimension(max_dimension) {};
   virtual ~SymmetricIndefiniteLinearSolver() = default;

   virtual void solve_indefinite_system(const SymmetricMatrix<IndexType, ElementType>& matrix, const std::vector<ElementType>& rhs,
         std::vector<ElementType>& result) = 0;

protected:
   const size_t max_dimension;
};

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H
