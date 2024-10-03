// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ITERATIVESYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_ITERATIVESYMMETRICINDEFINITELINEARSOLVER_H

#include "solvers/SymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   template <typename IndexType, typename NumericalType, typename LinearOperator>
   class IterativeSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<IndexType, NumericalType, LinearOperator> {
   public:
      explicit IterativeSymmetricIndefiniteLinearSolver(size_t dimension) :
         SymmetricIndefiniteLinearSolver<IndexType, NumericalType, LinearOperator>(dimension) { };
      ~IterativeSymmetricIndefiniteLinearSolver() = default;

      void solve_indefinite_system(const SymmetricMatrix<IndexType, NumericalType>& /*matrix*/, const Vector<NumericalType>& /*rhs*/,
            Vector<NumericalType>& /*result*/) override {
         throw std::runtime_error("IterativeSymmetricIndefiniteLinearSolver: solve_indefinite_system with matrix not implemented yet.");
      }

      virtual void solve_indefinite_system(const LinearOperator& linear_operator, const Vector<NumericalType>& rhs, Vector<NumericalType>& result) = 0;
   };
} // namespace

#endif // UNO_ITERATIVESYMMETRICINDEFINITELINEARSOLVER_H