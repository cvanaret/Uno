// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
#define UNO_SYMMETRICINDEFINITELINEARSYSTEM_H

#include "SymmetricMatrix.hpp"
#include "SparseStorageFactory.hpp"
#include "RectangularMatrix.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class SymmetricIndefiniteLinearSystem {
   public:
      SymmetricMatrix<size_t, ElementType> matrix{};
      Vector<ElementType> rhs{};
      Vector<ElementType> solution{};

      SymmetricIndefiniteLinearSystem() = default;
      SymmetricIndefiniteLinearSystem& operator=(SymmetricIndefiniteLinearSystem&& other) = default;

      void assemble_matrix(const SymmetricMatrix<size_t, double>& hessian, const RectangularMatrix<double>& constraint_jacobian,
         size_t number_variables, size_t number_constraints);
      void solve(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver);
   };

   template <typename ElementType>
   void SymmetricIndefiniteLinearSystem<ElementType>::assemble_matrix(const SymmetricMatrix<size_t, double>& hessian,
         const RectangularMatrix<double>& constraint_jacobian, size_t number_variables, size_t number_constraints) {
      this->matrix.set_dimension(number_variables + number_constraints);
      this->matrix.reset();
      // copy the Lagrangian Hessian in the top left block
      //size_t current_column = 0;
      for (const auto [row_index, column_index, element]: hessian) {
         // finalize all empty columns
         /*for (size_t column: Range(current_column, column_index)) {
            this->matrix.finalize_column(column);
            current_column++;
         }*/
         this->matrix.insert(element, row_index, column_index);
      }

      // Jacobian of general constraints
      for (size_t column_index: Range(number_constraints)) {
         for (const auto [row_index, derivative]: constraint_jacobian[column_index]) {
            this->matrix.insert(derivative, row_index, number_variables + column_index);
         }
         this->matrix.finalize_column(column_index);
      }
   }

   template <typename ElementType>
   void SymmetricIndefiniteLinearSystem<ElementType>::solve(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver) {
      linear_solver.solve_indefinite_system(this->matrix, this->rhs, this->solution);
   }
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
