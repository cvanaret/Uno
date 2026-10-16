// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIAGONALMATRIX_H
#define UNO_DIAGONALMATRIX_H

#include <cstddef>
#include <utility>
#include "symbolic/symbolic_traits.hpp"

namespace uno {
   // square diagonal matrix of a given dimension, given as a dematerialized vector (indices, component function):
   // the diagonal entry (i, i) is component_function(i) for i in indices, all other entries are zero.
   // The entries are never stored: they are computed when the matrix is evaluated or applied
   template <typename Indices, typename ComponentFunction>
   class DiagonalMatrix: public SymbolicExpression {
   public:
      using value_type = double;

      DiagonalMatrix(size_t dimension, Indices&& indices, ComponentFunction&& component_function):
            dimension(dimension), indices(std::forward<Indices>(indices)),
            component_function(std::forward<ComponentFunction>(component_function)) { }

      [[nodiscard]] size_t number_rows() const { return this->dimension; }
      [[nodiscard]] size_t number_columns() const { return this->dimension; }
      [[nodiscard]] const auto& get_indices() const { return this->indices; }
      [[nodiscard]] double entry(size_t index) const { return this->component_function(index); }

   protected:
      const size_t dimension;
      storage_t<Indices> indices;
      storage_t<ComponentFunction> component_function;
   };

   // free function
   template <typename Indices, typename ComponentFunction>
   DiagonalMatrix<Indices, ComponentFunction> diagonal_matrix(size_t dimension, Indices&& indices, ComponentFunction&& component_function) {
      return {dimension, std::forward<Indices>(indices), std::forward<ComponentFunction>(component_function)};
   }
} // namespace

#endif // UNO_DIAGONALMATRIX_H
