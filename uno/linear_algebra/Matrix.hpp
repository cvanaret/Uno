// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIX_H
#define UNO_MATRIX_H

#include <cstddef>
#include <tuple>

namespace uno {
   template <typename IndexType>
   class Matrix {
   public:
      Matrix() = default;
      virtual ~Matrix() = default;

      [[nodiscard]] virtual std::tuple<IndexType, IndexType, double> operator[](size_t nonzero_index) const = 0;
   };
} // namespace

#endif // UNO_MATRIX_H