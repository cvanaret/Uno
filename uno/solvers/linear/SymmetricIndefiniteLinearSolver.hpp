// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include <tuple>

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

   template <typename ElementType>
   class Vector;

   template <typename IndexType, typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      explicit SymmetricIndefiniteLinearSolver(size_t dimension) : dimension(dimension) {
      };
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void factorize(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
      virtual void do_symbolic_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
      virtual void do_numerical_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
      virtual void solve_indefinite_system(const SymmetricMatrix<IndexType, ElementType>& matrix, const Vector<ElementType>& rhs,
            Vector<ElementType>& result) = 0;

      [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;

   protected:
      const size_t dimension;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H