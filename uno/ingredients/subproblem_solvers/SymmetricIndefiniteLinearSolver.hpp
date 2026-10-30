// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include "LinearSystem.hpp"

namespace uno {
   // forward declarations
   class LinearSystem;

   template <typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      SymmetricIndefiniteLinearSolver() = default;
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void initialize_memory() = 0;

      // solve with a single right-hand side: the right-hand side is read from the linear system
      // (get_linear_system().rhs) and the solution is written to `solution`
      virtual void solve_indefinite_system(ElementType* solution) = 0;

      // solve with `number_of_rhs` right-hand sides stored column-major (each column has length
      // get_linear_system().dimension): `rhs` is the block of right-hand sides and `solution` receives
      // the corresponding solutions (the two blocks may alias)
      virtual void solve_indefinite_system(const ElementType* rhs, ElementType* solution, size_t number_of_rhs) = 0;

      [[nodiscard]] virtual LinearSystem& get_linear_system() = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H