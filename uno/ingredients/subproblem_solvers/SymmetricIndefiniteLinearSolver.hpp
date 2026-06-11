// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include "LinearSystem.hpp"
#include "linear_algebra/View.hpp"

namespace uno {
   template <typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      SymmetricIndefiniteLinearSolver() = default;
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void initialize_memory() = 0;

      // solve with a single right-hand side: the right-hand side is read from the linear system
      // (get_linear_system().rhs) and the solution is written to `solution`
      virtual void solve_indefinite_system(View<ElementType> solution) = 0;

      // solve with `number_of_rhs` right-hand sides stored column-major (each column has length
      // get_linear_system().dimension): `rhs` is the block of right-hand sides and `solution` receives
      // the corresponding solutions (the two blocks may alias).
      // The default implementation loops over single right-hand side solves; solvers with native
      // multiple-RHS support (e.g. MA57, MUMPS, SSIDS) override this for efficiency
      virtual void solve_indefinite_system(const ElementType* rhs, ElementType* solution, size_t number_of_rhs) {
         LinearSystem& linear_system = this->get_linear_system();
         // solve for each right-hand side
         for (size_t column_index: Range(number_of_rhs)) {
            const size_t column_offset = column_index * linear_system.dimension;
            const auto rhs_column = view(rhs, column_offset, column_offset + linear_system.dimension);
            // copy the column into linear_system.rhs, then solve into the matching solution column
            linear_system.rhs.view() = rhs_column;
            auto rhs_solution = view(solution, column_offset, column_offset + linear_system.dimension);
            this->solve_indefinite_system(rhs_solution);
         }
      }

      [[nodiscard]] virtual LinearSystem& get_linear_system() = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H