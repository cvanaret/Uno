// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include "linear_algebra/View.hpp"

namespace uno {
   // forward declarations
   class LinearSystem;

   template <typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      SymmetricIndefiniteLinearSolver() = default;
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void initialize_memory() = 0;
      virtual void solve_indefinite_system(View<double> result) = 0;
      [[nodiscard]] virtual LinearSystem& get_linear_system() = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H