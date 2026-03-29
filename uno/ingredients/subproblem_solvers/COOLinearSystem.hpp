// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOLINEARSYSTEM_H
#define UNO_COOLINEARSYSTEM_H

#include <vector>
#include "LinearSystem.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   class COOLinearSystem: public LinearSystem {
   public:
      explicit COOLinearSystem(int solver_indexing);
      ~COOLinearSystem() override = default;

      void initialize_hessian(const Subproblem& subproblem) override;
      void initialize_augmented_system(const Subproblem& subproblem) override;

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const override;

      // symmetric matrix (Hessian or augmented system)
      std::vector<uno_int> matrix_row_indices{};
      std::vector<uno_int> matrix_column_indices{};

   protected:
      const int solver_indexing;
   };
} // namespace

#endif // UNO_COOLINEARSYSTEM_H