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
      void initialize(const NewtonSystem& newton_system) override;
      using LinearSystem::compute_matrix_vector_product;
      void compute_matrix_vector_product(const double* vector, double* result, size_t excluded_begin,
         size_t excluded_end) const override;

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& vector) const override;

      // symmetric matrix (Hessian or augmented system)
      std::vector<uno_int> matrix_row_indices{};
      std::vector<uno_int> matrix_column_indices{};

   protected:
      const int solver_indexing;
   };
} // namespace

#endif // UNO_COOLINEARSYSTEM_H