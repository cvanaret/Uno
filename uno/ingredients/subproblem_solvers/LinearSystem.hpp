// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSYSTEM_H
#define UNO_LINEARSYSTEM_H

#include <cstddef>
#include "SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class NewtonSystem;

   class LinearSystem: public SolverWorkspace {
   public:
      LinearSystem() = default;
      ~LinearSystem() override = default;

      size_t dimension{0};
      size_t number_nonzeros{0};

      Vector<double> matrix_values{};
      Vector<double> rhs{};
      Vector<double> solution{};

      virtual void initialize_hessian(const Subproblem& subproblem) = 0;
      virtual void initialize_augmented_system(const Subproblem& subproblem) = 0;
      virtual void initialize(const NewtonSystem& newton_system) = 0;

      // result = matrix * vector (symmetric matrix), skipping the nonzeros in [excluded_begin, excluded_end) (e.g. the Schur
      // complement of a condensed Newton system). Both arrays have (at least) `dimension` entries
      virtual void compute_matrix_vector_product(const double* vector, double* result, size_t excluded_begin,
         size_t excluded_end) const = 0;
      void compute_matrix_vector_product(const double* vector, double* result) const {
         this->compute_matrix_vector_product(vector, result, 0, 0);
      }
   };
} // namespace

#endif // UNO_LINEARSYSTEM_H