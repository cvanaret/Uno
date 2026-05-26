// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONS_H
#define UNO_EVALUATIONS_H

#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // forward declarations
   class COOSparsity;
   class Model;

   class Evaluations {
   public:
      double objective{INF<double>};
      Vector<double> constraints;
      Vector<double> objective_gradient;
      Vector<double> jacobian_values;
      const COOSparsity* jacobian_sparsity; // Jacobian sparsity computed in EvaluationCache

      // flags to perform lazy computations
      bool is_objective_computed{false};
      bool are_constraints_computed{false};
      bool is_objective_gradient_computed{false};
      bool is_jacobian_computed{false};

      Evaluations(const Model& model, const COOSparsity* jacobian_sparsity);
      Evaluations(const Evaluations& other) = default;
      Evaluations(Evaluations&& other) noexcept = default;
      Evaluations& operator=(const Evaluations& other) = default;
      Evaluations& operator=(Evaluations&& other) noexcept = default;

      void evaluate_objective(const Model& model, const Vector<double>& primals);
      void evaluate_constraints(const Model& model, const Vector<double>& primals);
      void evaluate_objective_gradient(const Model& model, const Vector<double>& primals);
      void evaluate_jacobian(const Model& model, const Vector<double>& primals);

      void compute_jacobian_vector_product(const Model& model, const double* vector, double* result) const;
      void compute_jacobian_transposed_vector_product(const Model& model, const double* vector, double* result) const;

      // reset the evaluation flags
      void reset();
   };
} // namespace

#endif // UNO_EVALUATIONS_H