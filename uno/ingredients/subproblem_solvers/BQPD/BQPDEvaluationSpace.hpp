// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDEVALUATIONSPACE_H
#define UNO_BQPDEVALUATIONSPACE_H

#include <cstddef>
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   class Subproblem;
   class WarmstartInformation;

   class BQPDEvaluationSpace: public EvaluationSpace {
   public:
      BQPDEvaluationSpace() = default;
      ~BQPDEvaluationSpace() override = default;

      void initialize(const Subproblem& subproblem);

      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector,
         Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate, const WarmstartInformation& warmstart_information);

      Vector<double> constraints{};
      Vector<double> gradients{};
      Vector<uno_int> gradient_sparsity{};
      // COO constraint Jacobian
      Vector<uno_int> jacobian_row_indices{};
      Vector<uno_int> jacobian_column_indices{};
      Vector<double> jacobian_values{};
      Vector<size_t> permutation_vector{};
      // COO Hessian
      Vector<uno_int> hessian_row_indices{};
      Vector<uno_int> hessian_column_indices{};
      Vector<double> hessian_values{};
      bool evaluate_hessian{false};

   protected:
      void compute_gradients_sparsity(const Subproblem& subproblem);
   };
} // namespace

#endif // UNO_BQPDEVALUATIONSPACE_H