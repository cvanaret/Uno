// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSEVALUATIONSPACE_H
#define UNO_HIGHSEVALUATIONSPACE_H

#include <cstddef>
#include "optimization/EvaluationSpace.hpp"
#include "Highs.h"
#include "linear_algebra/Vector.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   class HiGHSEvaluationSpace: public EvaluationSpace {
   public:
      HiGHSEvaluationSpace() = default;
      ~HiGHSEvaluationSpace() override = default;

      void initialize_memory(const Subproblem& subproblem);

      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector,
         Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void evaluate_functions(Statistics& statistics, const Subproblem& subproblem, const WarmstartInformation& warmstart_information);

      HighsModel model;
      Vector<double> constraints{};
      Vector<double> linear_objective{};
      // constraint Jacobian in COO format
      Vector<uno_int> jacobian_row_indices{};
      Vector<uno_int> jacobian_column_indices{};
      Vector<double> jacobian_values{};
      Vector<size_t> jacobian_permutation_vector{};
      // Lagrangian Hessian in COO format
      Vector<uno_int> hessian_row_indices{};
      Vector<uno_int> hessian_column_indices{};
      Vector<double> hessian_values{};
      Vector<size_t> hessian_permutation_vector{};

   protected:
      void compute_jacobian_sparsity(const Subproblem& subproblem);
      void compute_hessian_sparsity(const Subproblem& subproblem);
   };
} // namespace

#endif // UNO_HIGHSEVALUATIONSPACE_H