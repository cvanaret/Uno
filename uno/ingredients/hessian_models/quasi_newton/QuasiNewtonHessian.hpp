// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QUASINEWTONHESSIAN_H
#define UNO_QUASINEWTONHESSIAN_H

#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;

   // express the Hessian approximation at iteration k by a low-rank update
   class QuasiNewtonHessian: public HessianModel {
   public:
      QuasiNewtonHessian(const std::string_view name, const Model& model, double objective_multiplier, const Options& options);
      ~QuasiNewtonHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const override;

      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;

   protected:
      const Model& model;
      const double fixed_objective_multiplier;
      const size_t memory_size; // user defined
      size_t number_entries_in_memory{0}; // 0 <= number_entries_in_memory <= memory_size
      size_t current_index{0}; // 0 <= current_index < memory_size
      // limited memory
      DenseMatrix<double> S;
      DenseMatrix<double> Y;
      Vector<double> current_lagrangian_gradient;
      Vector<double> trial_lagrangian_gradient;
      bool hessian_recomputation_required{false};

      void update_limited_memory(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache);
      void validate_update();
      void update_S(const Iterate& current_iterate, const Iterate& trial_iterate);
      void update_Y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache);
   };
} // namespace

#endif // UNO_QUASINEWTONHESSIAN_H