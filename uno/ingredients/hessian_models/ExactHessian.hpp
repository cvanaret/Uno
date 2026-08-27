// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXACTHESSIAN_H
#define UNO_EXACTHESSIAN_H

#include "HessianModel.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Model;

   class HessianBuffer {
   public:
      HessianBuffer(const Model& model);
      ~HessianBuffer() = default;

      void evaluate_hessian(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, View<double> hessian_values);

   protected:
      Vector<double> hessian_buffer;
      bool first_evaluation{true};
   };

   class ExactHessian : public HessianModel {
   public:
      explicit ExactHessian(const Model& model);
      ~ExactHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize_statistics(Statistics& statistics) const override;
      void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, View<double> hessian_values) override;
      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, View<double> result) override;

   protected:
      const Model& model;
      HessianBuffer hessian_buffer;
   };
} // namespace

#endif // UNO_EXACTHESSIAN_H