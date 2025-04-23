// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

namespace uno {
   // zero Hessian
   class ZeroHessian : public HessianModel {
   public:
      ZeroHessian() = default;
      ~ZeroHessian() override = default;

      void initialize(const Model& model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      [[nodiscard]] size_t number_nonzeros(const Model& model) const override;
      void notify_accepted_iterate(const Iterate& current_iterate, const Iterate& trial_iterate) override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;
   };
} // namespace