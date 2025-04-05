// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include "HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   class LBFGSHessian : public HessianModel {
   public:
      LBFGSHessian(size_t dimension, size_t memory_size);
      ~LBFGSHessian() override = default;

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void notify_accepted_iterate(const Iterate& current_iterate, const Iterate& trial_iterate) override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;

   protected:
      const size_t dimension;
      const size_t memory_capacity;
      size_t current_memory_size{0};
      size_t current_available_slot{0};
      DenseMatrix<double> S_matrix;
      DenseMatrix<double> Y_matrix;
      DenseMatrix<double> L_matrix;
      std::vector<double> D_matrix; // D is diagonal
      DenseMatrix<double> M_matrix;
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H