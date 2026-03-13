// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LSR1Hessian.hpp"
#include "linear_algebra/LAPACK_extension.hpp"
#include "model/Model.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LSR1Hessian::LSR1Hessian(const Model& model, double objective_multiplier, const Options& options):
            QuasiNewtonHessian("L-SR1", model, objective_multiplier, options),
         N(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size) {
   }

   bool LSR1Hessian::is_positive_definite() const {
      return false;
   }

   void LSR1Hessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|SR1|", Statistics::double_width - 3, 2, Statistics::column_order.at("|SR1|"));
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   void LSR1Hessian::notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      DEBUG << "\n*** Updating the SR1 memory at slot " << this->current_index << '\n';
      // update the matrices S and Y
      this->update_limited_memory(current_iterate, trial_iterate, evaluation_cache);

      // build the N matrix
      // TODO
      auto Nk = this->N.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);
      // compute an LDL' (signed Cholesky) factorization
      // TODO because the factorization may fail, make a copy
      const double near_zero_pivot_tolerance = 0.; // TODO
      const bool ldlt_success = ldlt_nopiv_lvl2_rightlooking(Nk.data(), Nk.number_rows, Nk.leading_dimension, near_zero_pivot_tolerance);
      if (ldlt_success) {
         this->validate_update();
      }
      else {
         DEBUG << "Skipping the update\n";
      }
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   // Hessian-vector product where the Hessian approximation is Bk = B0 + U P⁻¹ Uᵀ and B0 = delta I
   // Bk v = (B0 + U P⁻¹ Uᵀ) v = delta v + U (P⁻¹ (Uᵀ v))
   void LSR1Hessian::compute_hessian_vector_product(const double* /*x*/, const double* vector,
         double objective_multiplier, const Vector<double>& /*constraint_multipliers*/, double* result) {
      if (objective_multiplier != this->fixed_objective_multiplier) {
         throw std::runtime_error("The L-BFGS Hessian model was initialized with a different objective multiplier");
      }

      // update_limited_memory() has updated the limited memory. A recomputation of the Hessian representation may be required
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // diagonal contribution delta*I
      for (size_t variable_index: Range(this->model.number_variables)) {
         result[variable_index] = this->delta * vector[variable_index];
      }

      // rank-1 contribution: U in R^{n x m}
      // work on each column of U (P⁻¹ (Uᵀ v))
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         double U_coefficient = dot(current_U_column, vector);
         U_coefficient *= 1.; // TODO
         // result += coefficient * column(P⁻¹) * current_column
         blas1::add(this->model.number_variables, U_coefficient, current_U_column.data(), result);
      }
   }

   // protected member functions

   void LSR1Hessian::recompute_hessian_representation() {
      // TODO
   }
} // namespace