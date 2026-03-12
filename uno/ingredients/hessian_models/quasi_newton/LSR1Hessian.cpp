// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LSR1Hessian.hpp"
#include "linear_algebra/LAPACK_extension.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LSR1Hessian::LSR1Hessian(const Model& model, double objective_multiplier, const Options& options):
            QuasiNewtonHessian("L-SR1", model, objective_multiplier, options) {
   }

   bool LSR1Hessian::is_positive_definite() const {
      return false;
   }

   void LSR1Hessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|SR1|", Statistics::double_width - 2, 2, Statistics::column_order.at("|SR1|"));
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   void LSR1Hessian::notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   void LSR1Hessian::compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) {

   }
} // namespace