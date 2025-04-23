// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"

namespace uno {
   void ExactHessian::initialize(const Model& /*model*/) {
   }

   void ExactHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
      // do nothing
   }

   size_t ExactHessian::number_nonzeros(const Model& model) const {
      return model.number_hessian_nonzeros();
   }
   
   void ExactHessian::notify_accepted_iterate(const Model& /*model*/, Iterate& /*current_iterate*/, Iterate& /*trial_iterate*/) {
      // do nothing
   }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(model.number_variables);
      model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian);
      this->evaluation_count++;
   }

   void ExactHessian::compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) {
      model.compute_hessian_vector_product(vector, objective_multiplier, constraint_multipliers, result);
      // TODO fix the counter
      this->evaluation_count++;
   }
} // namespace