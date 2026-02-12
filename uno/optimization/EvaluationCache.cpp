#include "EvaluationCache.hpp"
#include "model/Model.hpp"

namespace uno {
   EvaluationCache::EvaluationCache(const Model& model):
         number_jacobian_nonzeros(model.number_jacobian_nonzeros()),
         jacobian_sparsity(this->number_jacobian_nonzeros),
         // pass (a const reference of) the Jacobian sparsity to both evaluations
         current_evaluations(model, &this->jacobian_sparsity),
         trial_evaluations(model, &this->jacobian_sparsity) {
      // compute the Jacobian sparsity
      this->jacobian_sparsity.row_indices.resize(this->number_jacobian_nonzeros);
      this->jacobian_sparsity.column_indices.resize(this->number_jacobian_nonzeros);
      model.compute_jacobian_sparsity(this->jacobian_sparsity.row_indices.data(), this->jacobian_sparsity.column_indices.data(),
         0, MatrixOrder::COLUMN_MAJOR); // TODO
   }
} // namespace