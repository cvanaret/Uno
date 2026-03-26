#include "EvaluationCache.hpp"
#include "linear_algebra/Indexing.hpp"
#include "model/Model.hpp"

namespace uno {
   EvaluationCache::EvaluationCache(const Model& model):
         number_jacobian_nonzeros(model.number_jacobian_nonzeros()),
         jacobian_sparsity(this->number_jacobian_nonzeros),
         // pass (a const reference of) the Jacobian sparsity to both evaluations
         current_evaluations(model, &this->jacobian_sparsity),
         trial_evaluations(model, &this->jacobian_sparsity) {
      model.compute_jacobian_sparsity(this->jacobian_sparsity.row_indices.data(), this->jacobian_sparsity.column_indices.data(),
         Indexing::C_indexing, MatrixOrder::ROW_MAJOR); // TODO
   }
} // namespace