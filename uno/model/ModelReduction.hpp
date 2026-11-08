// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODELREDUCTION_H
#define UNO_MODELREDUCTION_H

#include <cstddef>
#include <vector>

namespace uno {
   // forward declaration
   class Model;

   // elimination of the fixed variables (equal finite bounds) of a model. The constraints that only depend on fixed
   // variables (according to the Jacobian sparsity) are constant: they are removed if satisfied, otherwise they are kept
   // and the solver reports the infeasibility. The reduction only uses the bounds and the Jacobian sparsity of the model,
   // and is therefore independent of the modeling language
   class ModelReduction {
   public:
      // feasibility_tolerance: largest violation of a removed constraint
      [[nodiscard]] static ModelReduction detect(const Model& model, double feasibility_tolerance);

      std::vector<size_t> fixed_variables{};
      std::vector<size_t> free_variables{}; // reduced index -> original index
      std::vector<size_t> kept_constraints{}; // reduced index -> original index
      std::vector<size_t> removed_constraints{}; // satisfied constraints of fixed variables only

      [[nodiscard]] bool is_trivial() const;
   };
} // namespace

#endif // UNO_MODELREDUCTION_H
