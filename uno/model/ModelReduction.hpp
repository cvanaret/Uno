// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODELREDUCTION_H
#define UNO_MODELREDUCTION_H

#include <cstddef>
#include <optional>
#include <vector>

namespace uno {
   // forward declarations
   class Model;
   class Options;

   // linear constraint cl <= a x_i + b <= cu that depends on a single non-fixed variable x_i: it is a bound on x_i
   struct BoundConstraint {
      size_t constraint_index;
      size_t variable_index;
      double coefficient; // a
   };

   // structural reduction of a model. It is computed from the variable bounds, the linearity of the constraints and the
   // Jacobian (sparsity pattern and values), and is therefore independent of the modeling language:
   // - nonlinear constraints whose Jacobian is constant at a few random points are considered linear (optional);
   // - linear constraints with a single non-fixed variable are turned into bounds on that variable, and satisfied linear
   //   constraints without non-fixed variables are removed (optional);
   // - fixed variables (equal bounds, possibly after bound tightening) are eliminated. Since a fixed variable may turn
   //   other constraints into bound constraints, the last two steps are repeated in rounds until no variable gets fixed
   class ModelReduction {
   public:
      [[nodiscard]] static ModelReduction detect(const Model& model, const Options& options);

      std::vector<size_t> fixed_variables{}; // in the order of elimination
      std::vector<size_t> elimination_rounds{}; // round at which each fixed variable was fixed (0: fixed in the model)
      std::vector<size_t> free_variables{}; // reduced index -> original index
      std::vector<size_t> kept_constraints{}; // reduced index -> original index
      std::vector<BoundConstraint> bound_constraints{};
      std::vector<size_t> redundant_constraints{}; // satisfied linear constraints without non-fixed variables
      std::vector<size_t> detected_linear_constraints{}; // declared nonlinear by the model, but most likely linear
      // variable bounds (original space), tightened by the bound constraints. The fixed values are the lower bounds
      std::vector<double> variables_lower_bounds{};
      std::vector<double> variables_upper_bounds{};
      // index (into bound_constraints) of the bound constraint that defines each variable bound, if any
      std::vector<std::optional<size_t>> lower_bound_origins{};
      std::vector<std::optional<size_t>> upper_bound_origins{};

      // linearity of the constraints of the model (original space), including the detected linear constraints
      [[nodiscard]] std::vector<bool> linear_constraints(const Model& model) const;
      [[nodiscard]] bool is_trivial() const;
   };
} // namespace

#endif // UNO_MODELREDUCTION_H
