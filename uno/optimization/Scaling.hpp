// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALING_H
#define UNO_SCALING_H

#include <cstddef>
#include <vector>

namespace uno {
   // forward declarations
   class EvaluationSpace;
   class Iterate;

   class Scaling {
   public:
      Scaling(const Iterate& initial_iterate, const EvaluationSpace& evaluation_space, double gradient_threshold);
      ~Scaling() = default;

      [[nodiscard]] double get_objective_scaling() const;
      [[nodiscard]] double get_constraint_scaling(size_t constraint_index) const;

   protected:
      const double gradient_threshold;
      double objective_scaling;
      std::vector<double> constraint_scaling;
   };
} // namespace

#endif // UNO_SCALING_H