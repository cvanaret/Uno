// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALING_H
#define UNO_SCALING_H

#include <vector>
#include <cstddef>

namespace uno {
   // forward declarations
   class COOSparsity;
   class Model;
   template <typename ElementType>
   class Vector;

   class Scaling {
   public:
      Scaling(size_t number_constraints, double gradient_threshold);

      void compute(const Model& model, const Vector<double>& objective_gradient, const Vector<double>& jacobian_values);

      [[nodiscard]] bool is_objective_scaled() const;
      [[nodiscard]] bool are_constraints_scaled() const;

      [[nodiscard]] double get_objective_scaling() const;
      [[nodiscard]] const std::vector<double>& get_constraint_scaling() const;

   protected:
      const double gradient_threshold;
      double objective_scaling;
      std::vector<double> constraint_scaling;
      // lazy flags
      bool is_objective_scaled_flag{false};
      bool are_constraints_scaled_flag{false};
   };
} // namespace

#endif // UNO_SCALING_H