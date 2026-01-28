// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGIANGRADIENT_H
#define UNO_LAGRANGIANGRADIENT_H

#include "linear_algebra/Vector.hpp"

namespace uno {
   // Gradient of the Lagrangian
   // Keep the objective and constraint contributions separate. This helps:
   // - computing the KKT and FJ stationarity conditions
   // - forming quasi-Newton matrices with an objective multiplier
   class LagrangianGradient {
   public:
      Vector<double> objective_contribution{};
      Vector<double> constraints_contribution{};

      using value_type = double;

      explicit LagrangianGradient(size_t number_variables);
      [[nodiscard]] size_t size() const;
      [[nodiscard]] double operator[](size_t variable_index) const;
      void resize(size_t number_variables);
   };

   std::ostream& operator<<(std::ostream& stream, const LagrangianGradient& gradient);
} // namespace

#endif // UNO_LAGRANGIANGRADIENT_H