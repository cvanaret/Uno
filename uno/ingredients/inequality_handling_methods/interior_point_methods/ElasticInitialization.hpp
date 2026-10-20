// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ELASTICINITIALIZATION_H
#define UNO_ELASTICINITIALIZATION_H

#include <cmath>
#include <utility>

namespace uno {
   // unique positive solution (a, b) of the system
   //    a - b = r
   //    1/a + 1/b = 2/epsilon      (epsilon > 0)
   // that is a = (epsilon + r + sqrt(r^2 + epsilon^2))/2 and b = (epsilon - r + sqrt(r^2 + epsilon^2))/2.
   // The smaller root is evaluated without cancellation: sqrt(r^2 + epsilon^2) - |r| = epsilon^2/(sqrt(r^2 + epsilon^2) + |r|).
   // Both values are bounded below by epsilon/2.
   // This is the closed form of the initialization of a pair of nonnegative variables (a, b) with duals (mu/a, mu/b)
   // linked by a linear equation and whose duals sum to 2 mu/epsilon (IPOPT's restoration phase, eq. (33))
   [[nodiscard]] inline std::pair<double, double> solve_elastic_pair(double r, double epsilon) {
      const double radical = std::hypot(r, epsilon);
      const double absolute_r = std::abs(r);
      const double large_root = (epsilon + absolute_r + radical) / 2.;
      const double small_root = (epsilon + epsilon * epsilon / (radical + absolute_r)) / 2.;
      return (0. <= r) ? std::pair{large_root, small_root} : std::pair{small_root, large_root};
   }
} // namespace

#endif // UNO_ELASTICINITIALIZATION_H
