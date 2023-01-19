// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGIANGRADIENT_H
#define UNO_LAGRANGIANGRADIENT_H

#include <iostream>
#include <vector>
#include "linear_algebra/Vector.hpp"

// Gradient of the Lagrangian
// Keep the objective and constraint contributions separate. This helps:
// - computing the KKT and FJ stationarity conditions
// - forming quasi-Newton matrices with an objective multiplier
class LagrangianGradient {
public:
   std::vector<double> objective_contribution{};
   std::vector<double> constraints_contribution{};

   explicit LagrangianGradient(size_t number_variables);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] double operator[](size_t i) const;
   //[[nodiscard]] double norm_1() const;
   void resize(size_t new_number_variables);
};

std::ostream& operator<<(std::ostream& stream, const LagrangianGradient& gradient);

#endif // UNO_LAGRANGIANGRADIENT_H


