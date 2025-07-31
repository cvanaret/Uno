// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGIANGRADIENT_H
#define UNO_LAGRANGIANGRADIENT_H

#include <ostream>
#include "linear_algebra/Vector.hpp"

namespace uno {
   // Gradient of the Lagrangian
   // Keep the objective and constraint contributions separate. This helps:
   // - computing the KKT and FJ stationarity conditions
   // - forming quasi-Newton matrices with an objective multiplier
   template <typename ElementType>
   class LagrangianGradient {
   public:
      Vector<ElementType> objective_contribution{};
      Vector<ElementType> constraints_contribution{};

      using value_type = ElementType;

      explicit LagrangianGradient(size_t number_variables);
      [[nodiscard]] size_t size() const;
      [[nodiscard]] ElementType operator[](size_t variable_index) const;
      void resize(size_t number_variables);
   };

   template <typename ElementType>
   LagrangianGradient<ElementType>::LagrangianGradient(size_t number_variables) :
         objective_contribution(number_variables),
         constraints_contribution(number_variables) {
   }

   template <typename ElementType>
   size_t LagrangianGradient<ElementType>::size() const {
      return this->objective_contribution.size();
   }

   // access i-th element
   template <typename ElementType>
   ElementType LagrangianGradient<ElementType>::operator[](size_t variable_index) const {
      return this->objective_contribution[variable_index] + this->constraints_contribution[variable_index];
   }

   template <typename ElementType>
   void LagrangianGradient<ElementType>::resize(size_t number_variables) {
      this->objective_contribution.resize(number_variables);
      this->constraints_contribution.resize(number_variables);
   }

   template <typename ElementType>
   std::ostream& operator<<(std::ostream& stream, const LagrangianGradient<ElementType>& gradient) {
      for (size_t variable_index: Range(gradient.constraints_contribution.size())) {
         stream << gradient[variable_index] << ' ';
      }
      stream << '\n';
      return stream;
   }
} // namespace

#endif // UNO_LAGRANGIANGRADIENT_H