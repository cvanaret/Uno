// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGIANGRADIENT_H
#define UNO_LAGRANGIANGRADIENT_H

#include <ostream>
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class LagrangianGradient;

   // The LagrangianGradient class represents the dense Lagrangian gradient broken into:
   // - the objective contribution \nabla f(x_k)
   // - the constraint contribution -\nabla c(x_k) y_k - z_k
   // The two contributions can be assembled with a given objective multiplier \rho (see assemble(double)):
   // \rho \nabla f(x_k) - \nabla c(x_k) y_k - z_k
   // The resulting object, an AssembledLagrangianGradient, is a wrapper around the LagrangianGradient and \rho

   template <typename ElementType>
   class AssembledLagrangianGradient {
   public:
      AssembledLagrangianGradient(const LagrangianGradient<ElementType>& lagrangian_gradient, ElementType objective_multiplier);

      [[nodiscard]] size_t size() const;
      [[nodiscard]] ElementType operator[](size_t variable_index) const;

   protected:
      const LagrangianGradient<ElementType>& lagrangian_gradient;
      const ElementType objective_multiplier;
   };

   template<typename ElementType>
   AssembledLagrangianGradient<ElementType>::AssembledLagrangianGradient(const LagrangianGradient<ElementType> &lagrangian_gradient,
      ElementType objective_multiplier): lagrangian_gradient(lagrangian_gradient), objective_multiplier(objective_multiplier) { }

   template <typename ElementType>
   size_t AssembledLagrangianGradient<ElementType>::size() const {
      return this->lagrangian_gradient.size();
   }

  // access i-th element
  template <typename ElementType>
  ElementType AssembledLagrangianGradient<ElementType>::operator[](size_t variable_index) const {
     return this->objective_multiplier * this->lagrangian_gradient.objective_contribution[variable_index] +
        this->lagrangian_gradient.constraints_contribution[variable_index];
  }

   template <typename ElementType>
   std::ostream& operator<<(std::ostream& stream, const AssembledLagrangianGradient<ElementType>& gradient) {
      for (size_t variable_index: Range(gradient.size())) {
         stream << gradient[variable_index] << ' ';
      }
      stream << '\n';
      return stream;
   }

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
      void resize(size_t number_variables);
      AssembledLagrangianGradient<ElementType> assemble(double objective_multiplier) const;
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

   template <typename ElementType>
   void LagrangianGradient<ElementType>::resize(size_t number_variables) {
      this->objective_contribution.resize(number_variables);
      this->constraints_contribution.resize(number_variables);
   }

   template <typename ElementType>
   AssembledLagrangianGradient<ElementType> LagrangianGradient<ElementType>::assemble(double objective_multiplier) const {
      // TODO use existing Expression
      return {*this, objective_multiplier};
   }
} // namespace

#endif // UNO_LAGRANGIANGRADIENT_H