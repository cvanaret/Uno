// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include "symbolic_traits.hpp"

namespace uno {
   template <typename Indices, typename Callable>
   class VectorExpression {
   public:
      class iterator {
      public:
         using value_type = std::pair<size_t, double>;

         iterator(const VectorExpression& expression, size_t index): expression(expression), index(index) { }

         [[nodiscard]] std::pair<size_t, double> operator*() const {
            const auto [_, expression_index] = expression.indices.dereference_iterator(this->index);
            return {expression_index, expression.component_function[index]};
         }

         iterator& operator++() {
            ++this->index;
            return *this;
         }

         friend bool operator!=(const iterator& a, const iterator& b) {
            return &a.expression != &b.expression || a.index != b.index;
         }

      protected:
         const VectorExpression& expression;
         size_t index;
      };

      // compatible with algorithms that query the type of the elements
      using value_type = double;

      VectorExpression(const Indices& indices, Callable&& component_function):
         indices(indices), component_function(std::forward<Callable>(component_function)) { }

      [[nodiscard]] constexpr size_t size() const noexcept {
         return this->indices.size();
      }

      [[nodiscard]] constexpr double operator[](size_t index) const noexcept {
         return this->component_function(index);
      }

      constexpr iterator begin() const noexcept {
         return iterator(*this, 0);
      }

      constexpr iterator end() const noexcept {
         return iterator(*this, this->size());
      }

   protected:
      const Indices& indices; // store const reference or rvalue (temporary)
      storage_t<Callable> component_function;
   };
} // namespace

#endif // UNO_VECTOREXPRESSION_H