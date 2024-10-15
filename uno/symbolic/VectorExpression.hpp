// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>
#include "Collection.hpp"

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
            this->index++;
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

      VectorExpression(const Indices& indices, Callable&& component_function);
      [[nodiscard]] size_t size() const { return this->indices.size(); }
      [[nodiscard]] double operator[](size_t index) const;

      iterator begin() const { return iterator(*this, 0); }
      iterator end() const { return iterator(*this, this->size()); }

   protected:
      const Indices& indices; // store const reference or rvalue (temporary)
      const Callable component_function;
   };

   template <typename Indices, typename Callable>
   VectorExpression<Indices, Callable>::VectorExpression(const Indices& indices, Callable&& component_function):
         indices(indices), component_function(std::forward<Callable>(component_function)) {
   }

   template <typename Indices, typename Callable>
   double VectorExpression<Indices, Callable>::operator[](size_t index) const {
      return this->component_function(index);
   }
} // namespace

#endif // UNO_VECTOREXPRESSION_H
