// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>
#include "tools/Collection.hpp"

template <typename ElementType, typename Indices>
class VectorExpression {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = ElementType;

   VectorExpression(Indices&& indices, const std::function<ElementType(size_t)>& ith_component);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] ElementType operator[](size_t index) const;

   void for_each(const std::function<void (size_t, size_t)>& f) const;

protected:
   Indices indices; // store const reference or rvalue (temporary)
   const std::function<ElementType(size_t)> ith_component;
};

template <typename ElementType, typename Indices>
VectorExpression<ElementType, Indices>::VectorExpression(Indices&& indices, const std::function<ElementType(size_t)>& ith_component):
      indices(std::forward<Indices>(indices)), ith_component(ith_component) {
}

template <typename ElementType, typename Indices>
size_t VectorExpression<ElementType, Indices>::size() const {
   return this->indices.size();
}

template <typename ElementType, typename Indices>
ElementType VectorExpression<ElementType, Indices>::operator[](size_t index) const {
   return this->ith_component(index);
}

template <typename ElementType, typename Indices>
void VectorExpression<ElementType, Indices>::for_each(const std::function<void (size_t, size_t)>& f) const {
   this->indices.for_each(f);
}

#endif // UNO_VECTOREXPRESSION_H