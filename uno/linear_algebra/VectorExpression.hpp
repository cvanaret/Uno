// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>
#include "tools/Collection.hpp"

template <typename ElementType>
class VectorExpression {
public:
   const Collection<size_t>& indices;

   // compatible with algorithms that query the type of the elements
   using value_type = ElementType;

   VectorExpression(const Collection<size_t>& indices, const std::function<ElementType(size_t)>& ith_component);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] ElementType operator[](size_t index) const;

protected:
   const std::function<ElementType(size_t)> ith_component;
};

template <typename ElementType>
VectorExpression<ElementType>::VectorExpression(const Collection<size_t>& indices, const std::function<ElementType(size_t)>& ith_component):
      indices(indices), ith_component(ith_component) {
}

template <typename ElementType>
size_t VectorExpression<ElementType>::size() const {
   return this->indices.size();
}

template <typename ElementType>
ElementType VectorExpression<ElementType>::operator[](size_t index) const {
   return this->ith_component(index);
}

#endif // UNO_VECTOREXPRESSION_H