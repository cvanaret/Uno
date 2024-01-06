// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>

template <typename ElementType>
class VectorExpression {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = ElementType;

   VectorExpression(size_t size, const std::function<ElementType(size_t)>& ith_component);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] ElementType operator[](size_t index) const;

protected:
   const size_t length;
   const std::function<ElementType(size_t)> ith_component;
};

template <typename ElementType>
VectorExpression<ElementType>::VectorExpression(size_t size, const std::function<ElementType(size_t)>& ith_component): length(size), ith_component(ith_component) {
}

template <typename ElementType>
size_t VectorExpression<ElementType>::size() const {
   return this->length;
}

template <typename ElementType>
ElementType VectorExpression<ElementType>::operator[](size_t index) const {
   return this->ith_component(index);
}

#endif // UNO_VECTOREXPRESSION_H