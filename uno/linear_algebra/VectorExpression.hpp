// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>

template <typename T>
class VectorExpression {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = T;

   VectorExpression(size_t size, const std::function<T(size_t)>& ith_component);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] T operator[](size_t i) const;

protected:
   const size_t length;
   const std::function<T (size_t)> ith_component;
};

template <typename T>
VectorExpression<T>::VectorExpression(size_t size, const std::function<T(size_t)>& ith_component): length(size), ith_component(ith_component) {
}

template <typename T>
size_t VectorExpression<T>::size() const {
   return this->length;
}

template <typename T>
T VectorExpression<T>::operator[](size_t i) const {
   return this->ith_component(i);
}

#endif // UNO_VECTOREXPRESSION_H