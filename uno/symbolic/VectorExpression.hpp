// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include <functional>
#include "Collection.hpp"

template <typename Indices, typename Callable>
class VectorExpression {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = double;

   VectorExpression(Indices&& indices, Callable&& component_function);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] double operator[](size_t index) const;

   void for_each(const std::function<void (size_t, size_t)>& f) const;

protected:
   Indices indices; // store const reference or rvalue (temporary)
   Callable ith_component;
};

template <typename Indices, typename Callable>
VectorExpression<Indices, Callable>::VectorExpression(Indices&& indices, Callable&& component_function):
      indices(std::forward<Indices>(indices)), ith_component(std::forward<Callable>(component_function)) {
}

template <typename Indices, typename Callable>
size_t VectorExpression<Indices, Callable>::size() const {
   return this->indices.size();
}

template <typename Indices, typename Callable>
double VectorExpression<Indices, Callable>::operator[](size_t index) const {
   return this->ith_component(index);
}

template <typename Indices, typename Callable>
void VectorExpression<Indices, Callable>::for_each(const std::function<void (size_t, size_t)>& f) const {
   this->indices.for_each(f);
}

#endif // UNO_VECTOREXPRESSION_H