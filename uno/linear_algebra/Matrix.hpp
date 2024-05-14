// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIX_H
#define UNO_MATRIX_H

template <typename IndexType, typename ElementType>
class Matrix {
public:
   using value_type = ElementType;

   Matrix() = default;
   virtual ~Matrix() = default;

   virtual void for_each(const std::function<void(IndexType, IndexType, ElementType)>& f) const = 0;
};

#endif // UNO_MATRIX_H