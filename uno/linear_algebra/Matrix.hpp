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
};

#endif // UNO_MATRIX_H