// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXPRESSION_H
#define UNO_EXPRESSION_H

#include <functional>

class Expression {
public:
   Expression() = default;
   virtual ~Expression() = default;

   virtual void for_each(const std::function<void (size_t row_index, size_t column_index, double element)>& f) const = 0;
};

#endif // UNO_EXPRESSION_H