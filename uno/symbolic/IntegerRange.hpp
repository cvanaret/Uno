// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTEGERRANGE_H
#define UNO_INTEGERRANGE_H

#include <stdexcept>
#include "Collection.hpp"

namespace uno {
   class IntegerRange: public Collection<size_t> {
   public:
      using value_type = size_t;

      explicit IntegerRange(size_t end_value): IntegerRange(0, end_value) { }

      IntegerRange(size_t start_value, size_t end_value): start_value(start_value), end_value(end_value) {
         if (end_value < start_value) {
            throw std::runtime_error("Forward range: end index is smaller than start index");
         }
      }

      // iterable functions
      [[nodiscard]] size_t size() const override {
         return this->end_value - this->start_value;
      }

      [[nodiscard]] size_t dereference_iterator(size_t index) const override {
         return this->start_value + index;
      }

      void increment_iterator(size_t& index) const override {
         ++index;
      }

   protected:
      const size_t start_value;
      const size_t end_value;
   };
} // namespace

#endif // UNO_INTEGERRANGE_H