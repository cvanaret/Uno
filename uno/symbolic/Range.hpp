// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RANGE_H
#define UNO_RANGE_H

#include <stdexcept>
#include "Collection.hpp"

// direction of the range (FORWARD = increasing or BACKWARD = decreasing)
enum RangeDirection {
   FORWARD, BACKWARD
};

// Default direction is FORWARD (increasing)
template <RangeDirection direction = FORWARD>
class Range: public Collection<size_t> {
public:
   // https://internalpointers.com/post/writing-custom-iterators-modern-cpp
   class iterator {
   public:
      explicit iterator(size_t value) : value(value) { }
      const size_t& operator*() const { return this->value; }
      const size_t& operator->() { return this->value; }
      // prefix increment
      iterator& operator++() {
         if constexpr (direction == FORWARD) {
            this->value++;
         }
         else {
            this->value--;
         }
         return *this;
      }

      friend bool operator==(const iterator& a, const iterator& b) { return a.value == b.value; };
      friend bool operator!=(const iterator& a, const iterator& b) { return a.value != b.value; };

   private:
      size_t value;
   };

   using value_type = size_t;

   explicit Range(size_t end_value);
   Range(size_t start_value, size_t end_value);

   // iterable functions
   [[nodiscard]] iterator begin() const { return iterator(this->start_value); }
   [[nodiscard]] iterator end() const { return iterator(this->end_value); }
   [[nodiscard]] size_t size() const override;

   [[nodiscard]] size_t dereference_iterator(size_t index) const override;
   void increment_iterator(size_t& index) const override;

protected:
   const size_t start_value;
   const size_t end_value;
};

template <RangeDirection direction>
inline Range<direction>::Range(size_t end_value): Range(0, end_value) {
   static_assert(direction == FORWARD);
}

template <RangeDirection direction>
inline Range<direction>::Range(size_t start_value, size_t end_value): Collection<size_t>(), start_value(start_value), end_value(end_value) {
   if (direction == FORWARD && end_value < start_value) {
      throw std::runtime_error("Forward range: end index is smaller than start index\n");
   }
   else if (direction == BACKWARD && end_value > start_value) {
      throw std::runtime_error("Backward range: end index is larger than start index\n");
   }
}

template <RangeDirection direction>
inline size_t Range<direction>::size() const {
   if constexpr (direction == FORWARD) {
      return this->end_value - this->start_value;
   }
   else {
      return this->start_value - this->end_value;
   }
}

template <RangeDirection direction>
size_t Range<direction>::dereference_iterator(size_t index) const {
   if constexpr (direction == FORWARD) {
      return this->start_value + index;
   }
   else {
      return this->start_value - index;
   }
}

template <RangeDirection direction>
void Range<direction>::increment_iterator(size_t& index) const {
   index++;
}

using ForwardRange = Range<FORWARD>;
using BackwardRange = Range<BACKWARD>;

#endif // UNO_RANGE_H
