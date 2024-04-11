// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RANGE_H
#define UNO_RANGE_H

#include "Collection.hpp"
#include <iostream>

// direction of the range (FORWARD = increasing or BACKWARD = decreasing)
enum RangeDirection {
   FORWARD, BACKWARD
};

// https://en.wikipedia.org/wiki/Generator_(computer_programming)#C++
// Default direction is FORWARD (increasing)
template <RangeDirection direction = FORWARD>
class Range: public Collection<size_t> {
public:
   // https://internalpointers.com/post/writing-custom-iterators-modern-cpp
   class Iterator {
   public:
      explicit Iterator(size_t value) : value(value) {
      }
      const size_t& operator*() const { return this->value; }
      const size_t& operator->() { return this->value; }
      // prefix increment
      Iterator& operator++() {
         if constexpr (direction == FORWARD) {
            this->value++;
         }
         else {
            this->value--;
         }
         return *this;
      }

      friend bool operator== (const Iterator& a, const Iterator& b) { return a.value == b.value; };
      friend bool operator!= (const Iterator& a, const Iterator& b) { return a.value != b.value; };

   private:
      size_t value;
   };

   using value_type = size_t;

   explicit Range(size_t end_index);
   Range(size_t start_index, size_t end_index);

   // iterable functions
   [[nodiscard]] Iterator begin() const { return Iterator(this->start_index); }
   [[nodiscard]] Iterator end() const { return Iterator(this->end_index); }
   [[nodiscard]] size_t size() const override;

   void for_each(const std::function<void (size_t, size_t)>& f) const override;

protected:
   const size_t start_index;
   const size_t end_index;
};

template <RangeDirection direction>
inline Range<direction>::Range(size_t end_index): Range(0, end_index) {
   static_assert(direction == FORWARD);
}

template <RangeDirection direction>
inline Range<direction>::Range(size_t start_index, size_t end_index): start_index(start_index), end_index(end_index) {
}

template <RangeDirection direction>
inline size_t Range<direction>::size() const {
   if constexpr (direction == FORWARD) {
      return this->end_index - this->start_index;
   }
   else {
      return this->start_index - this->end_index;
   }
}

template <RangeDirection direction>
inline void Range<direction>::for_each(const std::function<void (size_t, size_t)>& f) const {
   size_t index_no_offset = 0;
   for (size_t index: *this) {
      f(index_no_offset, index);
      index_no_offset++;
   }
}

#endif // UNO_RANGE_H