// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RANGE_H
#define UNO_RANGE_H

// direction of the range (FORWARD = increasing or BACKWARD = decreasing)
enum RangeDirection {
   FORWARD, BACKWARD
};

// https://en.wikipedia.org/wiki/Generator_(computer_programming)#C++
// Default direction is FORWARD (increasing)
template <RangeDirection direction = FORWARD>
class Range {
public:
   explicit Range(size_t end_index);
   Range(size_t start_index, size_t end_index);

   // iterable functions
   [[nodiscard]] const Range& begin() const;
   [[nodiscard]] const Range& end() const;

   // iterator functions
   [[nodiscard]] bool operator!=(const Range&) const;
   void operator++();
   [[nodiscard]] size_t operator*() const;

   // size
   [[nodiscard]] size_t size() const;

protected:
   size_t end_index{};
   size_t current_index{};
};

template <RangeDirection direction>
inline Range<direction>::Range(size_t end_index): end_index(end_index), current_index(0) {
   static_assert(direction == FORWARD);
}

template <RangeDirection direction>
inline Range<direction>::Range(size_t start_index, size_t end_index): end_index(end_index), current_index(start_index) {
}

template <RangeDirection direction>
inline const Range<direction>& Range<direction>::begin() const {
   return *this;
}

template <RangeDirection direction>
inline const Range<direction>& Range<direction>::end() const {
   return *this;
}

template <RangeDirection direction>
inline bool Range<direction>::operator!=(const Range<direction>&) const {
   if constexpr (direction == FORWARD) {
      return (this->current_index < this->end_index);
   }
   else {
      return (this->current_index > this->end_index);
   }
}

template <RangeDirection direction>
inline void Range<direction>::operator++() {
   if constexpr (direction == FORWARD) {
      this->current_index++;
   }
   else {
      this->current_index--;
   }
}

template <RangeDirection direction>
inline size_t Range<direction>::operator*() const {
   return this->current_index;
}

template <RangeDirection direction>
inline size_t Range<direction>::size() const {
   if constexpr (direction == FORWARD) {
      return this->end_index - this->current_index;
   }
   else {
      return this->current_index - this->end_index;
   }
}

#endif // UNO_RANGE_H