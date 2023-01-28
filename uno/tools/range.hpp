// Copyright (c) 2018-2023 Charlie Vanaret
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
class range {
public:
   explicit range(size_t end_index);
   range(size_t start_index, size_t end_index);

   // iterable functions
   [[nodiscard]] const range& begin() const;
   [[nodiscard]] const range& end() const;

   // iterator functions
   [[nodiscard]] bool operator!=(const range&) const;
   void operator++();
   [[nodiscard]] size_t operator*() const;

   // size
   [[nodiscard]] size_t size() const;

protected:
   size_t end_index{};
   size_t current_index{};
};

template <RangeDirection direction>
inline range<direction>::range(size_t end_index): end_index(end_index), current_index(0) {
   static_assert(direction == FORWARD);
}

template <RangeDirection direction>
inline range<direction>::range(size_t start_index, size_t end_index): end_index(end_index), current_index(start_index) {
}

template <RangeDirection direction>
inline const range<direction>& range<direction>::begin() const {
   return *this;
}

template <RangeDirection direction>
inline const range<direction>& range<direction>::end() const {
   return *this;
}

template <RangeDirection direction>
inline bool range<direction>::operator!=(const range<direction>&) const {
   if constexpr (direction == FORWARD) {
      return (this->current_index < this->end_index);
   }
   else {
      return (this->current_index > this->end_index);
   }
}

template <RangeDirection direction>
inline void range<direction>::operator++() {
   if constexpr (direction == FORWARD) {
      this->current_index++;
   }
   else {
      this->current_index--;
   }
}

template <RangeDirection direction>
inline size_t range<direction>::operator*() const {
   return this->current_index;
}

template <RangeDirection direction>
inline size_t range<direction>::size() const {
   if constexpr (direction == FORWARD) {
      return this->end_index - this->current_index;
   }
   else {
      return this->current_index - this->end_index;
   }
}

#endif // UNO_RANGE_H