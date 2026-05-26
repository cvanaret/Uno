// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RANGE_H
#define UNO_RANGE_H

#include <cstddef>
#include <stdexcept>

namespace uno {
   // direction of the range (FORWARD = increasing or BACKWARD = decreasing)
   enum RangeDirection {
      FORWARD, BACKWARD
   };

   // Default direction is FORWARD (increasing)
   template <RangeDirection direction = FORWARD>
   class Range {
   public:
      using value_type = size_t;

      class iterator {
      public:
         using value_type = size_t;

         iterator(size_t start): index(start) { }

         [[nodiscard]] constexpr value_type operator*() const {
            return this->index;
         }

         constexpr iterator& operator++() {
            if constexpr (direction == FORWARD) {
               ++this->index;
            }
            else {
               --this->index;
            }
            return *this;
         }

         friend constexpr bool operator!=(const iterator& a, const iterator& b) {
            return a.index != b.index;
         }

      protected:
         size_t index;
      };

      explicit Range(size_t end_value): Range(0, end_value) {
         static_assert(direction == FORWARD);
      }

      Range(size_t start_value, size_t end_value): start_value(start_value), end_value(end_value) {
         if (direction == FORWARD && end_value < start_value) {
            throw std::invalid_argument("Forward range: end index is smaller than start index");
         }
         else if (direction == BACKWARD && end_value > start_value) {
            throw std::invalid_argument("Backward range: end index is larger than start index");
         }
      }

      [[nodiscard]] constexpr size_t size() const noexcept {
         if constexpr (direction == FORWARD) {
            return this->end_value - this->start_value;
         }
         else {
            return this->start_value - this->end_value;
         }
      }

      constexpr iterator begin() const noexcept {
         return iterator(this->start_value);
      }

      constexpr iterator end() const noexcept {
         return iterator(this->end_value);
      }

   protected:
      const size_t start_value;
      const size_t end_value;
   };

   using ForwardRange = Range<FORWARD>;
   using BackwardRange = Range<BACKWARD>;
} // namespace

#endif // UNO_RANGE_H