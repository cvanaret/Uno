// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INERTIA_H
#define UNO_INERTIA_H

#include <functional>
#include <ostream>

namespace uno {
   class Inertia {
   public:
      const size_t positive;
      const size_t negative;
      const size_t zero;

      Inertia(size_t positive, size_t negative, size_t zero):
         positive(positive), negative(negative), zero(zero) { }

      [[nodiscard]] bool operator==(const Inertia& other) const {
         return this->positive == other.positive && this->negative == other.negative && this->zero == other.zero;
      }

      friend std::ostream& operator<<(std::ostream& stream, const Inertia& inertia);
   };

   // Evaluates the inertia of the full matrix after each numerical factorization.
   using InertiaEvaluator = std::function<Inertia()>;

   inline std::ostream& operator<<(std::ostream& stream, const Inertia& inertia) {
      stream << '(' << inertia.positive << ", " << inertia.negative << ", " << inertia.zero << ')';
      return stream;
   }
} // namespace

#endif // UNO_INERTIA_H