// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLASVECTOR_H
#define UNO_BLASVECTOR_H

#include <cstddef>
#include "BLAS.hpp"

namespace uno {
   // constant contiguous array in memory on which BLAS can be called
   template <typename T>
   class BLASVector {
   public:
      using value_type = T;

      BLASVector() = default;
      virtual ~BLASVector() = default;

      [[nodiscard]] virtual size_t size() const = 0;
      [[nodiscard]] bool empty() const {
         return (this->size() == 0);
      }
      [[nodiscard]] virtual const T* data() const = 0;
      [[nodiscard]] virtual const T& operator[](size_t index) const = 0;
   };

   // non-constant contiguous array in memory on which BLAS can be called
   template <typename T>
   class MutableBLASVector: public BLASVector<T> {
   public:
      MutableBLASVector() = default;

      template <typename Vector>
      MutableBLASVector& operator=(const Vector& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr int increment = 1;
         BLAS_copy_vector(&size, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      template <typename Vector>
      MutableBLASVector& operator+=(const Vector& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr double factor = 1.;
         constexpr int increment = 1;
         BLAS_add_vectors(&size, &factor, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      template <typename Vector>
      MutableBLASVector& operator-=(const Vector& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr double factor = -1.;
         constexpr int increment = 1;
         BLAS_add_vectors(&size, &factor, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      [[nodiscard]] virtual T* data() = 0;
      [[nodiscard]] virtual T& operator[](size_t index) = 0;

      void scale(T factor) {
         const int size = static_cast<int>(this->size());
         constexpr int increment = 1;
         BLAS_scale_vector(&size, &factor, this->data(), &increment);
      }

      void operator*=(T factor) {
         this->scale(factor);
      }

      BLASVector<T>& operator+=(const BLASVector<T>& other) {
         const int size = static_cast<int>(this->size());
         constexpr double factor = 1.;
         constexpr int increment = 1;
         // this = this + factor * other
         BLAS_add_vectors(&size, &factor, other.data(), &increment, this->data(), &increment);
         return *this;
      }
   };

   inline double dot(const BLASVector<double>& x, const BLASVector<double>& y) {
      const int size = static_cast<int>(std::min(x.size(), y.size()));
      constexpr int increment = 1;
      return BLAS_dot_product(&size, x.data(), &increment, y.data(), &increment);
   }
} // namespace

#endif //UNO_BLASVECTOR_H