// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_POINTERWRAPPER_H
#define UNO_POINTERWRAPPER_H

#include <cassert>

namespace uno {
   template <typename T>
   class PointerWrapper {
   public:
      PointerWrapper(T* pointer): pointer(pointer) {
      }

      T* operator->() const {
         return this->pointer;
      }

      T& operator*() const {
         assert(this->pointer != nullptr && "PointerWrapper wrapping NULL, cannot be dereferenced");
         return *this->pointer;
      }

   protected:
      T* pointer;
   };

   template <typename T>
   PointerWrapper<T> wrap_pointer(T* pointer) {
      return {pointer};
   }
} // namespace

#endif // UNO_POINTERWRAPPER_H