// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSEVECTOR_WRAPPER_H
#define UNO_SPARSEVECTOR_WRAPPER_H

#include "linear_algebra/SparseVector.hpp"

#ifdef __cplusplus
namespace uno {
   extern "C" {
#endif

   // T = double
   SparseVector<double>* SparseVector_new(size_t capacity);
   void SparseVector_delete(SparseVector<double>* vector);
   void SparseVector_insert(SparseVector<double>* vector, size_t key, double value);
   void SparseVector_clear(SparseVector<double>* vector);
   void SparseVector_display(SparseVector<double>* vector);

#ifdef __cplusplus
   }
} // namespace
#endif

#endif // UNO_SPARSEVECTOR_WRAPPER_H
